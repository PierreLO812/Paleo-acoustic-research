"""M1 — Audio-only acoustic measurement module.

Reads only WAV waveforms / numpy arrays (never accesses cavity filter targets or
model priors) to measure:
  1. Fundamental frequency f0 via YIN (20–400 Hz).
  2. Formant candidates via LPC root-solving with harmonic-locking checks.
  3. Primary cavity transfer function H(f) = P_sf(f) / P_ss(f) and coherence gamma^2(f).
  4. Contiguous spectral band energy fractions (summing to 1.0) and harmonic band power.
  5. Glottal source harmonic spectral tilt (dB/octave in 100–2000 Hz).
"""

from __future__ import annotations

import argparse
import json
from math import gcd
from pathlib import Path
from typing import Iterable

import numpy as np
from scipy import linalg, signal
from scipy.io import wavfile


def read_wav_float(path: Path) -> tuple[int, np.ndarray, bool]:
    """Read WAV file as float64 mono signal in [-1, 1] and detect digital clipping."""
    if not path.exists():
        raise FileNotFoundError(f"Missing audio file: {path}")
    sr, data = wavfile.read(path)
    x = np.asarray(data, dtype=np.float64)
    if x.ndim > 1:
        x = np.mean(x, axis=1)
    if np.issubdtype(data.dtype, np.integer):
        info = np.iinfo(data.dtype)
        scale = float(max(abs(info.min), info.max))
        x = x / scale
    clipped = bool(np.any(np.abs(x) >= 0.999))
    return int(sr), x, clipped


def measure_f0(
    x: np.ndarray,
    fs: int,
    *,
    fmin: float = 20.0,
    fmax: float = 400.0,
    frame_ms: float = 60.0,
    hop_ms: float = 15.0,
    yin_threshold: float = 0.15,
) -> dict:
    """Estimate f0 over voiced frames using the YIN algorithm (de Cheveigné & Kawahara 2002)."""
    x = np.asarray(x, dtype=np.float64)
    if x.ndim > 1:
        x = np.mean(x, axis=1)
    if x.size == 0:
        return {"measured_hz": 0.0, "measured_std_hz": 0.0, "method": "yin", "voiced_frames": 0}

    tau_min = max(2, int(np.floor(fs / fmax)))
    tau_max = max(tau_min + 2, int(np.ceil(fs / fmin)))
    win_len = max(tau_max * 2 + 4, int(round(frame_ms * 1e-3 * fs)))
    hop_len = max(1, int(round(hop_ms * 1e-3 * fs)))

    if len(x) < win_len:
        win_len = len(x)
        tau_max = min(tau_max, max(tau_min + 1, win_len // 2 - 1))

    rms_global = float(np.sqrt(np.mean(x**2)))
    rms_thresh = max(1e-5, 0.05 * rms_global)

    f0_estimates: list[float] = []
    # Integrate-length W for difference function
    w_len = win_len - tau_max
    if w_len < 16:
        return {"measured_hz": 0.0, "measured_std_hz": 0.0, "method": "yin", "voiced_frames": 0}

    for start in range(0, len(x) - win_len + 1, hop_len):
        frame = x[start : start + win_len]
        frame = frame - np.mean(frame)
        if float(np.sqrt(np.mean(frame**2))) < rms_thresh:
            continue

        # Fast vectorized YIN difference function via FFT cross-correlation + cumulative energy
        seg = frame[:w_len]
        energy0 = float(np.dot(seg, seg))
        csum = np.concatenate(([0.0], np.cumsum(frame**2)))
        energy_tau = csum[w_len : w_len + tau_max + 1] - csum[: tau_max + 1]

        n_fft = 1 << int(np.ceil(np.log2(2 * win_len)))
        r_fft = np.fft.irfft(
            np.fft.rfft(frame, n=n_fft) * np.conj(np.fft.rfft(seg, n=n_fft)),
            n=n_fft,
        )
        cross = r_fft[: tau_max + 1]
        d = np.maximum(0.0, energy0 + energy_tau - 2.0 * cross)

        cmndf = np.ones(tau_max + 1, dtype=np.float64)
        running = np.cumsum(d[1:])
        denom = np.where(running > 1e-18, running, 1e-18)
        cmndf[1:] = d[1:] * np.arange(1, tau_max + 1) / denom

        search = cmndf[tau_min : tau_max + 1]
        below = np.where(search < yin_threshold)[0]
        if below.size > 0:
            idx = int(below[0])
            while idx + 1 < len(search) and search[idx + 1] < search[idx]:
                idx += 1
            tau = tau_min + idx
        else:
            tau = tau_min + int(np.argmin(search))
            if cmndf[tau] > 0.5:
                continue

        # Parabolic refinement around tau
        if 1 <= tau < tau_max:
            y0, y1, y2 = cmndf[tau - 1], cmndf[tau], cmndf[tau + 1]
            denom_p = y0 - 2.0 * y1 + y2
            if abs(denom_p) > 1e-12:
                delta = 0.5 * (y0 - y2) / denom_p
                tau_refined = tau + float(np.clip(delta, -0.5, 0.5))
            else:
                tau_refined = float(tau)
        else:
            tau_refined = float(tau)

        if tau_refined > 0:
            f_est = fs / tau_refined
            if fmin <= f_est <= fmax:
                f0_estimates.append(float(f_est))

    if not f0_estimates:
        return {"measured_hz": 0.0, "measured_std_hz": 0.0, "method": "yin", "voiced_frames": 0}

    return {
        "measured_hz": round(float(np.median(f0_estimates)), 4),
        "measured_std_hz": round(float(np.std(f0_estimates)), 4),
        "method": "yin",
        "voiced_frames": len(f0_estimates),
    }


def measure_formants_lpc(
    x: np.ndarray,
    fs: int,
    *,
    fs_lpc: int = 8000,
    order: int | None = None,
    pre_emphasis: float | None = 0.5,
    frame_ms: float = 50.0,
    hop_ms: float = 10.0,
    f0_hz: float | None = None,
) -> dict:
    """Estimate F1..F3 via frame-based LPC root solving and flag harmonic locking."""
    x = np.asarray(x, dtype=np.float64)
    if x.ndim > 1:
        x = np.mean(x, axis=1)
    if order is None:
        order = 2 + fs_lpc // 1000

    if f0_hz is None:
        f0_hz = measure_f0(x, fs)["measured_hz"]

    if fs != fs_lpc:
        g = gcd(int(fs), int(fs_lpc))
        xr = signal.resample_poly(x, fs_lpc // g, fs // g)
    else:
        xr = x.copy()

    if pre_emphasis is not None and pre_emphasis != 0.0 and len(xr) > 1:
        xr = np.append(xr[0], xr[1:] - float(pre_emphasis) * xr[:-1])

    flen = max(order + 8, int(round(frame_ms * 1e-3 * fs_lpc)))
    hop = max(1, int(round(hop_ms * 1e-3 * fs_lpc)))
    win = np.hamming(flen)
    rms_thresh = max(1e-6, 0.03 * float(np.sqrt(np.mean(xr**2))))

    f_tracks: list[list[float]] = [[], [], []]
    b_tracks: list[list[float]] = [[], [], []]
    lock_counts = [0, 0, 0]
    n_voiced = 0

    for start in range(0, max(1, len(xr) - flen + 1), hop):
        seg = xr[start : start + flen]
        if len(seg) < flen:
            break
        if float(np.sqrt(np.mean(seg**2))) < rms_thresh:
            continue
        frame = seg * win
        r_full = np.correlate(frame, frame, mode="full")
        r = r_full[flen - 1 : flen + order]
        if r[0] <= 1e-18:
            continue
        try:
            a_coef = linalg.solve_toeplitz((r[:-1], r[:-1]), -r[1:])
        except np.linalg.LinAlgError:
            continue

        poly = np.concatenate(([1.0], a_coef))
        roots = np.roots(poly)
        roots = roots[np.imag(roots) > 0]
        if roots.size == 0:
            continue

        freqs = np.angle(roots) * fs_lpc / (2.0 * np.pi)
        bws = -np.log(np.maximum(np.abs(roots), 1e-12)) * fs_lpc / np.pi
        mask = (freqs > 90.0) & (freqs < 3500.0) & (bws > 0.0) & (bws < 600.0)
        cands = sorted(zip(freqs[mask], bws[mask]), key=lambda item: item[0])
        if not cands:
            continue

        n_voiced += 1
        for idx, (cf, cb) in enumerate(cands[:3]):
            f_tracks[idx].append(float(cf))
            b_tracks[idx].append(float(cb))
            if f0_hz and f0_hz > 0.0:
                k = round(cf / f0_hz)
                if k >= 1 and abs(cf - k * f0_hz) <= 5.0:
                    lock_counts[idx] += 1

    f_med = [round(float(np.median(trk)), 2) for trk in f_tracks if trk]
    b_med = [round(float(np.median(trk)), 2) for trk in b_tracks if trk]
    lock_fracs = [
        round(float(cnt / max(1, len(trk))), 4) if trk else 0.0
        for cnt, trk in zip(lock_counts, f_tracks)
    ]

    low_reasons: list[str] = []
    if any(frac > 0.50 for frac in lock_fracs):
        low_reasons.append("harmonic_locking_gt_50pct")
    if f0_hz is not None and 0.0 < f0_hz <= 75.0:
        low_reasons.append("low_f0_dense_harmonics_le_75hz")

    return {
        "f_hz": f_med,
        "bw_hz": b_med,
        "order": int(order),
        "fs_lpc": int(fs_lpc),
        "pre_emphasis": None if pre_emphasis is None else float(pre_emphasis),
        "reliability": "low" if low_reasons else "ok",
        "reliability_reasons": low_reasons,
        "harmonic_lock_fractions": lock_fracs,
        "voiced_frames": n_voiced,
    }


def _parabolic_peak(freqs: np.ndarray, y_db: np.ndarray, idx: int) -> tuple[float, float]:
    if 1 <= idx < len(freqs) - 1:
        y0, y1, y2 = y_db[idx - 1], y_db[idx], y_db[idx + 1]
        denom = y0 - 2.0 * y1 + y2
        if abs(denom) > 1e-12:
            delta = float(np.clip(0.5 * (y0 - y2) / denom, -0.5, 0.5))
            df = float(freqs[idx + 1] - freqs[idx])
            f_peak = float(freqs[idx] + delta * df)
            y_peak = float(y1 - 0.25 * (y0 - y2) * delta)
            return f_peak, y_peak
    return float(freqs[idx]), float(y_db[idx])


def _bandwidth_3db(freqs: np.ndarray, y_db: np.ndarray, idx: int, peak_db: float) -> float:
    target = peak_db - 3.0
    df = float(freqs[1] - freqs[0]) if len(freqs) > 1 else 1.0

    left = idx
    while left > 0 and y_db[left] > target:
        left -= 1
    if left < idx and y_db[left] <= target < y_db[left + 1]:
        t = (target - y_db[left]) / max(y_db[left + 1] - y_db[left], 1e-12)
        f_left = float(freqs[left] + t * (freqs[left + 1] - freqs[left]))
    else:
        f_left = float(freqs[left])

    right = idx
    n = len(freqs)
    while right < n - 1 and y_db[right] > target:
        right += 1
    if right > idx and y_db[right] <= target < y_db[right - 1]:
        t = (target - y_db[right]) / max(y_db[right - 1] - y_db[right], 1e-12)
        f_right = float(freqs[right] - t * (freqs[right] - freqs[right - 1]))
    else:
        f_right = float(freqs[right])

    return float(max(f_right - f_left, df))


def measure_transfer_function(
    source: np.ndarray,
    filtered: np.ndarray,
    fs: int,
    *,
    nperseg: int = 16384,
    min_coherence: float = 0.8,
    fmin_hz: float = 90.0,
    fmax_hz: float = 3800.0,
    savgol_window: int = 9,
    prominence_db: float = 1.5,
) -> dict:
    """Measure H(f) = P_sf / P_ss and coherence gamma^2 via Welch cross-spectra."""
    s = np.asarray(source, dtype=np.float64)
    y = np.asarray(filtered, dtype=np.float64)
    n = min(len(s), len(y))
    s = s[:n]
    y = y[:n]

    nseg = min(int(nperseg), n)
    if nseg < 64:
        nseg = n
    noverlap = nseg // 2

    freqs, p_ss = signal.welch(s, fs=fs, window="hann", nperseg=nseg, noverlap=noverlap)
    _, p_ff = signal.welch(y, fs=fs, window="hann", nperseg=nseg, noverlap=noverlap)
    _, p_sf = signal.csd(s, y, fs=fs, window="hann", nperseg=nseg, noverlap=noverlap)

    h_complex = p_sf / (p_ss + 1e-30)
    coherence = np.real(np.abs(p_sf) ** 2 / (p_ss * p_ff + 1e-30))
    coherence = np.clip(coherence, 0.0, 1.0)
    h_db = 20.0 * np.log10(np.abs(h_complex) + 1e-12)

    df = float(freqs[1] - freqs[0]) if len(freqs) > 1 else 1.0
    win_len = savgol_window if savgol_window % 2 == 1 else savgol_window + 1
    win_len = min(win_len, len(h_db) - (1 - len(h_db) % 2))
    if win_len >= 5:
        h_db_smooth = signal.savgol_filter(h_db, window_length=win_len, polyorder=2)
    else:
        h_db_smooth = h_db.copy()

    min_dist_bins = max(2, int(round(100.0 / max(df, 1e-6))))
    band_mask = (freqs >= fmin_hz) & (freqs <= fmax_hz)
    h_band = h_db_smooth[band_mask]
    f_band = freqs[band_mask]
    c_band = coherence[band_mask]

    peak_indices, _ = signal.find_peaks(h_band, prominence=prominence_db, distance=min_dist_bins)
    # Also map back to full-grid indices for clean bandwidth calculation
    band_offset = int(np.argmax(band_mask))

    peaks_out: list[dict] = []
    for p_local in peak_indices:
        p_global = band_offset + int(p_local)
        lo_c = max(0, p_global - 1)
        hi_c = min(len(coherence), p_global + 2)
        # Use 3-bin local median coherence to avoid single-bin finite-window ringing bias at high-Q peaks
        coh_val = float(np.median(coherence[lo_c:hi_c]))
        if coh_val < min_coherence:
            continue
        f_peak, mag_peak = _parabolic_peak(freqs, h_db_smooth, p_global)
        bw_3db = _bandwidth_3db(freqs, h_db_smooth, p_global, mag_peak)
        q_val = float(f_peak / bw_3db) if bw_3db > 0 else 0.0
        peaks_out.append(
            {
                "f_hz": round(f_peak, 3),
                "bw_3db_hz": round(bw_3db, 3),
                "q": round(q_val, 3),
                "coherence": round(coh_val, 4),
                "gain_db": round(mag_peak, 3),
            }
        )

    return {
        "peaks": peaks_out,
        "freqs_hz": freqs,
        "h_db": h_db,
        "h_db_smooth": h_db_smooth,
        "coherence": coherence,
        "nperseg": int(nseg),
        "df_hz": df,
        "min_coherence": float(min_coherence),
    }


def default_frequency_bands(fs: int) -> list[tuple[float, float]]:
    nyq = float(fs) / 2.0
    return [
        (0.0, 100.0),
        (100.0, 300.0),
        (300.0, 420.0),
        (420.0, 800.0),
        (800.0, 1500.0),
        (1500.0, 3000.0),
        (3000.0, nyq),
    ]


def _band_label(lo: float, hi: float) -> str:
    lo_s = str(int(round(lo))) if abs(lo - round(lo)) < 1e-6 else f"{lo:.1f}"
    hi_s = str(int(round(hi))) if abs(hi - round(hi)) < 1e-6 else f"{hi:.1f}"
    return f"{lo_s}-{hi_s}"


def band_energy_fractions(
    x: np.ndarray,
    fs: int,
    bands: Iterable[tuple[float, float]] | None = None,
    *,
    nperseg: int = 16384,
) -> dict[str, float]:
    """Compute energy fraction in each frequency band using Welch PSD.

    When default bands are used, they partition [0, fs/2] so fractions sum to 1.0 (+/-1e-6).
    """
    x = np.asarray(x, dtype=np.float64)
    if x.ndim > 1:
        x = np.mean(x, axis=1)
    band_list = list(bands) if bands is not None else default_frequency_bands(fs)

    nseg = min(int(nperseg), len(x))
    freqs, psd = signal.welch(x, fs=fs, window="hann", nperseg=nseg, noverlap=nseg // 2)
    total = float(np.sum(psd))
    if total <= 0.0:
        return {_band_label(lo, hi): 0.0 for lo, hi in band_list}

    out: dict[str, float] = {}
    last_idx = len(band_list) - 1
    for idx, (lo, hi) in enumerate(band_list):
        if idx == last_idx:
            mask = (freqs >= lo) & (freqs <= hi)
        else:
            mask = (freqs >= lo) & (freqs < hi)
        out[_band_label(lo, hi)] = float(np.sum(psd[mask]) / total)
    return out


def harmonic_band_fractions(
    x: np.ndarray,
    fs: int,
    f0_hz: float,
    *,
    half_width_hz: float = 8.0,
    nperseg: int = 16384,
) -> dict[str, float]:
    """Measure energy fraction within f0 +/- 8 Hz and 2*f0 +/- 8 Hz."""
    x = np.asarray(x, dtype=np.float64)
    if x.ndim > 1:
        x = np.mean(x, axis=1)
    nseg = min(int(nperseg), len(x))
    freqs, psd = signal.welch(x, fs=fs, window="hann", nperseg=nseg, noverlap=nseg // 2)
    total = float(np.sum(psd)) or 1.0
    m1 = (freqs >= max(0.0, f0_hz - half_width_hz)) & (freqs <= f0_hz + half_width_hz)
    m2 = (freqs >= max(0.0, 2.0 * f0_hz - half_width_hz)) & (freqs <= 2.0 * f0_hz + half_width_hz)
    return {
        "f0_pm_8": float(np.sum(psd[m1]) / total),
        "2f0_pm_8": float(np.sum(psd[m2]) / total),
    }


def band_fraction_for_window(
    x: np.ndarray,
    fs: int,
    f_low: float,
    f_high: float,
    *,
    nperseg: int = 16384,
) -> float:
    """Compute energy fraction in [f_low, f_high]."""
    x = np.asarray(x, dtype=np.float64)
    if x.ndim > 1:
        x = np.mean(x, axis=1)
    nseg = min(int(nperseg), len(x))
    freqs, psd = signal.welch(x, fs=fs, window="hann", nperseg=nseg, noverlap=nseg // 2)
    total = float(np.sum(psd))
    if total <= 0.0:
        return 0.0
    mask = (freqs >= f_low) & (freqs <= f_high)
    return float(np.sum(psd[mask]) / total)


def compute_shaping_gain_db(
    source: np.ndarray,
    filtered: np.ndarray,
    fs: int,
    f_center_hz: float,
    bandwidth_hz: float,
    *,
    nperseg: int = 16384,
    min_half_width_hz: float | None = None,
) -> tuple[float, float, float]:
    """Return (shaping_gain_db, filtered_band_frac, source_band_frac) for [f_center - half, f_center + half]."""
    n_avail = max(1, min(len(source), len(filtered)))
    nseg = min(int(nperseg), n_avail)
    min_hw = (2.0 * float(fs) / max(nseg, 1)) if min_half_width_hz is None else float(min_half_width_hz)
    half = max(min_hw, 0.5 * float(bandwidth_hz))
    f_lo = max(0.0, float(f_center_hz) - half)
    f_hi = float(f_center_hz) + half
    e_flt = band_fraction_for_window(filtered, fs, f_lo, f_hi, nperseg=nperseg)
    e_src = band_fraction_for_window(source, fs, f_lo, f_hi, nperseg=nperseg)
    gain_db = 10.0 * np.log10(max(e_flt, 1e-12) / max(e_src, 1e-12))
    return float(gain_db), float(e_flt), float(e_src)


def compute_shaping_metrics(
    source: np.ndarray,
    filtered: np.ndarray,
    fs: int,
    f1_hz: float,
    b1_hz: float,
    *,
    nperseg: int = 16384,
) -> dict:
    """Compute both F1 +/- B1/2 (min 2 Welch bins) and legacy F1 +/- 60 Hz shaping metrics."""
    n_avail = max(1, min(len(source), len(filtered)))
    nseg = min(int(nperseg), n_avail)
    min_hw = 2.0 * float(fs) / max(nseg, 1)

    half_b3db = max(min_hw, 0.5 * float(b1_hz))
    lo_b3db = max(0.0, float(f1_hz) - half_b3db)
    hi_b3db = float(f1_hz) + half_b3db
    g_b3db, ef_b3db, es_b3db = compute_shaping_gain_db(
        source, filtered, fs, f1_hz, b1_hz, nperseg=nperseg, min_half_width_hz=min_hw
    )

    lo_pm60 = max(0.0, float(f1_hz) - 60.0)
    hi_pm60 = float(f1_hz) + 60.0
    g_pm60, ef_pm60, es_pm60 = compute_shaping_gain_db(
        source, filtered, fs, f1_hz, 120.0, nperseg=nperseg, min_half_width_hz=60.0
    )

    return {
        "b3db": {
            "shaping_gain_db_F1_b3db": round(g_b3db, 3),
            "band_hz": [round(lo_b3db, 3), round(hi_b3db, 3)],
            "f1_band_energy_fraction_measured_b3db": round(ef_b3db, 6),
            "source_band_energy_fraction_measured_b3db": round(es_b3db, 6),
            "min_half_width_hz": round(min_hw, 4),
        },
        "pm60": {
            "shaping_gain_db_F1_pm60": round(g_pm60, 3),
            "band_hz": [round(lo_pm60, 3), round(hi_pm60, 3)],
            "f1_band_energy_fraction_measured_pm60": round(ef_pm60, 6),
            "source_band_energy_fraction_measured_pm60": round(es_pm60, 6),
        },
        "shaping_gain_db_F1_b3db": round(g_b3db, 3),
        "shaping_gain_db_F1_pm60": round(g_pm60, 3),
        "f1_band_energy_fraction_measured_pm60": round(ef_pm60, 6),
    }


def measure_source_spectral_tilt(
    source: np.ndarray,
    fs: int,
    f0_hz: float,
    *,
    fmin_hz: float = 100.0,
    fmax_hz: float = 2000.0,
) -> float:
    """Measure harmonic spectral tilt (dB/octave) via linear regression on log2(f) in [fmin, fmax]."""
    s = np.asarray(source, dtype=np.float64)
    if s.ndim > 1:
        s = np.mean(s, axis=1)
    if len(s) < 64 or f0_hz <= 0.0:
        return 0.0

    win = np.hanning(len(s))
    spec = np.abs(np.fft.rfft(s * win))
    df = float(fs) / len(s)

    harms: list[float] = []
    amps_db: list[float] = []
    k_start = max(2, int(np.ceil(fmin_hz / f0_hz)))
    k_end = int(np.floor(fmax_hz / f0_hz))
    half_bins = max(1, int(round(4.0 / df)))

    for k in range(k_start, k_end + 1):
        hf = k * f0_hz
        idx = int(round(hf / df))
        lo = max(0, idx - half_bins)
        hi = min(len(spec), idx + half_bins + 1)
        if lo >= hi:
            continue
        pk = float(np.max(spec[lo:hi]))
        harms.append(hf)
        amps_db.append(float(20.0 * np.log10(pk + 1e-12)))

    if len(harms) < 2:
        return 0.0
    slope, _ = np.polyfit(np.log2(np.asarray(harms)), np.asarray(amps_db), 1)
    return round(float(slope), 3)


def plot_transfer_function(tf_data: dict, out_png: Path, *, fmax_plot: float = 3500.0) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    freqs = tf_data["freqs_hz"]
    h_db = tf_data["h_db"]
    h_smooth = tf_data["h_db_smooth"]
    coh = tf_data["coherence"]
    peaks = tf_data["peaks"]
    mask = freqs <= fmax_plot

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6), sharex=True, gridspec_kw={"height_ratios": [2.2, 1.0]})
    ax1.plot(freqs[mask], h_db[mask], color="0.75", lw=0.8, label="|H(f)| raw")
    ax1.plot(freqs[mask], h_smooth[mask], color="navy", lw=1.6, label="|H(f)| smoothed")
    for idx, pk in enumerate(peaks):
        if pk["f_hz"] <= fmax_plot:
            ax1.axvline(pk["f_hz"], color="crimson", ls="--", lw=1.0, alpha=0.7)
            ax1.plot(pk["f_hz"], pk["gain_db"], "ro", ms=5)
            ax1.annotate(
                f"F{idx+1}: {pk['f_hz']:.1f} Hz\n(Q={pk['q']:.1f})",
                xy=(pk["f_hz"], pk["gain_db"]),
                xytext=(6, 6),
                textcoords="offset points",
                fontsize=8,
                color="crimson",
            )
    ax1.set_ylabel("|H(f)| = |P_sf / P_ss| (dB)")
    ax1.set_title("Measured Acoustic Transfer Function H(f) & Coherence")
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc="upper right", fontsize=8)

    ax2.plot(freqs[mask], coh[mask], color="teal", lw=1.2, label="Coherence γ²(f)")
    ax2.axhline(tf_data["min_coherence"], color="darkorange", ls=":", lw=1.2, label=f"Threshold γ²={tf_data['min_coherence']}")
    ax2.set_ylim(0.0, 1.05)
    ax2.set_xlabel("Frequency (Hz)")
    ax2.set_ylabel("γ²(f)")
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc="lower right", fontsize=8)

    fig.tight_layout()
    fig.savefig(out_png, dpi=96)
    plt.close(fig)


def measure_run_dir(
    run_dir: Path,
    *,
    lpc_order: int | None = None,
    pre_emphasis: float | None = 0.5,
    write_outputs: bool = True,
) -> dict:
    """Measure source.wav and filtered.wav in a run directory without reading target priors."""
    src_path = run_dir / "source.wav"
    flt_path = run_dir / "filtered.wav"
    sr_s, source, clip_s = read_wav_float(src_path)
    sr_f, filtered, clip_f = read_wav_float(flt_path)
    if sr_s != sr_f:
        raise ValueError(f"Sample rate mismatch: source={sr_s} Hz vs filtered={sr_f} Hz")
    sr = sr_s

    f0_src = measure_f0(source, sr)
    f0_flt = measure_f0(filtered, sr)
    f0_val = f0_src["measured_hz"] or f0_flt["measured_hz"]

    tf_data = measure_transfer_function(source, filtered, sr)
    lpc_res = measure_formants_lpc(filtered, sr, order=lpc_order, pre_emphasis=pre_emphasis, f0_hz=f0_val)

    bands_flt = band_energy_fractions(filtered, sr)
    bands_src = band_energy_fractions(source, sr)
    harm_flt = harmonic_band_fractions(filtered, sr, f0_val) if f0_val > 0 else {"f0_pm_8": 0.0, "2f0_pm_8": 0.0}
    harm_src = harmonic_band_fractions(source, sr, f0_val) if f0_val > 0 else {"f0_pm_8": 0.0, "2f0_pm_8": 0.0}
    tilt = measure_source_spectral_tilt(source, sr, f0_val)

    peak_src_dbfs = round(float(20.0 * np.log10(max(float(np.max(np.abs(source))), 1e-12))), 3)
    peak_flt_dbfs = round(float(20.0 * np.log10(max(float(np.max(np.abs(filtered))), 1e-12))), 3)

    shaping_measured = None
    if tf_data["peaks"]:
        p1 = tf_data["peaks"][0]
        sh_all = compute_shaping_metrics(source, filtered, sr, p1["f_hz"], p1["bw_3db_hz"])
        shaping_measured = {
            "f1_hz": p1["f_hz"],
            "shaping_gain_db": sh_all["shaping_gain_db_F1_pm60"],
            "filtered_band_fraction": sh_all["f1_band_energy_fraction_measured_pm60"],
            "source_band_fraction": sh_all["pm60"]["source_band_energy_fraction_measured_pm60"],
            "shaping_gain_db_F1_b3db": sh_all["shaping_gain_db_F1_b3db"],
            "b3db": sh_all["b3db"],
            "pm60": sh_all["pm60"],
        }

    measured = {
        "run_id": run_dir.name,
        "sample_rate": sr,
        "f0": f0_src if f0_src["measured_hz"] > 0 else f0_flt,
        "f0_filtered": f0_flt,
        "formants_measured": {
            "transfer_function": tf_data["peaks"],
            "lpc": lpc_res,
        },
        "band_energy_fractions": {
            "filtered": bands_flt,
            "source": bands_src,
        },
        "harmonic_band_fractions": {
            "filtered": harm_flt,
            "source": harm_src,
        },
        "shaping_at_measured_f1": shaping_measured,
        "source_spectral_tilt_db_per_oct": tilt,
        "audio": {
            "source_peak_dbfs": peak_src_dbfs,
            "filtered_peak_dbfs": peak_flt_dbfs,
            "clipping_detected": bool(clip_s or clip_f),
        },
    }

    if write_outputs:
        (run_dir / "measured.json").write_text(json.dumps(measured, indent=2) + "\n", encoding="utf-8")
        plot_transfer_function(tf_data, run_dir / "transfer_function.png")

    return measured


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="M1: Audio-only acoustic measurement")
    parser.add_argument("--run-dir", required=True, help="Run output directory containing source.wav and filtered.wav")
    parser.add_argument("--lpc-order", type=int, default=None, help="LPC order (default: 2 + fs_lpc/1000 = 10)")
    parser.add_argument("--no-pre-emphasis", action="store_true", help="Disable LPC pre-emphasis")
    args = parser.parse_args(argv)

    run_dir = Path(args.run_dir)
    pre = None if args.no_pre_emphasis else 0.5
    result = measure_run_dir(run_dir, lpc_order=args.lpc_order, pre_emphasis=pre, write_outputs=True)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

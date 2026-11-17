"""M4 — Interactive two-mass glottal source + Kelly–Lochbaum waveguide synthesis.

Implements:
  - Interactive Steinecke–Herzel (1995) two-mass glottal model in SI units with
    vocal-fold collision and supraglottal pressure feedback P1(t).
  - Kelly–Lochbaum (KL) digital waveguide built from the same A(s) profile as M2,
    with Levine–Schwinger / flanged end correction and radiation low-pass damping.
  - Mandatory KL impulse-response vs M2 TMM verification (< 2% error below 2 kHz).
  - Legacy IIR synthesis path (--synth legacy_iir) for honest baseline comparison.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
from scipy import signal
from scipy.optimize import nnls

from pipeline.tmm_lossy import Air, WallModel

logger = logging.getLogger(__name__)


def normalize_audio(x: np.ndarray, target_dbfs: float = -1.0) -> tuple[np.ndarray, float, bool]:
    """Peak-normalize signal to target_dbfs (-1.0 dBFS => 0.89125) and measure clipping."""
    arr = np.asarray(x, dtype=np.float64)
    if arr.size == 0:
        return arr, -np.inf, False
    arr = arr - float(np.mean(arr))
    peak = float(np.max(np.abs(arr)))
    target_amp = 10.0 ** (target_dbfs / 20.0)
    if peak > 0.0:
        out = arr * (target_amp / peak)
    else:
        out = arr.copy()
    max_abs = float(np.max(np.abs(out)))
    clipped = bool(max_abs >= 0.999)
    assert max_abs < 0.999, f"Normalized audio exceeds safety ceiling 0.999 (got {max_abs:.6f})"
    peak_dbfs = float(20.0 * np.log10(max(max_abs, 1e-12)))
    return out, round(peak_dbfs, 3), clipped


def check_clipping_raw(x: np.ndarray, threshold: float = 0.999) -> bool:
    """Measure whether a raw or PCM-scaled waveform exhibits clipping."""
    arr = np.asarray(x, dtype=np.float64)
    return bool(np.any(np.abs(arr) >= threshold))


def _radiation_cutoff_hz(a_out_m: float, c: float, radiation: str) -> tuple[float, float, float]:
    """Return (end_correction_m, f_rad_cutoff_hz, R_L0) from outlet radius and radiation model."""
    if radiation == "ideal_open":
        return 0.0, 20000.0, -0.998
    if radiation == "flanged":
        return 0.8216 * a_out_m, (c / (2.0 * np.pi * max(a_out_m, 1e-4))) * 1.2, -0.994
    # Unflanged Levine–Schwinger default
    return 0.6133 * a_out_m, (c / (2.0 * np.pi * max(a_out_m, 1e-4))) * 1.4, -0.995


def _fit_kirchhoff_section_biquads(
    r_sec_m: np.ndarray,
    dx_m: float,
    fs_wg: float,
    air: Air,
    *,
    extra_damping_factor: float = 1.0,
) -> dict:
    """Fit stable minimum-phase 2-pole loss biquads H_j(z) to exp(-alpha_j(f)*dx) over 50 Hz-4 kHz."""
    f_grid = np.geomspace(50.0, 4000.0, 64)
    omega = 2.0 * np.pi * f_grid
    wall_factor = (
        float(extra_damping_factor)
        * np.sqrt(omega * air.mu / (2.0 * air.rho))
        * (1.0 + (air.gamma - 1.0) / np.sqrt(air.Pr))
    )
    p1 = float(np.exp(-2.0 * np.pi * 250.0 / fs_wg))
    p2 = float(np.exp(-2.0 * np.pi * 1800.0 / fs_wg))
    z_inv = np.exp(-1j * omega / fs_wg)
    h_hp1 = (1.0 - z_inv) / (1.0 - p1 * z_inv)
    h_hp2 = (1.0 - z_inv) / (1.0 - p2 * z_inv)

    # Fit attenuation per unit (1/r) on the log-frequency grid using non-negative basis weights
    target_unit = (wall_factor / air.c) * dx_m
    design = np.column_stack([np.ones_like(f_grid), np.real(h_hp1), np.real(h_hp2)])
    c_unit, _ = nnls(design, target_unit)

    c_sec = c_unit[None, :] / r_sec_m[:, None]
    # Keep filter strictly causal, minimum-phase, and passive if extreme narrow section occurs
    row_sum = np.sum(c_sec, axis=1)
    scale = np.where(row_sum > 0.90, 0.90 / np.maximum(row_sum, 1e-12), 1.0)
    c_sec = c_sec * scale[:, None]

    c0, c1, c2 = c_sec[:, 0], c_sec[:, 1], c_sec[:, 2]
    b0 = 1.0 - c0 - c1 - c2
    b1 = -(p1 + p2) * (1.0 - c0) + c1 * (1.0 + p2) + c2 * (1.0 + p1)
    b2 = (p1 * p2) * (1.0 - c0) - c1 * p2 - c2 * p1
    a1 = -(p1 + p2)
    a2 = p1 * p2

    num = b0[:, None] + b1[:, None] * z_inv[None, :] + b2[:, None] * (z_inv[None, :] ** 2)
    den = 1.0 + a1 * z_inv + a2 * (z_inv ** 2)
    mag_db = 20.0 * np.log10(np.maximum(np.abs(num / den[None, :]), 1e-12))
    targ_db = -20.0 / np.log(10.0) * (wall_factor[None, :] / (r_sec_m[:, None] * air.c)) * dx_m
    rms_db = np.sqrt(np.mean((mag_db - targ_db) ** 2, axis=1))

    return {
        "b0": b0,
        "b1": b1,
        "b2": b2,
        "a1": a1,
        "a2": a2,
        "fit_rms_db": {
            "median": round(float(np.median(rms_db)), 6),
            "max": round(float(np.max(rms_db)), 6),
            "f_grid_hz": [50.0, 4000.0],
        },
    }


def _radiation_reflection_biquad(
    a_out_m: float,
    c: float,
    fs_wg: float,
    radiation: str,
    area_ratio_beta: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Passive 2-pole bilinear transform of radiation reflection R_L(s) = -(D - beta*N)/(D + beta*N)."""
    if radiation == "ideal_open":
        return np.array([-1.0, 0.0, 0.0], dtype=np.float64), np.array([1.0, 0.0, 0.0], dtype=np.float64)

    tau = max(a_out_m, 1e-4) / c
    beta = max(float(area_ratio_beta), 1e-6)
    if radiation == "flanged":
        c1_n, c2_n = 0.8216 * beta, 0.50 * beta
        c1_d, c2_d = 1.217, 0.50
    elif radiation == "unflanged":
        # Silva et al. (2009) passive 2nd-order rational fit matching Levine-Schwinger at low ka
        c1_n, c2_n = 0.6133 * beta, 0.2406 * beta
        c1_d, c2_d = 0.80, 0.2406
    else:
        raise ValueError(f"Unsupported radiation model '{radiation}'.")

    b_s = np.array([-(c2_d - c2_n) * (tau ** 2), -(c1_d - c1_n) * tau, -1.0], dtype=np.float64)
    a_s = np.array([(c2_d + c2_n) * (tau ** 2), (c1_d + c1_n) * tau, 1.0], dtype=np.float64)
    b_z, a_z = signal.bilinear(b_s, a_s, fs=fs_wg)
    return np.asarray(b_z, dtype=np.float64), np.asarray(a_z, dtype=np.float64)


def build_kl_waveguide(
    s_mm: np.ndarray,
    area_mm2: np.ndarray,
    air: Air,
    *,
    dh_mm: np.ndarray | None = None,
    perimeter_mm: np.ndarray | None = None,
    radiation: str = "unflanged",
    sr: int = 44100,
    oversample: int = 2,
    kl_loss_model: str = "kirchhoff_fit",
    extra_damping_factor: float = 1.0,
    wall: WallModel | None = None,
    f1_hint_hz: float = 400.0,
) -> dict:
    """Discretize A(s) into a passive Kelly–Lochbaum waveguide at fs_wg >= oversample*sr."""
    if kl_loss_model not in ("legacy", "kirchhoff_fit"):
        raise ValueError(f"Unsupported kl_loss_model '{kl_loss_model}'. Expected 'legacy' or 'kirchhoff_fit'.")

    s_m = np.asarray(s_mm, dtype=np.float64) * 1e-3
    area_m2 = np.asarray(area_mm2, dtype=np.float64) * 1e-6
    l_phys = float(s_m[-1] - s_m[0])
    a_out_m = float(np.sqrt(area_m2[-1] / np.pi))
    fs_target = sr * max(2, int(oversample))

    if kl_loss_model == "legacy":
        end_corr_m, f_rad_hz, r_l0 = _radiation_cutoff_hz(a_out_m, air.c, radiation)
        n0 = max(16, int(round(2.0 * l_phys * fs_target / air.c)))
        fs0 = n0 * air.c / max(l_phys, 1e-6)
        lam0 = float(np.exp(-2.0 * np.pi * f_rad_hz / fs0))
        tau_lp = lam0 / max((1.0 - lam0) * fs0, 1e-9)

        l_wg = l_phys + max(0.0, end_corr_m - 0.5 * air.c * tau_lp)
        n_sec = max(16, int(round(2.0 * l_wg * fs_target / air.c)))
        fs_wg = n_sec * air.c / l_wg
        dx = l_wg / n_sec

        s_sec = s_m[0] + (np.arange(n_sec) + 0.5) * dx
        area_sec = np.interp(np.minimum(s_sec, s_m[-1]), s_m, area_m2)
        r_junc = (area_sec[:-1] - area_sec[1:]) / (area_sec[:-1] + area_sec[1:])
        z0_glottis = float(air.rho * air.c / area_sec[0])
        lam = float(np.exp(-2.0 * np.pi * f_rad_hz / fs_wg))

        alpha_roundtrip = 0.024 * float(extra_damping_factor)
        loss_sec = float(np.exp(-alpha_roundtrip / (2.0 * n_sec)))
        r_g_closed = 0.999

        return {
            "N": int(n_sec),
            "fs_wg": float(fs_wg),
            "dx_m": float(dx),
            "L_wg_m": float(l_wg),
            "r_junc": r_junc,
            "area_sec": area_sec,
            "Z0_glottis": z0_glottis,
            "loss_sec": loss_sec,
            "R_g": r_g_closed,
            "R_L0": float(r_l0),
            "lam": lam,
            "f_rad_hz": float(f_rad_hz),
            "kl_loss_model": "legacy",
            "extra_damping_factor": float(extra_damping_factor),
            "fit_rms_db": None,
            "radiation_notes": [],
        }

    # kirchhoff_fit: exact physical length (radiation biquad models Levine-Schwinger/Silva end correction directly)
    n_sec = max(16, int(np.ceil(l_phys * fs_target / air.c)))
    fs_wg = n_sec * air.c / max(l_phys, 1e-6)
    dx = l_phys / n_sec
    s_sec = s_m[0] + (np.arange(n_sec) + 0.5) * dx
    area_sec = np.interp(np.minimum(s_sec, s_m[-1]), s_m, area_m2)

    if dh_mm is not None:
        dh_m = np.asarray(dh_mm, dtype=np.float64) * 1e-3
        dh_sec = np.interp(np.minimum(s_sec, s_m[-1]), s_m, dh_m)
    elif perimeter_mm is not None:
        p_m = np.asarray(perimeter_mm, dtype=np.float64) * 1e-3
        p_sec = np.interp(np.minimum(s_sec, s_m[-1]), s_m, p_m)
        dh_sec = 4.0 * area_sec / np.maximum(p_sec, 1e-6)
    else:
        dh_sec = 2.0 * np.sqrt(area_sec / np.pi)

    r_sec = np.maximum(0.5 * dh_sec, 1e-4)
    r_junc = (area_sec[:-1] - area_sec[1:]) / (area_sec[:-1] + area_sec[1:])
    z0_glottis = float(air.rho * air.c / area_sec[0])

    loss_biquad = _fit_kirchhoff_section_biquads(
        r_sec, dx, fs_wg, air, extra_damping_factor=extra_damping_factor
    )
    beta_rad = float(area_sec[-1] / area_m2[-1])
    b_rad, a_rad = _radiation_reflection_biquad(a_out_m, air.c, fs_wg, radiation, beta_rad)
    rad_note = (
        f"Refined KL outlet radiation reflection filter (2-pole bilinear Silva/flanged fit, "
        f"beta_out={beta_rad:.4f}) to match TMM F1-F3 radiation damping."
    )
    logger.info(rad_note)

    _, f_rad_hz, r_l0 = _radiation_cutoff_hz(a_out_m, air.c, radiation)
    lam_rad = float(np.exp(-2.0 * np.pi * f_rad_hz / fs_wg))
    loss_sec_f1 = float((1.0 - 0.024) ** (1.0 / (2.0 * n_sec)))

    wall_branch: dict | None = None
    if wall is not None and not wall.is_rigid():
        # Prewarped bilinear 3-port parallel wall shunt (zero delay-free loop, strictly passive)
        f_warp = float(np.clip(f1_hint_hz, 80.0, 0.25 * fs_wg))
        kappa = (2.0 * np.pi * f_warp) / float(np.tan(np.pi * f_warp / fs_wg))
        r0 = float(wall.resistance_per_area)
        m0 = kappa * float(wall.mass_per_area)
        k0 = float(wall.stiffness_per_area) / kappa
        d0 = max(r0 + m0 + k0, 1e-12)

        perim_sec = 2.0 * area_sec / r_sec
        perim_junc = 0.5 * (perim_sec[:-1] + perim_sec[1:]) * (n_sec / max(n_sec - 1, 1))
        y1 = area_sec[:-1] / (air.rho * air.c)
        y2 = area_sec[1:] / (air.rho * air.c)
        gw = (perim_junc * dx) / d0
        ytot = y1 + y2 + gw
        wall_branch = {
            "alpha1": 2.0 * y1 / ytot,
            "alpha2": 2.0 * y2 / ytot,
            "alphaw": 2.0 * gw / ytot,
            "g1": float((k0 - m0) / d0),
            "g2": float((m0 + k0) / d0),
            "g3": float(r0 / d0),
            "f_warp_hz": f_warp,
            "has_stiffness": bool(wall.stiffness_per_area > 0.0),
        }

    return {
        "N": int(n_sec),
        "fs_wg": float(fs_wg),
        "dx_m": float(dx),
        "L_wg_m": float(l_phys),
        "r_junc": r_junc,
        "area_sec": area_sec,
        "r_sec": r_sec,
        "Z0_glottis": z0_glottis,
        "loss_sec": loss_sec_f1,
        "loss_biquad": loss_biquad,
        "b_rad": b_rad,
        "a_rad": a_rad,
        "R_g": 1.0,
        "R_L0": float(r_l0),
        "lam": lam_rad,
        "f_rad_hz": float(f_rad_hz),
        "kl_loss_model": "kirchhoff_fit",
        "extra_damping_factor": float(extra_damping_factor),
        "wall_branch": wall_branch,
        "fit_rms_db": loss_biquad["fit_rms_db"],
        "radiation_notes": [rad_note],
    }


def _run_kl_loop(
    wg: dict,
    n_steps: int,
    *,
    impulse: bool = False,
    glottal: dict | None = None,
) -> tuple[np.ndarray, np.ndarray | None]:
    """Execute passive Kelly–Lochbaum time-stepping loop for IR or coupled glottal synthesis."""
    n_sec = wg["N"]
    r_junc = wg["r_junc"]
    one_m_rj = 1.0 - r_junc
    one_p_rj = 1.0 + r_junc
    r_g = wg["R_g"]
    is_legacy = wg.get("kl_loss_model", "legacy") == "legacy"

    fwd = np.zeros(n_sec, dtype=np.float64)
    bwd = np.zeros(n_sec, dtype=np.float64)
    out = np.zeros(n_steps, dtype=np.float64)

    if glottal is not None:
        z0 = wg["Z0_glottis"]
        p_sub_env = glottal["p_sub_env"]
        p_sub_max = glottal["p_sub_max"]
        l_g = glottal["l_g"]
        d1 = glottal["d1"]
        x0 = glottal["x0"]
        m1, m2 = glottal["m1"], glottal["m2"]
        r1, r2 = glottal["r1"], glottal["r2"]
        k1, k2, kc = glottal["k1"], glottal["k2"], glottal["kc"]
        c1, c2 = glottal["c1"], glottal["c2"]
        dt = glottal["dt"]
        rho = glottal["rho"]
        x1, v1, x2, v2 = 0.0001, 0.0, 0.0, 0.0
        ug_out = np.zeros(n_steps, dtype=np.float64)
    else:
        ug_out = None

    if is_legacy or glottal is not None:
        loss_sec = wg["loss_sec"]
        r_l0 = wg["R_L0"]
        lam = wg["lam"]
        one_m_lam = 1.0 - lam
        lp = 0.0
        for k in range(n_steps):
            if impulse:
                inp = 1.0 if k == 0 else 0.0
            elif glottal is not None:
                p_supra = z0 * (fwd[0] + bwd[0])
                lim = 0.35 * p_sub_max
                if p_supra > lim:
                    p_supra = lim
                elif p_supra < -lim:
                    p_supra = -lim
                p_sub = float(p_sub_env[k])
                dp = p_sub - p_supra if p_sub > p_supra else 0.0
                a1_g = 2.0 * l_g * (x0 + x1)
                a2_g = 2.0 * l_g * (x0 + x2)
                a_min = a1_g if a1_g < a2_g else a2_g
                if a_min > 0.0 and a1_g > 0.0:
                    ug = a_min * ((2.0 * dp / rho) ** 0.5)
                    p1_g = p_supra + dp * (1.0 - (a_min / a1_g) ** 2) if a1_g > a_min else p_supra
                else:
                    ug = 0.0
                    p1_g = p_sub if a1_g > 0.0 else 0.0
                f1_force = p1_g * l_g * d1
                fc1 = -c1 * (x0 + x1) if (x0 + x1) < 0.0 else 0.0
                fc2 = -c2 * (x0 + x2) if (x0 + x2) < 0.0 else 0.0
                acc1 = (-r1 * v1 - k1 * x1 - kc * (x1 - x2) + fc1 + f1_force) / m1
                acc2 = (-r2 * v2 - k2 * x2 - kc * (x2 - x1) + fc2) / m2
                v1 += acc1 * dt
                x1 += v1 * dt
                v2 += acc2 * dt
                x2 += v2 * dt
                ug_out[k] = ug
                inp = ug
            else:
                inp = 0.0

            f_in = fwd[:-1] * loss_sec
            b_in = bwd[1:] * loss_sec
            nf_inner = one_m_rj * f_in - r_junc * b_in
            nb_inner = r_junc * f_in + one_p_rj * b_in
            f_end = fwd[-1] * loss_sec
            b_start = bwd[0] * loss_sec

            fwd[0] = inp + r_g * b_start
            fwd[1:] = nf_inner
            bwd[:-1] = nb_inner
            lp = one_m_lam * f_end + lam * lp
            bwd[-1] = r_l0 * lp
            out[k] = f_end - bwd[-1]
        return out, ug_out

    lb = wg["loss_biquad"]
    b0, b1, b2 = lb["b0"], lb["b1"], lb["b2"]
    a1, a2 = float(lb["a1"]), float(lb["a2"])
    waves = np.zeros((2, n_sec), dtype=np.float64)
    z1 = np.zeros((2, n_sec), dtype=np.float64)
    z2 = np.zeros((2, n_sec), dtype=np.float64)
    fwd_v = waves[0]
    bwd_v = waves[1]

    b_rad, a_rad = wg["b_rad"], wg["a_rad"]
    br0, br1, br2 = float(b_rad[0]), float(b_rad[1]), float(b_rad[2])
    ar1, ar2 = float(a_rad[1]), float(a_rad[2])
    zr1, zr2 = 0.0, 0.0

    wb = wg.get("wall_branch")
    if wb is not None:
        alpha1, alpha2, alphaw = wb["alpha1"], wb["alpha2"], wb["alphaw"]
        g1, g2, g3 = wb["g1"], wb["g2"], wb["g3"]
        w_in_1 = np.zeros(n_sec - 1, dtype=np.float64)
        w_in_2 = np.zeros(n_sec - 1, dtype=np.float64)
        w_out_1 = np.zeros(n_sec - 1, dtype=np.float64)
        w_out_2 = np.zeros(n_sec - 1, dtype=np.float64)

    for k in range(n_steps):
        if impulse:
            inp = 1.0 if k == 0 else 0.0
        elif glottal is not None:
            p_supra = z0 * (fwd_v[0] + bwd_v[0])
            lim = 0.35 * p_sub_max
            if p_supra > lim:
                p_supra = lim
            elif p_supra < -lim:
                p_supra = -lim
            p_sub = float(p_sub_env[k])
            dp = p_sub - p_supra if p_sub > p_supra else 0.0
            a1_g = 2.0 * l_g * (x0 + x1)
            a2_g = 2.0 * l_g * (x0 + x2)
            a_min = a1_g if a1_g < a2_g else a2_g
            if a_min > 0.0 and a1_g > 0.0:
                ug = a_min * ((2.0 * dp / rho) ** 0.5)
                p1_g = p_supra + dp * (1.0 - (a_min / a1_g) ** 2) if a1_g > a_min else p_supra
            else:
                ug = 0.0
                p1_g = p_sub if a1_g > 0.0 else 0.0
            f1_force = p1_g * l_g * d1
            fc1 = -c1 * (x0 + x1) if (x0 + x1) < 0.0 else 0.0
            fc2 = -c2 * (x0 + x2) if (x0 + x2) < 0.0 else 0.0
            acc1 = (-r1 * v1 - k1 * x1 - kc * (x1 - x2) + fc1 + f1_force) / m1
            acc2 = (-r2 * v2 - k2 * x2 - kc * (x2 - x1) + fc2) / m2
            v1 += acc1 * dt
            x1 += v1 * dt
            v2 += acc2 * dt
            x2 += v2 * dt
            ug_out[k] = ug
            inp = ug
        else:
            inp = 0.0

        # TDF-II section Kirchhoff loss filters for forward & backward waves in one batch
        filt = waves * b0 + z1
        z1 = waves * b1 - filt * a1 + z2
        z2 = waves * b2 - filt * a2
        ff_v = filt[0]
        bf_v = filt[1]

        f_in = ff_v[:-1]
        b_in = bf_v[1:]
        if wb is None:
            nf_inner = one_m_rj * f_in - r_junc * b_in
            nb_inner = r_junc * f_in + one_p_rj * b_in
        else:
            w_in = g1 * (w_out_1 - w_in_1) + g2 * w_out_2 + g3 * w_in_2
            pj = alpha1 * f_in + alpha2 * b_in + alphaw * w_in
            nb_inner = pj - f_in
            nf_inner = pj - b_in
            w_out = pj - w_in
            w_in_2 = w_in_1
            w_in_1 = w_in
            w_out_2 = w_out_1
            w_out_1 = w_out
        f_end = float(ff_v[-1])
        b_start = float(bf_v[0])

        b_ref = br0 * f_end + zr1
        zr1 = br1 * f_end - ar1 * b_ref + zr2
        zr2 = br2 * f_end - ar2 * b_ref

        fwd_v[0] = inp + r_g * b_start
        fwd_v[1:] = nf_inner
        bwd_v[:-1] = nb_inner
        bwd_v[-1] = b_ref
        out[k] = f_end - b_ref

    return out, ug_out


def kl_impulse_response(wg: dict, *, ir_duration_s: float = 0.8) -> np.ndarray:
    """Compute Kelly–Lochbaum volume-velocity impulse response (closed glottis, impulse input)."""
    n_ir = max(256, int(round(ir_duration_s * wg["fs_wg"])))
    ir, _ = _run_kl_loop(wg, n_ir, impulse=True)
    if wg.get("wall_branch") is not None and len(ir) >= 2000:
        head_e = float(np.sum(ir[:1000] ** 2))
        tail_e = float(np.sum(ir[-1000:] ** 2))
        if head_e > 0.0 and tail_e > 0.05 * head_e:
            raise FloatingPointError(
                f"KL 3-port wall junction IR failed passivity energy decay check (tail/head={tail_e / head_e:.3e})."
            )
    return ir


def extract_kl_peaks(
    ir: np.ndarray,
    fs_wg: float,
    *,
    fmax_hz: float = 3500.0,
    zero_pad_s: float | None = None,
    interpolate_bw: bool = False,
) -> tuple[np.ndarray, np.ndarray, list[dict]]:
    n_fft = len(ir) if zero_pad_s is None else max(len(ir), int(round(float(zero_pad_s) * fs_wg)))
    spec = np.abs(np.fft.rfft(ir, n=n_fft))
    freqs = np.fft.rfftfreq(n_fft, 1.0 / fs_wg)
    mag_db = 20.0 * np.log10(np.maximum(spec, 1e-12))
    mask = (freqs >= 30.0) & (freqs <= fmax_hz)
    f_band = freqs[mask]
    m_band = mag_db[mask]
    pks, _ = signal.find_peaks(m_band, prominence=2.0)

    out: list[dict] = []
    df = float(freqs[1] - freqs[0])
    for idx in pks:
        if 1 <= idx < len(f_band) - 1:
            y0, y1, y2 = m_band[idx - 1], m_band[idx], m_band[idx + 1]
            denom = y0 - 2.0 * y1 + y2
            delta = float(np.clip(0.5 * (y0 - y2) / denom, -0.5, 0.5)) if abs(denom) > 1e-12 else 0.0
            f_pk = float(f_band[idx] + delta * df)
            m_pk = float(y1 - 0.25 * (y0 - y2) * delta)
        else:
            f_pk = float(f_band[idx])
            m_pk = float(m_band[idx])

        target = m_pk - 3.0
        l_idx, r_idx = int(idx), int(idx)
        while l_idx > 0 and m_band[l_idx] > target:
            l_idx -= 1
        while r_idx < len(m_band) - 1 and m_band[r_idx] > target:
            r_idx += 1

        if interpolate_bw:
            if l_idx < idx and m_band[l_idx] <= target < m_band[l_idx + 1]:
                tl = (target - m_band[l_idx]) / max(m_band[l_idx + 1] - m_band[l_idx], 1e-12)
                f_l = float(f_band[l_idx] + tl * (f_band[l_idx + 1] - f_band[l_idx]))
            else:
                f_l = float(f_band[l_idx])
            if r_idx > idx and m_band[r_idx] <= target < m_band[r_idx - 1]:
                tr = (target - m_band[r_idx]) / max(m_band[r_idx - 1] - m_band[r_idx], 1e-12)
                f_r = float(f_band[r_idx] - tr * (f_band[r_idx] - f_band[r_idx - 1]))
            else:
                f_r = float(f_band[r_idx])
            bw = max(float(f_r - f_l), df)
        else:
            bw = max(float(f_band[r_idx] - f_band[l_idx]), df)

        out.append({"f_hz": round(f_pk, 3), "bw_3db_hz": round(bw, 3), "q": round(f_pk / bw, 3), "mag_db": round(m_pk, 3)})
    return freqs, mag_db, out


def compare_kl_vs_tmm(
    wg: dict,
    ir: np.ndarray,
    tmm_summary: dict,
    *,
    max_freq_hz: float = 2000.0,
    max_rel_err: float = 0.02,
    out_png: Path | None = None,
    tmm_raw: dict | None = None,
) -> dict:
    """Verify KL impulse-response resonances match M2 TMM peaks (F1-F3) in frequency and log Q error."""
    is_legacy = wg.get("kl_loss_model", "legacy") == "legacy"
    freqs_kl, mag_kl_db, kl_peaks_raw = extract_kl_peaks(
        ir,
        wg["fs_wg"],
        fmax_hz=3500.0,
        zero_pad_s=None if is_legacy else 4.0,
        interpolate_bw=not is_legacy,
    )
    wb = wg.get("wall_branch")
    if wb is not None and wb.get("has_stiffness"):
        kl_peaks = [kp for kp in kl_peaks_raw if not (kp["f_hz"] < 120.0 and kp["q"] < 15.0)]
        if not kl_peaks:
            kl_peaks = kl_peaks_raw
    else:
        kl_peaks = kl_peaks_raw
    tmm_sub = [pk for pk in tmm_summary["peaks"] if pk["f_hz"] <= max_freq_hz]
    if len(tmm_sub) < 3 and len(tmm_summary["peaks"]) >= 3:
        tmm_sub = list(tmm_summary["peaks"][:3])

    comparisons: list[dict] = []
    rel_errors: list[float] = []
    q_rel_errors: list[float] = []
    for idx, t_pk in enumerate(tmm_sub):
        if not kl_peaks:
            raise AssertionError("KL waveguide impulse response produced no resonance peaks.")
        nearest = min(kl_peaks, key=lambda kp: abs(kp["f_hz"] - t_pk["f_hz"]))
        rel_err = abs(nearest["f_hz"] - t_pk["f_hz"]) / max(t_pk["f_hz"], 1e-9)
        q_rel_err = abs(nearest["q"] - t_pk["q"]) / max(t_pk["q"], 1e-9)
        rel_errors.append(round(float(rel_err), 6))
        q_rel_errors.append(round(float(q_rel_err), 6))
        comparisons.append(
            {
                "formant": f"F{idx + 1}",
                "tmm_hz": t_pk["f_hz"],
                "kl_hz": nearest["f_hz"],
                "rel_err": round(float(rel_err), 6),
                "tmm_q": t_pk["q"],
                "kl_q": nearest["q"],
                "q_rel_err": round(float(q_rel_err), 6),
            }
        )

    passed = bool(comparisons) and all(err <= max_rel_err for err in rel_errors)
    if not passed:
        raise AssertionError(
            f"KL waveguide resonances deviate > {max_rel_err*100:.1f}% from M2 TMM below {max_freq_hz:.0f} Hz: {comparisons}"
        )

    if out_png is not None:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(10, 4.5))
        mask_kl = (freqs_kl >= 20.0) & (freqs_kl <= 3500.0)
        kl_norm_db = mag_kl_db[mask_kl] - np.max(mag_kl_db[mask_kl])
        ax.plot(freqs_kl[mask_kl], kl_norm_db, color="darkorange", lw=1.6, label="Guide d'onde KL |H_ir(f)| (dB norm.)")
        if tmm_raw is not None:
            f_t = tmm_raw["freqs_hz"]
            m_t = 20.0 * np.log10(np.maximum(np.abs(tmm_raw["h_transfer"]), 1e-12))
            mt_mask = f_t <= 3500.0
            ax.plot(f_t[mt_mask], m_t[mt_mask] - np.max(m_t[mt_mask]), color="navy", lw=1.3, ls="--", label="TMM avec pertes |H(f)| (dB norm.)")
        for c_item in comparisons:
            ax.axvline(c_item["tmm_hz"], color="teal", ls=":", lw=1.0, alpha=0.7)
        ax.set_xlabel("Fréquence (Hz)")
        ax.set_ylabel("Magnitude normalisée (dB)")
        max_q_err = max(q_rel_errors) * 100.0 if q_rel_errors else 0.0
        ax.set_title(
            f"Comparaison KL ({wg.get('kl_loss_model', 'legacy')}) vs TMM — "
            f"err max F={max(rel_errors)*100:.2f} %, Q={max_q_err:.1f} %"
        )
        ax.grid(True, alpha=0.3)
        ax.legend(loc="upper right", fontsize=8)
        fig.tight_layout()
        fig.savefig(out_png, dpi=96)
        plt.close(fig)

    return {
        "passed": passed,
        "kl_peaks": kl_peaks,
        "comparisons": comparisons,
        "rel_errors": rel_errors,
        "q_rel_errors": q_rel_errors,
        "kl_vs_tmm_q_rel_err": q_rel_errors,
        "loss_model": wg.get("kl_loss_model", "legacy"),
        "fit_rms_db": wg.get("fit_rms_db"),
        "notes": list(wg.get("radiation_notes", [])),
    }


def _interactive_two_mass_flow(
    wg: dict,
    air: Air,
    *,
    duration_s: float,
    f0_prior_hz: float,
    seed: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Run SI two-mass glottal oscillator (Steinecke & Herzel 1995) coupled with KL supraglottal pressure."""
    fs_wg = wg["fs_wg"]
    dt = 1.0 / fs_wg
    n_wg = int(round(duration_s * fs_wg))
    z0 = wg["Z0_glottis"]

    # Biomechanical SI two-mass parameters calibrated around f0_prior_hz (baseline 65.2 Hz at k1=210 N/m)
    m1, m2 = 1.5e-3, 0.4e-3
    scale_k = (max(f0_prior_hz, 20.0) / 65.2) ** 2
    k1, k2, kc = 210.0 * scale_k, 25.0 * scale_k, 40.0 * scale_k
    r1 = 2.0 * 0.08 * np.sqrt(m1 * k1)
    r2 = 2.0 * 0.12 * np.sqrt(m2 * k2)
    c1, c2 = 3.0 * k1, 3.0 * k2
    l_g, d1, x0 = 0.025, 0.005, 0.0002
    p_sub_max = 1200.0

    t_wg = np.arange(n_wg) * dt
    p_sub_env = p_sub_max * np.interp(t_wg, [0.0, 0.08, duration_s * 0.85, duration_s], [0.0, 1.0, 0.9, 0.0])

    if wg.get("kl_loss_model", "legacy") == "legacy":
        glottal = {
            "p_sub_env": p_sub_env,
            "p_sub_max": p_sub_max,
            "l_g": l_g,
            "d1": d1,
            "x0": x0,
            "m1": m1,
            "m2": m2,
            "r1": r1,
            "r2": r2,
            "k1": k1,
            "k2": k2,
            "kc": kc,
            "c1": c1,
            "c2": c2,
            "dt": dt,
            "rho": air.rho,
        }
        _, ug_wg = _run_kl_loop(wg, n_wg, glottal=glottal)
    else:
        p_sub_list = p_sub_env.tolist()
        ug_list = [0.0] * n_wg
        lim = 0.35 * p_sub_max
        inv_rho_2 = 2.0 / air.rho
        rj0 = float(wg["r_junc"][0]) if len(wg["r_junc"]) > 0 else 0.0
        ls = float(wg["loss_sec"])
        r_g = float(wg["R_g"])
        x1, v1, x2, v2 = 0.0001, 0.0, 0.0, 0.0
        fw0, bw0 = 0.0, 0.0
        for k in range(n_wg):
            p_supra = z0 * (fw0 + bw0)
            if p_supra > lim:
                p_supra = lim
            elif p_supra < -lim:
                p_supra = -lim
            p_sub = p_sub_list[k]
            dp = p_sub - p_supra if p_sub > p_supra else 0.0
            a1_g = 2.0 * l_g * (x0 + x1)
            a2_g = 2.0 * l_g * (x0 + x2)
            a_min = a1_g if a1_g < a2_g else a2_g
            if a_min > 0.0 and a1_g > 0.0:
                ug = a_min * ((dp * inv_rho_2) ** 0.5)
                p1_g = p_supra + dp * (1.0 - (a_min / a1_g) ** 2) if a1_g > a_min else p_supra
            else:
                ug = 0.0
                p1_g = p_sub if a1_g > 0.0 else 0.0
            f1_force = p1_g * l_g * d1
            fc1 = -c1 * (x0 + x1) if (x0 + x1) < 0.0 else 0.0
            fc2 = -c2 * (x0 + x2) if (x0 + x2) < 0.0 else 0.0
            v1 += ((-r1 * v1 - k1 * x1 - kc * (x1 - x2) + fc1 + f1_force) / m1) * dt
            x1 += v1 * dt
            v2 += ((-r2 * v2 - k2 * x2 - kc * (x2 - x1) + fc2) / m2) * dt
            x2 += v2 * dt
            ug_list[k] = ug
            bw0 = rj0 * (fw0 * ls)
            fw0 = ug + r_g * (bw0 * ls)
        ug_wg = np.asarray(ug_list, dtype=np.float64)

    # Physiological glottal flow derivative (return-phase lowpass ~ 240 Hz gives ~ -10 dB/oct tilt)
    rng = np.random.default_rng(seed)
    dug_wg = np.gradient(ug_wg, dt)
    b_ret, a_ret = signal.butter(1, 240.0 / (fs_wg / 2.0), btype="low")
    dug_ret = signal.lfilter(b_ret, a_ret, dug_wg)

    # Controlled broadband aspiration noise (up to 4 kHz) modulated by respiratory envelope
    b_as, a_as = signal.butter(2, 4000.0 / (fs_wg / 2.0), btype="low")
    asp = signal.lfilter(b_as, a_as, rng.normal(0.0, 1.0, size=n_wg))
    asp = asp / (np.std(asp) + 1e-12) * (0.09 * np.std(dug_ret)) * (p_sub_env / p_sub_max)
    src_wg = dug_ret + asp
    return src_wg, ug_wg


def synthesize_kl(
    s_mm: np.ndarray,
    area_mm2: np.ndarray,
    air: Air,
    tmm_summary: dict,
    *,
    dh_mm: np.ndarray | None = None,
    perimeter_mm: np.ndarray | None = None,
    duration_s: float = 2.5,
    sr: int = 44100,
    f0_hz: float = 60.0,
    radiation: str = "unflanged",
    seed: int = 0,
    kl_plot_path: Path | None = None,
    tmm_raw: dict | None = None,
    kl_loss_model: str = "kirchhoff_fit",
    extra_damping_factor: float = 1.0,
    wall: WallModel | None = None,
) -> dict:
    """Synthesize source.wav (glottal flow derivative) and filtered.wav (radiated KL output)."""
    tmm_peaks = tmm_summary.get("peaks") or []
    f1_hint = float(tmm_peaks[0]["f_hz"]) if tmm_peaks else 400.0
    wg = build_kl_waveguide(
        s_mm,
        area_mm2,
        air,
        dh_mm=dh_mm,
        perimeter_mm=perimeter_mm,
        radiation=radiation,
        sr=sr,
        oversample=2,
        kl_loss_model=kl_loss_model,
        extra_damping_factor=extra_damping_factor,
        wall=wall,
        f1_hint_hz=f1_hint,
    )
    ir_dur = 0.8 if kl_loss_model == "legacy" else 0.45
    ir_wg = kl_impulse_response(wg, ir_duration_s=ir_dur)
    kl_cmp = compare_kl_vs_tmm(wg, ir_wg, tmm_summary, out_png=kl_plot_path, tmm_raw=tmm_raw)

    src_wg, _ = _interactive_two_mass_flow(wg, air, duration_s=duration_s, f0_prior_hz=f0_hz, seed=seed)
    n_wg = len(src_wg)
    flt_wg = signal.fftconvolve(src_wg, ir_wg, mode="full")[:n_wg]

    n_sr = int(round(duration_s * sr))
    source_raw = signal.resample(src_wg, n_sr)
    filtered_raw = signal.resample(flt_wg, n_sr)

    source, src_peak_dbfs, src_clip = normalize_audio(source_raw, target_dbfs=-1.0)
    filtered, flt_peak_dbfs, flt_clip = normalize_audio(filtered_raw, target_dbfs=-1.0)

    return {
        "source": source,
        "filtered": filtered,
        "source_description": "interactive_two_mass_glottal_flow_derivative_with_aspiration",
        "filtered_description": (
            "kelly_lochbaum_radiated_pressure_yielding_walls"
            if wg.get("wall_branch") is not None
            else "kelly_lochbaum_radiated_pressure"
        ),
        "independent_of_tmm": True,
        "kl_check": kl_cmp,
        "waveguide": {
            "N_sections": wg["N"],
            "fs_wg": round(wg["fs_wg"], 2),
            "dx_mm": round(wg["dx_m"] * 1e3, 4),
            "L_wg_m": round(wg["L_wg_m"], 6),
            "kl_loss_model": wg["kl_loss_model"],
            "fit_rms_db": wg["fit_rms_db"],
            "has_yielding_walls": wg.get("wall_branch") is not None,
        },
        "audio": {
            "source_peak_dbfs": src_peak_dbfs,
            "filtered_peak_dbfs": flt_peak_dbfs,
            "clipping_detected": bool(src_clip or flt_clip),
        },
    }


def synthesize_tmm_fir(
    tmm_raw: dict,
    air: Air,
    *,
    duration_s: float = 2.5,
    sr: int = 44100,
    f0_hz: float = 60.0,
    seed: int = 0,
    ir_len_s: float = 0.5,
) -> dict:
    """Non-independent sensitivity audio fallback: two-mass source convolved with min-phase FIR of H_TMM.

    Marked independent_of_tmm=False; peer_review refuses it for independent validation.
    """
    n_sr = int(round(duration_s * sr))
    dt = 1.0 / sr
    # Non-interactive two-mass glottal source at sample rate sr
    m1, m2 = 1.5e-3, 0.4e-3
    scale_k = (max(f0_hz, 20.0) / 65.2) ** 2
    k1, k2, kc = 210.0 * scale_k, 25.0 * scale_k, 40.0 * scale_k
    r1 = 2.0 * 0.08 * np.sqrt(m1 * k1)
    r2 = 2.0 * 0.12 * np.sqrt(m2 * k2)
    c1, c2 = 3.0 * k1, 3.0 * k2
    l_g, d1, x0 = 0.025, 0.005, 0.0002
    p_sub_max = 1200.0
    t_arr = np.arange(n_sr) * dt
    p_sub_env = p_sub_max * np.interp(t_arr, [0.0, 0.08, duration_s * 0.85, duration_s], [0.0, 1.0, 0.9, 0.0])

    x1, v1, x2, v2 = 0.0001, 0.0, 0.0, 0.0
    ug = np.zeros(n_sr, dtype=np.float64)
    for i in range(n_sr):
        dp = float(p_sub_env[i])
        a1_g = 2.0 * l_g * (x0 + x1)
        a2_g = 2.0 * l_g * (x0 + x2)
        a_min = min(a1_g, a2_g)
        if a_min > 0.0 and a1_g > 0.0:
            flow = a_min * np.sqrt(2.0 * dp / air.rho)
            p1_g = dp * (1.0 - (a_min / a1_g) ** 2) if a1_g > a_min else 0.0
        else:
            flow = 0.0
            p1_g = dp if a1_g > 0.0 else 0.0
        fc1 = -c1 * (x0 + x1) if (x0 + x1) < 0.0 else 0.0
        fc2 = -c2 * (x0 + x2) if (x0 + x2) < 0.0 else 0.0
        acc1 = (-r1 * v1 - k1 * x1 - kc * (x1 - x2) + fc1 + p1_g * l_g * d1) / m1
        acc2 = (-r2 * v2 - k2 * x2 - kc * (x2 - x1) + fc2) / m2
        v1 += acc1 * dt
        x1 += v1 * dt
        v2 += acc2 * dt
        x2 += v2 * dt
        ug[i] = flow

    rng = np.random.default_rng(seed)
    dug = np.gradient(ug, dt)
    b_ret, a_ret = signal.butter(1, 240.0 / (sr / 2.0), btype="low")
    dug_ret = signal.lfilter(b_ret, a_ret, dug)
    b_as, a_as = signal.butter(2, 4000.0 / (sr / 2.0), btype="low")
    asp = signal.lfilter(b_as, a_as, rng.normal(0.0, 1.0, size=n_sr))
    asp = asp / (np.std(asp) + 1e-12) * (0.09 * np.std(dug_ret)) * (p_sub_env / p_sub_max)
    src_raw = dug_ret + asp

    # Minimum-phase FIR via real cepstrum of |H_TMM(f)|
    n_fft = 1 << int(np.ceil(np.log2(max(2048, int(round(ir_len_s * sr))))))
    f_fft = np.fft.rfftfreq(n_fft, 1.0 / sr)
    f_tmm = np.asarray(tmm_raw["freqs_hz"], dtype=np.float64)
    mag_tmm = np.abs(np.asarray(tmm_raw["h_transfer"], dtype=np.complex128))
    mag_interp = np.interp(f_fft, f_tmm, mag_tmm, left=float(mag_tmm[0]), right=float(mag_tmm[-1]) * 0.1)
    log_mag = np.log(np.maximum(mag_interp, 1e-12))
    sym_log = np.concatenate([log_mag, log_mag[-2:0:-1]])
    cep = np.fft.ifft(sym_log).real
    win = np.zeros(n_fft, dtype=np.float64)
    win[0] = 1.0
    win[1 : n_fft // 2] = 2.0
    win[n_fft // 2] = 1.0
    h_min = np.fft.ifft(np.exp(np.fft.fft(cep * win))).real
    ir_fir = h_min[: int(round(ir_len_s * sr))]
    flt_raw = signal.fftconvolve(src_raw, ir_fir, mode="full")[:n_sr]

    source, src_peak_dbfs, src_clip = normalize_audio(src_raw, target_dbfs=-1.0)
    filtered, flt_peak_dbfs, flt_clip = normalize_audio(flt_raw, target_dbfs=-1.0)

    return {
        "synth": "tmm_fir",
        "independent_of_tmm": False,
        "source": source,
        "filtered": filtered,
        "source_description": "non_interactive_two_mass_glottal_flow_derivative",
        "filtered_description": "tmm_minimum_phase_fir_non_independent",
        "kl_check": None,
        "audio": {
            "source_peak_dbfs": src_peak_dbfs,
            "filtered_peak_dbfs": flt_peak_dbfs,
            "clipping_detected": bool(src_clip or flt_clip),
        },
    }


def synthesize_legacy_iir(
    formants_hz: list[float],
    q_factors: list[float],
    *,
    duration_s: float = 2.5,
    sr: int = 44100,
    f0_hz: float = 60.0,
    seed: int = 0,
) -> dict:
    """Legacy baseline source + IIR peak filter cascade kept for honest regression comparison."""
    rng = np.random.default_rng(seed)
    n = int(round(duration_s * sr))
    t = np.arange(n) / sr
    phase = 2.0 * np.pi * f0_hz * t + 0.15 * np.sin(2.0 * np.pi * 3.0 * t)
    pulse = np.sin(phase) + 0.35 * np.sin(2.0 * phase)
    env = np.interp(t, [0.0, 0.08, duration_s * 0.7, duration_s], [0.0, 1.0, 0.85, 0.0])
    turb = rng.normal(0.0, 1.0, size=n)
    b_lp, a_lp = signal.butter(2, min(0.45, 1000.0 / (sr / 2.0)), btype="low")
    turb = signal.lfilter(b_lp, a_lp, turb)
    src_raw = (pulse + 0.03 * turb) * env

    out_raw = np.zeros_like(src_raw)
    gains = [1.0, 0.45, 0.2][: len(formants_hz)]
    for f, q, g in zip(formants_hz, q_factors, gains):
        if 0.0 < f < sr / 2.0:
            b, a = signal.iirpeak(f, Q=max(0.5, q), fs=sr)
            out_raw += signal.lfilter(b, a, src_raw) * g

    source, src_peak_dbfs, src_clip = normalize_audio(src_raw, target_dbfs=-1.0)
    filtered, flt_peak_dbfs, flt_clip = normalize_audio(out_raw, target_dbfs=-1.0)

    return {
        "source": source,
        "filtered": filtered,
        "source_description": "legacy_two_sine_sketch_plus_turbulence",
        "filtered_description": "legacy_parallel_iir_peak_filter",
        "audio": {
            "source_peak_dbfs": src_peak_dbfs,
            "filtered_peak_dbfs": flt_peak_dbfs,
            "clipping_detected": bool(src_clip or flt_clip),
        },
    }

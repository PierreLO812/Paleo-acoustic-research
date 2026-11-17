"""M2 — Independent lossy Transfer Matrix Method (TMM) acoustic solver.

Implements cylindrical/conical 1D waveguide propagation with Kirchhoff
thermoviscous wall losses (Kirchhoff 1868; Keefe 1984, J. Acoust. Soc. Am. 75(1))
and Levine–Schwinger (1948) / Silva et al. (2009) unflanged (or flanged) outlet
radiation impedance.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
from scipy import signal

from pipeline.geometry_qc import load_area_profile


@dataclass(frozen=True)
class WallModel:
    """Lumped yielding-wall mechanical impedance Z_w(w) = R_w + j*w*M_w + K_w/(j*w) (Sondhi 1974; Ishizaka et al. 1975)."""

    mass_per_area: float = 0.0
    stiffness_per_area: float = 0.0
    resistance_per_area: float = 0.0

    def is_rigid(self) -> bool:
        return (
            self.mass_per_area == 0.0
            and self.stiffness_per_area == 0.0
            and self.resistance_per_area == 0.0
        )

    def impedance(self, omega: np.ndarray) -> np.ndarray:
        """Specific wall impedance Z_w(omega) in Pa*s/m."""
        return (
            float(self.resistance_per_area)
            + 1j * omega * float(self.mass_per_area)
            + float(self.stiffness_per_area) / (1j * omega)
        )


@dataclass(frozen=True)
class Air:
    """Thermodynamic and transport properties of air derived from temperature."""

    T_c: float = 20.0
    c: float = 343.0
    rho: float = 1.204
    mu: float = 1.81e-5
    gamma: float = 1.4
    Pr: float = 0.71

    @classmethod
    def from_temperature(cls, T_c: float = 20.0) -> "Air":
        if abs(T_c - 20.0) < 1e-6:
            return cls(T_c=20.0, c=343.0, rho=1.204, mu=1.81e-5, gamma=1.4, Pr=0.71)
        t_k = float(T_c) + 273.15
        c = 331.3 * np.sqrt(t_k / 273.15)
        rho = 101325.0 / (287.05 * t_k)
        # Sutherland's formula for dynamic viscosity of air
        t0 = 291.15
        mu0 = 1.827e-5
        s_const = 120.0
        mu = mu0 * ((t0 + s_const) / (t_k + s_const)) * (t_k / t0) ** 1.5
        return cls(
            T_c=float(T_c),
            c=round(float(c), 4),
            rho=round(float(rho), 6),
            mu=float(mu),
            gamma=1.4,
            Pr=0.71,
        )


def _subdivide_profile(
    s_m: np.ndarray,
    area_m2: np.ndarray,
    dh_m: np.ndarray,
    max_seg_len_m: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Subdivide geometry intervals so every segment satisfies ds <= lambda_min / 20."""
    lengths: list[float] = []
    areas: list[float] = []
    radii: list[float] = []

    for i in range(len(s_m) - 1):
        ds = float(s_m[i + 1] - s_m[i])
        n_sub = max(1, int(np.ceil(ds / max_seg_len_m)))
        sub_len = ds / n_sub
        a0, a1 = float(area_m2[i]), float(area_m2[i + 1])
        r0, r1 = 0.5 * float(dh_m[i]), 0.5 * float(dh_m[i + 1])
        for k in range(n_sub):
            t_mid = (k + 0.5) / n_sub
            lengths.append(sub_len)
            areas.append(a0 + t_mid * (a1 - a0))
            radii.append(max(r0 + t_mid * (r1 - r0), 1e-4))

    return np.asarray(lengths), np.asarray(areas), np.asarray(radii)


def radiation_impedance(
    omega: np.ndarray,
    a_out_m: float,
    zc_out: float,
    air: Air,
    model: str = "unflanged",
) -> tuple[np.ndarray, list[str]]:
    """Compute outlet radiation impedance Zr(omega) using Levine–Schwinger / Silva et al. (2009)."""
    warnings: list[str] = []
    ka = (omega / air.c) * a_out_m

    if model == "ideal_open":
        return np.zeros_like(omega, dtype=np.complex128), warnings

    if model == "flanged":
        if np.any(ka > 1.0):
            warnings.append(
                f"ka > 1 inside band of interest (max ka={float(np.max(ka)):.2f}); "
                "using flanged high-ka saturation."
            )
        low = ka < 1.0
        z_norm = np.empty_like(ka, dtype=np.complex128)
        z_norm[low] = 0.5 * (ka[low] ** 2) + 1j * 0.8216 * ka[low]
        # Passive 2nd-order rational continuation for ka >= 1
        s = 1j * ka[~low]
        z_norm[~low] = (0.8216 * s + 0.50 * (s ** 2)) / (1.0 + 1.217 * s + 0.50 * (s ** 2))
        return zc_out * z_norm, warnings

    if model != "unflanged":
        raise ValueError(f"Unsupported radiation model '{model}'. Expected 'unflanged', 'flanged', or 'ideal_open'.")

    z_norm = np.empty_like(ka, dtype=np.complex128)
    low = ka < 1.0
    # Levine & Schwinger (1948) low-ka approximation
    z_norm[low] = 0.25 * (ka[low] ** 2) + 1j * 0.6133 * ka[low]

    if np.any(~low):
        warnings.append(
            f"ka > 1 inside band of interest (max ka={float(np.max(ka)):.2f}); "
            "switched to Silva et al. (2009) rational fit above ka=1."
        )
        # Silva et al. (2009) passive 2nd-order rational fit for unflanged circular pipe
        # Matches 0.6133*j*ka + 0.25*(ka)^2 at low ka and converges to 1.0 at high ka.
        s = 1j * ka[~low]
        z_norm[~low] = (0.6133 * s + 0.2406 * (s ** 2)) / (1.0 + 0.80 * s + 0.2406 * (s ** 2))

    return zc_out * z_norm, warnings


def _parabolic_peak(freqs: np.ndarray, mag_db: np.ndarray, idx: int) -> tuple[float, float]:
    if 1 <= idx < len(freqs) - 1:
        y0, y1, y2 = mag_db[idx - 1], mag_db[idx], mag_db[idx + 1]
        denom = y0 - 2.0 * y1 + y2
        if abs(denom) > 1e-12:
            delta = float(np.clip(0.5 * (y0 - y2) / denom, -0.5, 0.5))
            df = float(freqs[idx + 1] - freqs[idx])
            f_pk = float(freqs[idx] + delta * df)
            m_pk = float(y1 - 0.25 * (y0 - y2) * delta)
            return f_pk, m_pk
    return float(freqs[idx]), float(mag_db[idx])


def _bw_3db(freqs: np.ndarray, mag_db: np.ndarray, idx: int, peak_db: float) -> float:
    target = peak_db - 3.0
    df = float(freqs[1] - freqs[0]) if len(freqs) > 1 else 0.5

    left = idx
    while left > 0 and mag_db[left] > target:
        left -= 1
    if left < idx and mag_db[left] <= target < mag_db[left + 1]:
        t = (target - mag_db[left]) / max(mag_db[left + 1] - mag_db[left], 1e-12)
        f_l = float(freqs[left] + t * (freqs[left + 1] - freqs[left]))
    else:
        f_l = float(freqs[left])

    right = idx
    n = len(freqs)
    while right < n - 1 and mag_db[right] > target:
        right += 1
    if right > idx and mag_db[right] <= target < mag_db[right - 1]:
        t = (target - mag_db[right]) / max(mag_db[right - 1] - mag_db[right], 1e-12)
        f_r = float(freqs[right] - t * (freqs[right] - freqs[right - 1]))
    else:
        f_r = float(freqs[right])

    return float(max(f_r - f_l, df))


def _extract_peaks(freqs: np.ndarray, mag: np.ndarray, *, min_prominence_db: float = 1.0) -> list[dict]:
    mag_db = 20.0 * np.log10(np.maximum(mag, 1e-15))
    idx_peaks, _ = signal.find_peaks(mag_db, prominence=min_prominence_db)
    peaks: list[dict] = []
    for idx in idx_peaks:
        f_pk, m_pk = _parabolic_peak(freqs, mag_db, int(idx))
        bw = _bw_3db(freqs, mag_db, int(idx), m_pk)
        q = float(f_pk / bw) if bw > 0 else 0.0
        peaks.append(
            {
                "f_hz": round(f_pk, 3),
                "bw_3db_hz": round(bw, 3),
                "q": round(q, 3),
                "mag_db": round(m_pk, 3),
            }
        )
    return peaks


def solve_tmm(
    s_mm: np.ndarray,
    area_mm2: np.ndarray,
    dh_mm: np.ndarray,
    *,
    air: Air | None = None,
    fmin: float = 10.0,
    fmax: float = 4000.0,
    df: float = 0.5,
    radiation: str = "unflanged",
    losses: bool = True,
    wall: WallModel | None = None,
    extra_damping_factor: float = 1.0,
) -> dict:
    """Compute Z_in(f), H(f) = U_out/U_in (or P_out/U_in if radiation='closed'), and resonance peaks."""
    air = air or Air()
    if df > 0.5:
        raise ValueError(f"Frequency step df={df} Hz exceeds required <= 0.5 Hz.")

    s_m = np.asarray(s_mm, dtype=np.float64) * 1e-3
    area_m2 = np.asarray(area_mm2, dtype=np.float64) * 1e-6
    dh_m = np.asarray(dh_mm, dtype=np.float64) * 1e-3

    max_seg_len = air.c / (20.0 * max(fmax, 1.0))
    seg_len, seg_area, seg_rad = _subdivide_profile(s_m, area_m2, dh_m, max_seg_len)

    freqs = np.arange(fmin, fmax + 0.5 * df, df, dtype=np.float64)
    omega = 2.0 * np.pi * freqs

    a_mat = np.ones_like(freqs, dtype=np.complex128)
    b_mat = np.zeros_like(freqs, dtype=np.complex128)
    c_mat = np.zeros_like(freqs, dtype=np.complex128)
    d_mat = np.ones_like(freqs, dtype=np.complex128)

    if losses:
        wall_factor = float(extra_damping_factor) * np.sqrt(omega * air.mu / (2.0 * air.rho)) * (
            1.0 + (air.gamma - 1.0) / np.sqrt(air.Pr)
        )
    else:
        wall_factor = np.zeros_like(omega)

    use_wall = wall is not None and not wall.is_rigid()
    inv_z_w = (1.0 / wall.impedance(omega)) if use_wall else None

    for lj, aj, rj in zip(seg_len, seg_area, seg_rad):
        alpha = wall_factor / (rj * air.c)
        c_eff = air.c * np.maximum(0.5, 1.0 - (alpha * air.c) / omega)
        # Complex wavenumber k = omega/c_eff - j*alpha  =>  propagation gamma_p = j*k = alpha + j*omega/c_eff
        k_complex = (omega / c_eff) - 1j * alpha
        gamma_p = 1j * k_complex
        zc = air.rho * air.c / aj

        ch = np.cosh(gamma_p * lj)
        sh = np.sinh(gamma_p * lj)

        if not use_wall:
            na = a_mat * ch + b_mat * (sh / zc)
            nb = a_mat * (zc * sh) + b_mat * ch
            nc = c_mat * ch + d_mat * (sh / zc)
            nd = c_mat * (zc * sh) + d_mat * ch
        else:
            # Lumped wall shunt Y_w = P*l / Z_w split as two half-shunts around segment matrix (Sondhi 1974)
            pj = 2.0 * aj / rj
            yh = (0.5 * pj * lj) * inv_z_w
            t12 = zc * sh
            t11 = ch + yh * t12
            t21 = (sh / zc) + yh * (ch + t11)
            na = a_mat * t11 + b_mat * t21
            nb = a_mat * t12 + b_mat * t11
            nc = c_mat * t11 + d_mat * t21
            nd = c_mat * t12 + d_mat * t11
        a_mat, b_mat, c_mat, d_mat = na, nb, nc, nd

    a_out_m = float(np.sqrt(area_m2[-1] / np.pi))
    zc_out = float(air.rho * air.c / area_m2[-1])

    if radiation == "closed":
        warnings: list[str] = []
        z_in = a_mat / np.where(np.abs(c_mat) > 1e-30, c_mat, 1e-30 + 0j)
        h_transfer = 1.0 / np.where(np.abs(c_mat) > 1e-30, c_mat, 1e-30 + 0j)
    else:
        z_r, warnings = radiation_impedance(omega, a_out_m, zc_out, air, model=radiation)
        denom = c_mat * z_r + d_mat
        z_in = (a_mat * z_r + b_mat) / denom
        h_transfer = 1.0 / denom

    if not (np.all(np.isfinite(z_in)) and np.all(np.isfinite(h_transfer))):
        raise FloatingPointError("Non-finite values (Inf/NaN) encountered in TMM solution.")

    h_mag = np.abs(h_transfer)
    zin_mag = np.abs(z_in)
    h_peaks_all = _extract_peaks(freqs, h_mag)
    wall_compliance_peaks: list[dict] = []
    if use_wall and wall is not None and wall.stiffness_per_area > 0.0 and radiation != "closed":
        h_peaks = []
        for pk in h_peaks_all:
            if pk["f_hz"] < 120.0 and pk["mag_db"] < 15.0:
                wall_compliance_peaks.append(pk)
            else:
                h_peaks.append(pk)
        if not h_peaks:
            h_peaks = h_peaks_all
    else:
        h_peaks = h_peaks_all
    zin_peaks = _extract_peaks(freqs, zin_mag)

    l_total_m = float(s_m[-1] - s_m[0])
    if radiation == "unflanged":
        end_corr = 0.6133 * a_out_m
    elif radiation == "flanged":
        end_corr = 0.8216 * a_out_m
    else:
        end_corr = 0.0
    l_eff_m = l_total_m + end_corr

    wall_res_hz = None
    if use_wall and wall is not None and wall.mass_per_area > 0.0:
        a_bar = float(np.sum(seg_len * seg_area) / np.sum(seg_len))
        p_bar = float(np.sum(seg_len * (2.0 * seg_area / seg_rad)) / np.sum(seg_len))
        wall_res_hz = float((1.0 / (2.0 * np.pi)) * np.sqrt(air.rho * (air.c ** 2) * p_bar / (a_bar * wall.mass_per_area)))

    return {
        "freqs_hz": freqs,
        "z_in": z_in,
        "h_transfer": h_transfer,
        "h_peaks": h_peaks,
        "wall_compliance_peaks": wall_compliance_peaks,
        "zin_peaks": zin_peaks,
        "air": asdict(air),
        "radiation_model": radiation,
        "losses_enabled": bool(losses),
        "extra_damping_factor": float(extra_damping_factor),
        "wall": asdict(wall) if use_wall and wall is not None else None,
        "wall_resonance_closed_hz": None if wall_res_hz is None else round(wall_res_hz, 3),
        "n_segments": int(len(seg_len)),
        "max_segment_mm": round(float(np.max(seg_len) * 1e3), 4),
        "df_hz": float(df),
        "L_m": l_total_m,
        "L_eff_m": l_eff_m,
        "a_out_m": a_out_m,
        "closed_open_hz": [((2 * n - 1) * air.c) / (4.0 * l_total_m) for n in range(1, 5)],
        "closed_open_eff_hz": [((2 * n - 1) * air.c) / (4.0 * l_eff_m) for n in range(1, 5)],
        "open_open_hz": [(n * air.c) / (2.0 * l_total_m) for n in range(1, 5)],
        "warnings": warnings,
    }


def check_grid_independence(
    s_mm: np.ndarray,
    area_mm2: np.ndarray,
    dh_mm: np.ndarray,
    *,
    air: Air | None = None,
    fmin: float = 10.0,
    fmax: float = 4000.0,
    df: float = 0.5,
    radiation: str = "unflanged",
    losses: bool = True,
    wall: WallModel | None = None,
    extra_damping_factor: float = 1.0,
    max_freq_shift_pct: float = 0.2,
    max_mag_shift_db: float = 0.5,
) -> dict:
    """Solve at df and df/2 and verify peak frequency (<0.2%) and magnitude (<0.5 dB) convergence."""
    coarse = solve_tmm(
        s_mm,
        area_mm2,
        dh_mm,
        air=air,
        fmin=fmin,
        fmax=fmax,
        df=df,
        radiation=radiation,
        losses=losses,
        wall=wall,
        extra_damping_factor=extra_damping_factor,
    )
    fine = solve_tmm(
        s_mm,
        area_mm2,
        dh_mm,
        air=air,
        fmin=fmin,
        fmax=fmax,
        df=df / 2.0,
        radiation=radiation,
        losses=losses,
        wall=wall,
        extra_damping_factor=extra_damping_factor,
    )

    p1 = coarse["h_peaks"]
    p2 = fine["h_peaks"]
    max_f_pct = 0.0
    max_m_db = 0.0
    details: list[dict] = []

    for pk_c in p1:
        if not p2:
            raise AssertionError("Fine grid found no H(f) peaks during grid-independence check.")
        nearest = min(p2, key=lambda q: abs(q["f_hz"] - pk_c["f_hz"]))
        f_pct = 100.0 * abs(nearest["f_hz"] - pk_c["f_hz"]) / max(pk_c["f_hz"], 1e-9)
        m_db = abs(nearest["mag_db"] - pk_c["mag_db"])
        max_f_pct = max(max_f_pct, f_pct)
        max_m_db = max(max_m_db, m_db)
        details.append({"coarse_hz": pk_c["f_hz"], "fine_hz": nearest["f_hz"], "df_pct": round(f_pct, 5), "dmag_db": round(m_db, 5)})

    if max_f_pct >= max_freq_shift_pct or max_m_db >= max_mag_shift_db:
        raise AssertionError(
            f"TMM grid-independence failed (df={df} vs {df/2}): "
            f"max_f_shift={max_f_pct:.4f}% (limit {max_freq_shift_pct}%), "
            f"max_mag_shift={max_m_db:.4f} dB (limit {max_mag_shift_db} dB)."
        )

    fine["grid_independence"] = {
        "passed": True,
        "df_coarse_hz": float(df),
        "df_fine_hz": float(df / 2.0),
        "max_freq_shift_pct": round(max_f_pct, 6),
        "max_mag_shift_db": round(max_m_db, 6),
        "comparisons": details,
    }
    return fine


def plot_tmm_lossy(res: dict, out_png: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    freqs = res["freqs_hz"]
    zin_db = 20.0 * np.log10(np.maximum(np.abs(res["z_in"]), 1e-12))
    h_db = 20.0 * np.log10(np.maximum(np.abs(res["h_transfer"]), 1e-12))

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
    ax1.plot(freqs, h_db, color="navy", lw=1.5, label="|H(f)| = |U_out / U_in| (dB)")
    for idx, pk in enumerate(res["h_peaks"][:4]):
        ax1.plot(pk["f_hz"], pk["mag_db"], "ro", ms=5)
        ax1.annotate(
            f"F{idx+1}={pk['f_hz']:.1f} Hz (Q={pk['q']:.1f})",
            xy=(pk["f_hz"], pk["mag_db"]),
            xytext=(5, 4),
            textcoords="offset points",
            fontsize=8,
            color="crimson",
        )
    for idx, fq in enumerate(res["closed_open_hz"][:4]):
        if fq <= freqs[-1]:
            ax1.axvline(fq, color="teal", ls="--", lw=1.0, alpha=0.65, label="(2n-1)c/4L" if idx == 0 else None)
    for idx, fh in enumerate(res["open_open_hz"][:4]):
        if fh <= freqs[-1]:
            ax1.axvline(fh, color="gray", ls=":", lw=1.0, alpha=0.65, label="nc/2L" if idx == 0 else None)
    ax1.set_ylabel("|H(f)| (dB)")
    ax1.set_title(
        f"Lossy TMM ({res['radiation_model']}, T={res['air']['T_c']}°C, "
        f"{res['n_segments']} segs, Δf={res['df_hz']} Hz)"
    )
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc="upper right", fontsize=8)

    ax2.plot(freqs, zin_db, color="indigo", lw=1.5, label="|Z_in(f)| (dB ref 1 Pa·s/m³)")
    ax2.set_xlabel("Frequency (Hz)")
    ax2.set_ylabel("|Z_in(f)| (dB)")
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc="upper right", fontsize=8)

    fig.tight_layout()
    fig.savefig(out_png, dpi=96)
    plt.close(fig)


def run_tmm_profile(
    profile_csv: Path,
    out_dir: Path,
    *,
    fmax: float = 4000.0,
    df: float = 0.5,
    temperature_c: float = 20.0,
    radiation: str = "unflanged",
    write_outputs: bool = True,
    prof: dict | None = None,
    wall: WallModel | None = None,
    extra_damping_factor: float = 1.0,
) -> dict:
    prof = prof if prof is not None else load_area_profile(profile_csv)
    air = Air.from_temperature(temperature_c)
    res = check_grid_independence(
        prof["s_mm"],
        prof["area_mm2"],
        prof["dh_mm"],
        air=air,
        fmin=10.0,
        fmax=fmax,
        df=df,
        radiation=radiation,
        losses=True,
        wall=wall,
        extra_damping_factor=extra_damping_factor,
    )

    summary = {
        "profile_csv": str(profile_csv),
        "air": res["air"],
        "radiation_model": res["radiation_model"],
        "extra_damping_factor": res["extra_damping_factor"],
        "wall": res["wall"],
        "wall_resonance_closed_hz": res["wall_resonance_closed_hz"],
        "n_segments": res["n_segments"],
        "max_segment_mm": res["max_segment_mm"],
        "df_hz": res["df_hz"],
        "L_m": round(res["L_m"], 6),
        "L_eff_m": round(res["L_eff_m"], 6),
        "a_out_m": round(res["a_out_m"], 6),
        "peaks": res["h_peaks"],
        "zin_peaks": res["zin_peaks"],
        "analytic_closed_open_hz": [round(v, 3) for v in res["closed_open_hz"]],
        "analytic_closed_open_eff_hz": [round(v, 3) for v in res["closed_open_eff_hz"]],
        "analytic_open_open_hz": [round(v, 3) for v in res["open_open_hz"]],
        "grid_independence": res["grid_independence"],
        "warnings": res["warnings"],
    }

    if write_outputs:
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "tmm_lossy.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
        plot_tmm_lossy(res, out_dir / "tmm_lossy.png")

    return {"summary": summary, "raw": res}


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="M2: Lossy Kirchhoff TMM solver")
    parser.add_argument("--profile", required=True, help="CSV with s_mm (or distance_mm), area_mm2[, perimeter_mm]")
    parser.add_argument("--fmax", type=float, default=4000.0)
    parser.add_argument("--df", type=float, default=0.5)
    parser.add_argument("--temperature", type=float, default=20.0)
    parser.add_argument("--radiation", choices=["unflanged", "flanged", "ideal_open"], default="unflanged")
    parser.add_argument("--out", required=True, help="Output directory for tmm_lossy.json and tmm_lossy.png")
    args = parser.parse_args(argv)

    out = run_tmm_profile(
        Path(args.profile),
        Path(args.out),
        fmax=args.fmax,
        df=args.df,
        temperature_c=args.temperature,
        radiation=args.radiation,
        write_outputs=True,
    )
    print(json.dumps(out["summary"], indent=2))


if __name__ == "__main__":
    main()

"""M2 (Step 3) — Deterministic LHS + OAT soft-tissue sensitivity and uncertainty analysis.

Sweeps all range parameters in soft_tissue_hypotheses.json via seeded Latin
Hypercube Sampling (scipy.stats.qmc.LatinHypercube) crossed with temperature in
{20, 37} C, using the yielding-wall Kirchhoff TMM solver (df=0.5 Hz, fmax=3000 Hz).
Computes OAT tornado sensitivities, Spearman/PRCC rankings, spaghetti/distribution
plots, and renders 3 scenario WAVs (min, median, max mean bandwidth).
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import time
from pathlib import Path

import numpy as np
from scipy import signal, stats
from scipy.stats import qmc

from pipeline.epb import apply_epb_table_to_config
from pipeline.geometry_qc import apply_end_section_trim, evaluate_geometry_qc
from pipeline.hypothesis import load_metadata
from pipeline.measure import (
    band_energy_fractions,
    measure_f0,
    measure_transfer_function,
    read_wav_float,
)
from pipeline.paths import MODELS, OUTPUT, ROOT
from pipeline.soft_tissue import (
    SoftTissueConfig,
    SoftTissueError,
    apply_lumen_narrowing,
    load_hypotheses,
)
from pipeline.source_tract import (
    _interactive_two_mass_flow,
    build_kl_waveguide,
    kl_impulse_response,
    normalize_audio,
    synthesize_kl,
)
from pipeline.tmm_lossy import Air, WallModel, solve_tmm

logger = logging.getLogger(__name__)

SWEPT_CONTINUOUS_PARAMS = (
    "mucosa_thickness_mm",
    "wall_mass_per_area",
    "wall_stiffness_per_area",
    "wall_resistance_per_area",
    "extra_damping_factor",
    "f0_prior_hz",
)
ALL_SWEPT_PARAMS = SWEPT_CONTINUOUS_PARAMS + ("air_temperature_c",)

SAMPLE_CSV_COLUMNS = [
    "sample_id",
    "mucosa_thickness_mm",
    "wall_mass_per_area",
    "wall_stiffness_per_area",
    "wall_resistance_per_area",
    "extra_damping_factor",
    "f0_prior_hz",
    "air_temperature_c",
    "valid",
    "invalid_reason",
    "F1_hz",
    "F2_hz",
    "F3_hz",
    "Q1",
    "Q2",
    "Q3",
    "B1_hz",
    "B2_hz",
    "B3_hz",
    "mean_bw_hz",
    "wall_resonance_closed_hz",
    "shaping_gain_db_F1_b3db_predicted",
]


def _step2_source_psd(prof_bone: dict) -> tuple[np.ndarray, np.ndarray, float]:
    """Load or deterministically synthesize the step-2 baseline glottal source PSD."""
    step2_src = ROOT / "output" / "parasaurolophus_ram14000_20261117T104254Z_19c2be" / "source.wav"
    if step2_src.exists():
        sr, x, _ = read_wav_float(step2_src)
        nseg = min(16384, len(x))
        freqs, psd = signal.welch(x, fs=sr, window="hann", nperseg=nseg, noverlap=nseg // 2)
        return freqs, psd, 60.7601

    air = Air.from_temperature(20.0)
    wg = build_kl_waveguide(
        prof_bone["s_mm"],
        prof_bone["area_mm2"],
        air,
        dh_mm=prof_bone["dh_mm"],
        radiation="unflanged",
        sr=44100,
        kl_loss_model="legacy",
    )
    src_wg, _ = _interactive_two_mass_flow(wg, air, duration_s=1.5, f0_prior_hz=60.0, seed=0)
    src_44 = signal.resample(src_wg, int(round(1.5 * 44100)))
    src_norm, _, _ = normalize_audio(src_44, target_dbfs=-1.0)
    nseg = min(16384, len(src_norm))
    freqs, psd = signal.welch(src_norm, fs=44100, window="hann", nperseg=nseg, noverlap=nseg // 2)
    return freqs, psd, 60.0


def _predicted_shaping_gain_db(
    tmm_res: dict,
    src_freqs: np.ndarray,
    src_psd: np.ndarray,
    f0_ref_hz: float,
    f0_sample_hz: float,
    f1_hz: float,
    b1_hz: float,
) -> float:
    """Compute predicted shaping_gain_db_F1_b3db by applying |H_TMM(f)|^2 to the baseline source PSD."""
    freqs = tmm_res["freqs_hz"]
    h_sq = np.abs(tmm_res["h_transfer"]) ** 2
    # Shift harmonic axis proportionally to f0_sample/f0_ref while preserving overall spectral range
    scale = float(f0_ref_hz) / max(float(f0_sample_hz), 1e-6)
    p_src = np.interp(freqs * scale, src_freqs, src_psd, left=float(src_psd[0]), right=float(src_psd[-1]))
    p_flt = p_src * h_sq
    tot_src = float(np.sum(p_src))
    tot_flt = float(np.sum(p_flt))
    if tot_src <= 0.0 or tot_flt <= 0.0:
        return 0.0
    min_hw = 2.0 * 44100.0 / 16384.0
    half = max(min_hw, 0.5 * float(b1_hz))
    mask = (freqs >= max(0.0, f1_hz - half)) & (freqs <= f1_hz + half)
    e_src = float(np.sum(p_src[mask]) / tot_src)
    e_flt = float(np.sum(p_flt[mask]) / tot_flt)
    return float(10.0 * np.log10(max(e_flt, 1e-12) / max(e_src, 1e-12)))


def _evaluate_sample_tmm(
    sample: dict,
    prof_bone: dict,
    min_area_fraction: float,
    src_freqs: np.ndarray,
    src_psd: np.ndarray,
    f0_ref_hz: float,
    *,
    fmax: float = 3000.0,
    df: float = 0.5,
) -> tuple[dict, np.ndarray | None, np.ndarray | None]:
    """Run lumen narrowing + yielding-wall TMM for one parameter vector; never clips occluded lumens."""
    row = dict(sample)
    try:
        prof_lumen = apply_lumen_narrowing(
            prof_bone,
            float(sample["mucosa_thickness_mm"]),
            min_area_fraction=min_area_fraction,
        )
        air = Air.from_temperature(float(sample["air_temperature_c"]))
        wall = WallModel(
            mass_per_area=float(sample["wall_mass_per_area"]),
            stiffness_per_area=float(sample["wall_stiffness_per_area"]),
            resistance_per_area=float(sample["wall_resistance_per_area"]),
        )
        tmm = solve_tmm(
            prof_lumen["s_mm"],
            prof_lumen["area_mm2"],
            prof_lumen["dh_mm"],
            air=air,
            fmin=10.0,
            fmax=fmax,
            df=df,
            radiation="unflanged",
            losses=True,
            wall=wall,
            extra_damping_factor=float(sample["extra_damping_factor"]),
        )
    except (SoftTissueError, FloatingPointError, ValueError) as exc:
        row.update(
            {
                "valid": False,
                "invalid_reason": str(exc),
                "F1_hz": None,
                "F2_hz": None,
                "F3_hz": None,
                "Q1": None,
                "Q2": None,
                "Q3": None,
                "B1_hz": None,
                "B2_hz": None,
                "B3_hz": None,
                "mean_bw_hz": None,
                "wall_resonance_closed_hz": None,
                "shaping_gain_db_F1_b3db_predicted": None,
            }
        )
        return row, None, None

    pks = tmm["h_peaks"]
    if len(pks) < 3:
        row.update(
            {
                "valid": False,
                "invalid_reason": f"fewer than 3 TMM resonance peaks detected ({len(pks)}) below {fmax:.0f} Hz",
                "F1_hz": pks[0]["f_hz"] if len(pks) > 0 else None,
                "F2_hz": pks[1]["f_hz"] if len(pks) > 1 else None,
                "F3_hz": None,
                "Q1": pks[0]["q"] if len(pks) > 0 else None,
                "Q2": pks[1]["q"] if len(pks) > 1 else None,
                "Q3": None,
                "B1_hz": pks[0]["bw_3db_hz"] if len(pks) > 0 else None,
                "B2_hz": pks[1]["bw_3db_hz"] if len(pks) > 1 else None,
                "B3_hz": None,
                "mean_bw_hz": None,
                "wall_resonance_closed_hz": tmm["wall_resonance_closed_hz"],
                "shaping_gain_db_F1_b3db_predicted": None,
            }
        )
        return row, None, None

    gain_pred = _predicted_shaping_gain_db(
        tmm,
        src_freqs,
        src_psd,
        f0_ref_hz,
        float(sample["f0_prior_hz"]),
        float(pks[0]["f_hz"]),
        float(pks[0]["bw_3db_hz"]),
    )
    mean_bw = float((pks[0]["bw_3db_hz"] + pks[1]["bw_3db_hz"] + pks[2]["bw_3db_hz"]) / 3.0)
    row.update(
        {
            "valid": True,
            "invalid_reason": "",
            "F1_hz": round(float(pks[0]["f_hz"]), 3),
            "F2_hz": round(float(pks[1]["f_hz"]), 3),
            "F3_hz": round(float(pks[2]["f_hz"]), 3),
            "Q1": round(float(pks[0]["q"]), 3),
            "Q2": round(float(pks[1]["q"]), 3),
            "Q3": round(float(pks[2]["q"]), 3),
            "B1_hz": round(float(pks[0]["bw_3db_hz"]), 3),
            "B2_hz": round(float(pks[1]["bw_3db_hz"]), 3),
            "B3_hz": round(float(pks[2]["bw_3db_hz"]), 3),
            "mean_bw_hz": round(mean_bw, 4),
            "wall_resonance_closed_hz": tmm["wall_resonance_closed_hz"],
            "shaping_gain_db_F1_b3db_predicted": round(gain_pred, 3),
        }
    )
    h_db = 20.0 * np.log10(np.maximum(np.abs(tmm["h_transfer"]), 1e-12))
    return row, tmm["freqs_hz"], h_db


def _generate_lhs_samples(cfg: SoftTissueConfig, n_samples: int, seed: int) -> list[dict]:
    """Generate deterministic Latin Hypercube samples crossed with temperature in {20, 37} C."""
    params = cfg.parameters()
    sampler = qmc.LatinHypercube(d=len(SWEPT_CONTINUOUS_PARAMS), seed=int(seed))
    u = sampler.random(n=int(n_samples))

    rng = np.random.default_rng(int(seed))
    temp_hyp = params["air_temperature_c"]
    disc_temps = list(temp_hyp.discrete_values) if temp_hyp.discrete_values else [20.0, 37.0]
    temps = np.array([disc_temps[i % len(disc_temps)] for i in range(n_samples)], dtype=np.float64)
    rng.shuffle(temps)

    samples: list[dict] = []
    for i in range(n_samples):
        item: dict = {"sample_id": i}
        for col_idx, name in enumerate(SWEPT_CONTINUOUS_PARAMS):
            h = params[name]
            lo, hi = float(h.range[0]), float(h.range[1])
            ui = float(u[i, col_idx])
            if h.scale == "log":
                val = float(np.exp(np.log(lo) + ui * (np.log(hi) - np.log(lo))))
            else:
                val = lo + ui * (hi - lo)
            item[name] = round(val, 6)
        item["air_temperature_c"] = float(temps[i])
        samples.append(item)
    return samples


def _max_valid_thickness(prof_bone: dict, lo_mm: float, hi_mm: float, min_area_fraction: float) -> float:
    """Deterministically find highest non-occluded mucosa thickness via bisection."""
    a, b = float(lo_mm), float(hi_mm)
    for _ in range(16):
        mid = 0.5 * (a + b)
        try:
            apply_lumen_narrowing(prof_bone, mid, min_area_fraction=min_area_fraction)
            a = mid
        except SoftTissueError:
            b = mid
    return round(a, 4)


def _compute_importance(
    valid_rows: list[dict],
    oat_records: list[dict],
) -> dict:
    """Compute Spearman rank correlations, PRCCs, and combined ranking for F1, Q1, and predicted gain."""
    targets = ("F1_hz", "Q1", "shaping_gain_db_F1_b3db_predicted")
    spearman_out: dict[str, dict] = {}
    prcc_out: dict[str, dict] = {}

    if len(valid_rows) >= 5:
        x_mat = np.array(
            [[float(r[p]) for p in ALL_SWEPT_PARAMS] for r in valid_rows],
            dtype=np.float64,
        )
        # Rank-transform predictors once for PRCC
        r_x = np.apply_along_axis(stats.rankdata, 0, x_mat)
        for tgt in targets:
            y_vec = np.array([float(r[tgt]) for r in valid_rows], dtype=np.float64)
            r_y = stats.rankdata(y_vec)
            sp_map: dict[str, dict] = {}
            pr_map: dict[str, float] = {}

            for j, pname in enumerate(ALL_SWEPT_PARAMS):
                if np.std(x_mat[:, j]) < 1e-12 or np.std(y_vec) < 1e-12:
                    sp_map[pname] = {"rho": 0.0, "p_value": 1.0}
                    pr_map[pname] = 0.0
                    continue
                rho, pval = stats.spearmanr(x_mat[:, j], y_vec)
                sp_map[pname] = {
                    "rho": round(float(0.0 if np.isnan(rho) else rho), 4),
                    "p_value": round(float(1.0 if np.isnan(pval) else pval), 6),
                }
                # Partial rank correlation controlling for remaining predictors
                other_cols = [c for c in range(r_x.shape[1]) if c != j]
                z_ctrl = np.column_stack([np.ones(len(r_y)), r_x[:, other_cols]])
                bx, *_ = np.linalg.lstsq(z_ctrl, r_x[:, j], rcond=None)
                by, *_ = np.linalg.lstsq(z_ctrl, r_y, rcond=None)
                res_x = r_x[:, j] - z_ctrl @ bx
                res_y = r_y - z_ctrl @ by
                denom = float(np.sqrt(np.sum(res_x ** 2) * np.sum(res_y ** 2)))
                prcc_val = float(np.dot(res_x, res_y) / denom) if denom > 1e-12 else 0.0
                pr_map[pname] = round(prcc_val, 4)

            spearman_out[tgt] = sp_map
            prcc_out[tgt] = pr_map

    rankings: dict[str, list[dict]] = {}
    for tgt, delta_key in (
        ("F1_hz", "delta_F1_hz"),
        ("Q1", "delta_Q1"),
        ("shaping_gain_db_F1_b3db_predicted", "delta_gain_db"),
    ):
        oat_by_param = {rec["parameter"]: rec for rec in oat_records}
        sp_map = spearman_out.get(tgt, {})
        pr_map = prcc_out.get(tgt, {})
        ranked = []
        for pname in ALL_SWEPT_PARAMS:
            rho_val = float((sp_map.get(pname) or {}).get("rho", 0.0))
            prcc_val = float(pr_map.get(pname, 0.0))
            oat_span = float((oat_by_param.get(pname) or {}).get(delta_key, 0.0))
            ranked.append(
                {
                    "parameter": pname,
                    "abs_spearman_rho": round(abs(rho_val), 4),
                    "spearman_rho": round(rho_val, 4),
                    "prcc": round(prcc_val, 4),
                    "oat_span": round(oat_span, 3),
                }
            )
        ranked.sort(key=lambda item: (item["abs_spearman_rho"], item["oat_span"]), reverse=True)
        rankings[tgt] = ranked

    return {
        "spearman": spearman_out,
        "prcc": prcc_out,
        "rankings": rankings,
    }


def _summary_intervals(valid_rows: list[dict]) -> dict[str, dict]:
    keys = [
        "F1_hz",
        "F2_hz",
        "F3_hz",
        "Q1",
        "Q2",
        "Q3",
        "B1_hz",
        "B2_hz",
        "B3_hz",
        "mean_bw_hz",
        "wall_resonance_closed_hz",
        "shaping_gain_db_F1_b3db_predicted",
    ]
    intervals: dict[str, dict] = {}
    for key in keys:
        vals = [float(r[key]) for r in valid_rows if r.get(key) is not None]
        if not vals:
            intervals[key] = {"n": 0, "median": None, "p05": None, "p95": None, "min": None, "max": None}
        else:
            arr = np.asarray(vals, dtype=np.float64)
            intervals[key] = {
                "n": int(len(arr)),
                "median": round(float(np.median(arr)), 4),
                "p05": round(float(np.percentile(arr, 5)), 4),
                "p95": round(float(np.percentile(arr, 95)), 4),
                "min": round(float(np.min(arr)), 4),
                "max": round(float(np.max(arr)), 4),
            }
    return intervals


def _plot_tornado(oat_records: list[dict], out_png: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    labels_fr = {
        "mucosa_thickness_mm": "Épaisseur muqueuse (mm)",
        "wall_mass_per_area": "Masse pariétale M_w",
        "wall_stiffness_per_area": "Raideur pariétale K_w",
        "wall_resistance_per_area": "Résistance pariétale R_w",
        "extra_damping_factor": "Facteur d'amortissement",
        "f0_prior_hz": "Fréquence f0 source (Hz)",
        "air_temperature_c": "Température air (°C)",
    }

    fig, axes = plt.subplots(1, 3, figsize=(15, 5.2), sharey=True)
    panels = [
        ("F1_hz", "delta_F1_hz", "Variation de F1 (Hz)"),
        ("Q1", "delta_Q1", "Variation du facteur de qualité Q1"),
        ("shaping_gain_db_F1_b3db_predicted", "delta_gain_db", "Variation du gain de filtrage F1 (dB)"),
    ]
    # Order parameters by F1 span
    ordered = sorted(oat_records, key=lambda r: r["delta_F1_hz"])
    y_pos = np.arange(len(ordered))
    y_labels = [labels_fr.get(r["parameter"], r["parameter"]) for r in ordered]

    for ax, (metric_k, _, title) in zip(axes, panels):
        lo_deltas = [r[f"{metric_k}_lo_delta"] for r in ordered]
        hi_deltas = [r[f"{metric_k}_hi_delta"] for r in ordered]
        ax.barh(y_pos, lo_deltas, height=0.45, color="steelblue", alpha=0.85, label="Borne basse")
        ax.barh(y_pos, hi_deltas, height=0.45, color="darkorange", alpha=0.85, label="Borne haute")
        ax.axvline(0.0, color="black", lw=1.0)
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("Écart vs valeur nominale")
        ax.grid(True, axis="x", alpha=0.3)
    axes[0].set_yticks(y_pos)
    axes[0].set_yticklabels(y_labels, fontsize=9)
    axes[2].legend(loc="lower right", fontsize=8)
    fig.suptitle("Analyse de sensibilité univariée (Tornado) autour des hypothèses nominales", fontsize=12)
    fig.tight_layout()
    fig.savefig(out_png, dpi=96)
    plt.close(fig)


def _plot_spaghetti(
    freqs: np.ndarray,
    valid_h_db: list[np.ndarray],
    nom_h_db: np.ndarray | None,
    bone_h_db: np.ndarray,
    out_png: Path,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(10.5, 5.2))
    for h_db in valid_h_db:
        ax.plot(freqs, h_db, color="0.78", lw=0.7, alpha=0.55)
    if valid_h_db:
        mat = np.vstack(valid_h_db)
        p05 = np.percentile(mat, 5, axis=0)
        p95 = np.percentile(mat, 95, axis=0)
        ax.fill_between(freqs, p05, p95, color="lightsteelblue", alpha=0.35, label="Enveloppe 5–95 % (hypothèses)")
    ax.plot(freqs, bone_h_db, color="navy", lw=1.4, ls="--", label="Référence os seul (étape 2, paroi rigide)")
    if nom_h_db is not None:
        ax.plot(freqs, nom_h_db, color="black", lw=1.8, label="Scénario nominal tissus mous (exploratoire)")
    ax.set_xlim(10.0, float(freqs[-1]))
    ax.set_xlabel("Fréquence (Hz)")
    ax.set_ylabel("|H_TMM(f)| (dB)")
    ax.set_title("Enveloppe spectrale de sensibilité : os seul vs hypothèses de tissus mous (RAM 14000)")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_png, dpi=96)
    plt.close(fig)


def _plot_distributions(
    valid_rows: list[dict],
    nom_row: dict,
    bone_row: dict,
    out_png: Path,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 3, figsize=(13.5, 7.0))
    specs = [
        ("F1_hz", "F1 (Hz)"),
        ("F2_hz", "F2 (Hz)"),
        ("F3_hz", "F3 (Hz)"),
        ("Q1", "Q1"),
        ("Q2", "Q2"),
        ("Q3", "Q3"),
    ]
    for ax, (key, label) in zip(axes.flat, specs):
        vals = [float(r[key]) for r in valid_rows if r.get(key) is not None]
        if vals:
            ax.hist(vals, bins=min(15, max(5, len(vals) // 2)), color="steelblue", alpha=0.75, edgecolor="white")
        if bone_row.get(key) is not None:
            ax.axvline(float(bone_row[key]), color="navy", ls="--", lw=1.6, label=f"Os seul: {bone_row[key]:.1f}")
        if nom_row.get(key) is not None:
            ax.axvline(float(nom_row[key]), color="crimson", ls="-", lw=1.6, label=f"Nominal: {nom_row[key]:.1f}")
        ax.set_xlabel(label)
        ax.set_ylabel("Effectif")
        ax.grid(True, alpha=0.25)
        ax.legend(fontsize=7.5)
    fig.suptitle("Distributions des fréquences formantiques et facteurs Q (échantillonnage LHS)", fontsize=12)
    fig.tight_layout()
    fig.savefig(out_png, dpi=96)
    plt.close(fig)


def _render_scenario_audio(
    tag: str,
    sample_row: dict,
    prof_bone: dict,
    min_area_fraction: float,
    out_dir: Path,
    seed: int,
    duration_s: float,
) -> dict:
    """Synthesize and measure one sensitivity damping scenario WAV + sidecar JSON."""
    prof_lumen = apply_lumen_narrowing(
        prof_bone,
        float(sample_row["mucosa_thickness_mm"]),
        min_area_fraction=min_area_fraction,
    )
    air = Air.from_temperature(float(sample_row["air_temperature_c"]))
    wall = WallModel(
        mass_per_area=float(sample_row["wall_mass_per_area"]),
        stiffness_per_area=float(sample_row["wall_stiffness_per_area"]),
        resistance_per_area=float(sample_row["wall_resistance_per_area"]),
    )
    tmm = solve_tmm(
        prof_lumen["s_mm"],
        prof_lumen["area_mm2"],
        prof_lumen["dh_mm"],
        air=air,
        fmax=3000.0,
        df=0.5,
        radiation="unflanged",
        losses=True,
        wall=wall,
        extra_damping_factor=float(sample_row["extra_damping_factor"]),
    )
    synth = synthesize_kl(
        prof_lumen["s_mm"],
        prof_lumen["area_mm2"],
        air,
        {"peaks": tmm["h_peaks"]},
        dh_mm=prof_lumen["dh_mm"],
        perimeter_mm=prof_lumen.get("perimeter_mm"),
        duration_s=duration_s,
        sr=44100,
        f0_hz=float(sample_row["f0_prior_hz"]),
        radiation="unflanged",
        seed=seed,
        kl_loss_model="kirchhoff_fit",
        extra_damping_factor=float(sample_row["extra_damping_factor"]),
        wall=wall,
    )

    from scipy.io import wavfile

    wav_name = f"scenario_damping_{tag}.wav"
    wav_path = out_dir / wav_name
    arr = np.asarray(synth["filtered"], dtype=np.float64)
    assert float(np.max(np.abs(arr))) < 0.999
    pcm = np.int16(np.round(np.clip(arr, -1.0, 1.0) * 32767.0))
    wavfile.write(wav_path, 44100, pcm)

    tf = measure_transfer_function(synth["source"], synth["filtered"], 44100)
    f0_meas = measure_f0(synth["source"], 44100)
    bands_flt = band_energy_fractions(synth["filtered"], 44100)

    sidecar = {
        "kind": "hypothesis_scenario",
        "scenario": tag,
        "wav_file": wav_name,
        "sample_id": sample_row["sample_id"],
        "parameters": {k: sample_row[k] for k in ALL_SWEPT_PARAMS},
        "predicted_tmm_peaks": tmm["h_peaks"][:3],
        "measured_f0_hz": f0_meas["measured_hz"],
        "measured_hf_peaks": tf["peaks"][:4],
        "band_energy_fractions_filtered": bands_flt,
        "audio_peak_dbfs": synth["audio"]["filtered_peak_dbfs"],
        "clipping_detected": synth["audio"]["clipping_detected"],
    }
    sidecar_path = out_dir / f"scenario_damping_{tag}.json"
    sidecar_path.write_text(json.dumps(sidecar, indent=2) + "\n", encoding="utf-8")
    return sidecar


def run_sensitivity(
    taxon_slug: str,
    *,
    n_samples: int = 64,
    seed: int = 0,
    out_dir: Path | None = None,
    soft_tissue_path: Path | None = None,
    epb_table: Path | None = None,
    trim_end_sections: str = "0",
    duration_s: float = 2.5,
    render_audio: bool = True,
    render_plots: bool = True,
) -> dict:
    """Execute deterministic soft-tissue sensitivity sweep, plots, and 3 damping scenario WAVs."""
    t_start = time.perf_counter()
    meta = load_metadata(taxon_slug)
    hyp_path = Path(soft_tissue_path) if soft_tissue_path is not None else (MODELS / taxon_slug / "soft_tissue_hypotheses.json")
    cfg = load_hypotheses(hyp_path)
    cfg, epb_info, epb_warnings = apply_epb_table_to_config(cfg, epb_table)

    profile_csv = ROOT / meta["geometry"]["airway_features_csv"]
    qc, prof_bone = evaluate_geometry_qc(
        profile_csv,
        source_mesh=meta["geometry"]["stl_path"],
        extraction_method="bone_mesh_skeleton_axis_slice",
        require_sidecar=False,
    )
    qc, prof_bone, trim_msgs = apply_end_section_trim(qc, prof_bone, trim_end_sections)

    if out_dir is None:
        out_dir = OUTPUT / f"{taxon_slug}_sensitivity_n{n_samples}_seed{seed}"
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    src_freqs, src_psd, f0_ref_hz = _step2_source_psd(prof_bone)
    nom_vals = cfg.nominal_values()
    min_area_frac = float(nom_vals["lumen_min_area_fraction"])

    # Baseline bone-only step-2 rigid curve and nominal soft-tissue curve
    bone_sample = {
        "sample_id": -2,
        "mucosa_thickness_mm": 0.0,
        "wall_mass_per_area": 0.0,
        "wall_stiffness_per_area": 0.0,
        "wall_resistance_per_area": 0.0,
        "extra_damping_factor": 1.0,
        "f0_prior_hz": 60.0,
        "air_temperature_c": 20.0,
    }
    bone_row, freqs_grid, bone_h_db = _evaluate_sample_tmm(
        bone_sample, prof_bone, min_area_frac, src_freqs, src_psd, f0_ref_hz
    )
    nom_sample = {"sample_id": -1, **{k: float(nom_vals[k]) for k in ALL_SWEPT_PARAMS}}
    nom_row, _, nom_h_db = _evaluate_sample_tmm(
        nom_sample, prof_bone, min_area_frac, src_freqs, src_psd, f0_ref_hz
    )

    # 1. Seeded LHS sweep
    lhs_inputs = _generate_lhs_samples(cfg, n_samples, seed)
    sample_rows: list[dict] = []
    valid_rows: list[dict] = []
    valid_h_db: list[np.ndarray] = []
    invalid_reasons: dict[str, int] = {}

    for smp in lhs_inputs:
        row, _, h_db = _evaluate_sample_tmm(
            smp, prof_bone, min_area_frac, src_freqs, src_psd, f0_ref_hz
        )
        sample_rows.append(row)
        if row["valid"]:
            valid_rows.append(row)
            if h_db is not None:
                valid_h_db.append(h_db)
        else:
            rsn = row["invalid_reason"] or "unknown"
            key = rsn.split(":")[0]
            invalid_reasons[key] = invalid_reasons.get(key, 0) + 1

    # Write deterministic samples.csv
    csv_path = out_dir / "samples.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=SAMPLE_CSV_COLUMNS)
        w.writeheader()
        for r in sample_rows:
            w.writerow({col: r.get(col) for col in SAMPLE_CSV_COLUMNS})

    # 2. OAT sensitivity around nominal
    params = cfg.parameters()
    oat_records: list[dict] = []
    for pname in ALL_SWEPT_PARAMS:
        h = params[pname]
        if pname == "air_temperature_c" and h.discrete_values:
            lo_val, hi_val = float(min(h.discrete_values)), float(max(h.discrete_values))
        else:
            lo_val, hi_val = float(h.range[0]), float(h.range[1])

        s_lo = dict(nom_sample)
        s_lo[pname] = lo_val
        r_lo, _, _ = _evaluate_sample_tmm(s_lo, prof_bone, min_area_frac, src_freqs, src_psd, f0_ref_hz)

        s_hi = dict(nom_sample)
        s_hi[pname] = hi_val
        r_hi, _, _ = _evaluate_sample_tmm(s_hi, prof_bone, min_area_frac, src_freqs, src_psd, f0_ref_hz)
        hi_clamped = None
        if not r_hi["valid"] and pname == "mucosa_thickness_mm":
            hi_clamped = _max_valid_thickness(prof_bone, nom_sample[pname], hi_val, min_area_frac)
            s_hi[pname] = hi_clamped
            r_hi, _, _ = _evaluate_sample_tmm(s_hi, prof_bone, min_area_frac, src_freqs, src_psd, f0_ref_hz)

        rec = {
            "parameter": pname,
            "lo_value": round(lo_val, 6),
            "hi_value": round(float(s_hi[pname]), 6),
            "hi_clamped_before_occlusion": hi_clamped,
            "lo_valid": bool(r_lo["valid"]),
            "hi_valid": bool(r_hi["valid"]),
        }
        for mk, dk in (
            ("F1_hz", "delta_F1_hz"),
            ("Q1", "delta_Q1"),
            ("shaping_gain_db_F1_b3db_predicted", "delta_gain_db"),
        ):
            v_nom = float(nom_row[mk]) if nom_row.get(mk) is not None else 0.0
            d_lo = (float(r_lo[mk]) - v_nom) if r_lo.get(mk) is not None else 0.0
            d_hi = (float(r_hi[mk]) - v_nom) if r_hi.get(mk) is not None else 0.0
            rec[f"{mk}_lo_delta"] = round(d_lo, 4)
            rec[f"{mk}_hi_delta"] = round(d_hi, 4)
            rec[dk] = round(float(max(abs(d_lo), abs(d_hi), abs(d_hi - d_lo))), 4)
        oat_records.append(rec)

    importance = _compute_importance(valid_rows, oat_records)
    intervals = _summary_intervals(valid_rows)

    # 3. Plots
    tornado_png = out_dir / "tornado.png"
    spaghetti_png = out_dir / "spaghetti_H.png"
    dist_png = out_dir / "formant_distributions.png"
    if render_plots:
        _plot_tornado(oat_records, tornado_png)
        if freqs_grid is not None and bone_h_db is not None:
            _plot_spaghetti(freqs_grid, valid_h_db, nom_h_db, bone_h_db, spaghetti_png)
        _plot_distributions(valid_rows, nom_row, bone_row, dist_png)

    # 4. Audio for 3 damping scenarios (min, median, max of mean B1-B3)
    scenarios: list[dict] = []
    if render_audio and valid_rows:
        ordered_bw = sorted(valid_rows, key=lambda r: (float(r["mean_bw_hz"]), int(r["sample_id"])))
        picks = [
            ("min", ordered_bw[0]),
            ("median", ordered_bw[len(ordered_bw) // 2]),
            ("max", ordered_bw[-1]),
        ]
        for tag, smp_row in picks:
            scenarios.append(
                _render_scenario_audio(tag, smp_row, prof_bone, min_area_frac, out_dir, seed, duration_s)
            )

    summary = {
        "schema_version": "1.0",
        "taxon_slug": taxon_slug,
        "geometry_status": qc["status"],
        "state_effect": "none (sensitivity scenarios are exploratory and cannot upgrade state)",
        "soft_tissue_file": str(hyp_path),
        "epb": epb_info,
        "seed": int(seed),
        "n_samples": int(n_samples),
        "n_valid": int(len(valid_rows)),
        "n_invalid": int(len(sample_rows) - len(valid_rows)),
        "invalid_reasons": invalid_reasons,
        "bone_only_baseline": bone_row,
        "nominal_soft_tissue": nom_row,
        "intervals_p05_median_p95": intervals,
        "oat_sensitivity": oat_records,
        "importance": importance,
        "scenarios": scenarios,
        "warnings": epb_warnings + trim_msgs,
        "artifacts": {
            "samples_csv": "samples.csv",
            "tornado_png": "tornado.png",
            "spaghetti_png": "spaghetti_H.png",
            "formant_distributions_png": "formant_distributions.png",
        },
        "runtime_log": "sensitivity_log.txt",
    }
    (out_dir / "sensitivity.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    elapsed = time.perf_counter() - t_start
    (out_dir / "sensitivity_log.txt").write_text(
        f"taxon_slug: {taxon_slug}\n"
        f"seed: {seed}\n"
        f"n_samples: {n_samples}\n"
        f"n_valid: {len(valid_rows)}\n"
        f"n_invalid: {len(sample_rows) - len(valid_rows)}\n"
        f"runtime_wall_s: {elapsed:.3f}\n",
        encoding="utf-8",
    )
    return summary


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="M2: Soft-tissue sensitivity / uncertainty sweep")
    parser.add_argument("--taxon", required=True, help="Slug under models/")
    parser.add_argument("--n-samples", type=int, default=64)
    parser.add_argument("--quick", action="store_true", help="Run fast 16-sample sensitivity check")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--out", type=Path, default=None)
    parser.add_argument("--soft-tissue", type=Path, default=None)
    parser.add_argument("--epb-table", type=Path, default=None)
    parser.add_argument("--trim-end-sections", default="0")
    parser.add_argument("--duration", type=float, default=2.5)
    args = parser.parse_args(argv)

    n = 16 if args.quick else int(args.n_samples)
    summary = run_sensitivity(
        args.taxon,
        n_samples=n,
        seed=args.seed,
        out_dir=args.out,
        soft_tissue_path=args.soft_tissue,
        epb_table=args.epb_table,
        trim_end_sections=args.trim_end_sections,
        duration_s=args.duration,
        render_audio=True,
    )
    print(
        json.dumps(
            {
                "n_samples": summary["n_samples"],
                "n_valid": summary["n_valid"],
                "n_invalid": summary["n_invalid"],
                "rankings": summary["importance"]["rankings"],
                "intervals_p05_median_p95": summary["intervals_p05_median_p95"],
                "scenarios": [
                    {
                        "scenario": s["scenario"],
                        "wav_file": s["wav_file"],
                        "audio_peak_dbfs": s["audio_peak_dbfs"],
                        "clipping_detected": s["clipping_detected"],
                    }
                    for s in summary["scenarios"]
                ],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

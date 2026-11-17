"""État 3 — PEER_REVIEW: non-circular 3-way physics review (measured vs TMM vs analytic)."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from pipeline.measure import (
    compute_shaping_gain_db,
    compute_shaping_metrics,
    measure_run_dir,
    read_wav_float,
)
from pipeline.paths import OUTPUT, ROOT
from pipeline.tmm_lossy import run_tmm_profile


REVIEW_CONFIG = {
    "f1_measured_vs_tmm_max_rel_err": 0.05,
    "f2_f3_measured_vs_tmm_max_rel_err": 0.08,
    "kl_vs_tmm_max_q_rel_err": 0.20,
    "tmm_vs_analytic_info_rel_err": 0.15,
    "tmm_vs_analytic_fail_rel_err": 0.30,
    "min_shaping_gain_db": 3.0,
    "min_f1_band_energy_fraction": 0.02,
    "provenance_eps": 1e-6,
    "max_compare_freq_hz": 2000.0,
}

REQUIRED = ("source.wav", "filtered.wav", "metrics.json", "spectrogram.png")


def normalize_legacy_state(state: str) -> str:
    """Map legacy VALIDATED state to CONTRACT_OK for backward compatibility."""
    return "CONTRACT_OK" if state == "VALIDATED" else state


def resolve_terminal_state(passed: bool, geometry_status: str | None) -> str:
    """Resolve terminal pipeline state from review pass/fail and geometry QC status."""
    if not passed:
        return "REVIEW_FAILED"
    if geometry_status == "AIRWAY_OK":
        return "SCIENTIFICALLY_VALIDATED"
    if geometry_status == "BONE_DERIVED_EXPLORATORY":
        return "EXPLORATORY"
    return "CONTRACT_OK"


def _load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _wav_issues(path: Path) -> list[str]:
    if not path.exists():
        return [f"missing audio: {path.name}"]
    _, x, clipped = read_wav_float(path)
    issues: list[str] = []
    if x.size == 0:
        return [f"{path.name}: empty wav"]
    if not np.isfinite(x).all():
        issues.append(f"{path.name}: non-finite samples")
    peak = float(np.max(np.abs(x)))
    if clipped or peak >= 0.999:
        issues.append(f"{path.name}: clipping detected (peak={peak:.4f} >= 0.999)")
    if peak < 1e-4:
        issues.append(f"{path.name}: near-silent peak={peak:.3g}")
    return issues


def _check_provenance(metrics: dict, issues: list[str], feedback: list[dict]) -> None:
    eps = REVIEW_CONFIG["provenance_eps"]
    targets = [float(t) for t in (metrics.get("formant_targets_hz") or [])]
    legacy_peaks = metrics.get("formant_peaks_hz")
    if legacy_peaks is not None:
        issues.append("deprecated field 'formant_peaks_hz' present at root (must be removed/migrated)")
        if targets and len(legacy_peaks) == len(targets) and all(
            abs(float(p) - t) <= eps for p, t in zip(legacy_peaks, targets)
        ):
            msg = "circular metric: formant_peaks_hz is an exact copy of formant_targets_hz"
            issues.append(msg)
            feedback.append(
                {
                    "issue": msg,
                    "observed": legacy_peaks,
                    "expected": "independent H(f)/LPC measurements in formants_measured",
                    "suggested_action": "Remove formant_peaks_hz copy and measure WAVs via pipeline.measure",
                }
            )

    fm = metrics.get("formants_measured") or {}
    tf_peaks = [float(p["f_hz"]) for p in (fm.get("transfer_function") or []) if "f_hz" in p]
    lpc_peaks = [float(v) for v in ((fm.get("lpc") or {}).get("f_hz") or [])]
    for src_name, meas in (("transfer_function", tf_peaks), ("lpc", lpc_peaks)):
        for m_val in meas:
            for t_val in targets:
                if abs(m_val - t_val) <= eps:
                    msg = f"circular metric: {src_name} measured peak {m_val} Hz equals target {t_val} Hz within {eps}"
                    issues.append(msg)
                    feedback.append(
                        {
                            "issue": msg,
                            "observed": m_val,
                            "expected": f"independent signal measurement != {t_val}",
                            "suggested_action": "Verify pipeline.measure never reads formant_targets_hz",
                        }
                    )


def _check_kl_q_errors(
    metrics: dict,
    issues: list[str],
    info_notes: list[str],
    feedback: list[dict],
) -> None:
    kl_block = metrics.get("kl_vs_tmm") or {}
    q_errs = kl_block.get("q_rel_err") or []
    loss_model = kl_block.get("loss_model", "legacy")
    tol = REVIEW_CONFIG["kl_vs_tmm_max_q_rel_err"]
    for idx, q_err in enumerate(q_errs[:3]):
        val = float(q_err)
        if val > tol:
            if loss_model == "legacy":
                info_notes.append(
                    f"KL(legacy) F{idx + 1} Q deviates {val * 100:.1f}% > {tol * 100:.0f}% from TMM (expected in legacy mode)."
                )
            else:
                msg = (
                    f"KL vs TMM F{idx + 1} Q relative error {val * 100:.1f}% exceeds {tol * 100:.0f}% "
                    f"(loss_model={loss_model})"
                )
                issues.append(msg)
                feedback.append(
                    {
                        "issue": msg,
                        "observed": round(val, 4),
                        "expected": f"<= {tol}",
                        "suggested_action": "Verify Kirchhoff biquad loss fit and radiation reflection filter in source_tract.py",
                    }
                )


def review_metrics_only(metrics: dict) -> dict:
    """Evaluate a standalone metrics.json dictionary (e.g. fixture test without WAVs)."""
    if metrics.get("synthetic") is True:
        raise ValueError("Refusing to peer-review a synthetic what-if metrics payload (synthetic=true).")
    issues: list[str] = []
    info_notes: list[str] = []
    feedback: list[dict] = []
    _check_provenance(metrics, issues, feedback)
    _check_kl_q_errors(metrics, issues, info_notes, feedback)
    if metrics.get("synth_mode") == "tmm_fir" or metrics.get("independent_of_tmm") is False:
        issues.append(
            "synth_mode='tmm_fir' uses the TMM transfer function itself (independent_of_tmm=false); refusing circular validation."
        )

    if metrics.get("clipping_detected") or (metrics.get("audio") or {}).get("clipping_detected"):
        issues.append("clipping_detected is true")

    if "formants_measured" not in metrics or "formants_predicted_tmm" not in metrics:
        issues.append("missing v2 measurement/prediction blocks (formants_measured / formants_predicted_tmm)")
        feedback.append(
            {
                "issue": "legacy circular metrics schema without independent M1/M2 blocks",
                "observed": sorted(metrics.keys()),
                "expected": ["formants_measured", "formants_predicted_tmm", "shaping_gain_db_F1"],
                "suggested_action": "Run pipeline.measure (M1) and pipeline.tmm_lossy (M2)",
            }
        )
    else:
        sh_block = metrics.get("shaping") or {}
        sg_b3db = sh_block.get("shaping_gain_db_F1_b3db")
        sg = sg_b3db if sg_b3db is not None else metrics.get("shaping_gain_db_F1")
        if sg is None or float(sg) < REVIEW_CONFIG["min_shaping_gain_db"]:
            issues.append("cavity filter does not shape the spectrum")
            feedback.append(
                {
                    "issue": "cavity filter does not shape the spectrum",
                    "observed": sg,
                    "expected": f">= {REVIEW_CONFIG['min_shaping_gain_db']} dB",
                    "suggested_action": "increase source spectral richness (M4)",
                }
            )
        e_pm60 = sh_block.get("f1_band_energy_fraction_measured_pm60")
        if e_pm60 is not None and float(e_pm60) < REVIEW_CONFIG["min_f1_band_energy_fraction"]:
            issues.append(
                f"F1 band energy fraction in filtered.wav too low "
                f"({float(e_pm60) * 100:.3f}% < {REVIEW_CONFIG['min_f1_band_energy_fraction'] * 100:.1f}%)"
            )

    passed = len(issues) == 0
    geom_status = (metrics.get("geometry_qc") or {}).get("status")
    return {
        "state": "PEER_REVIEW",
        "terminal_state": resolve_terminal_state(passed, geom_status),
        "run_id": metrics.get("run_id", "unknown"),
        "passed": passed,
        "issues": issues,
        "tolerances": dict(REVIEW_CONFIG),
        "retry_feedback": feedback,
    }


def review_run(run_id_or_dir: str | Path, *, profile_csv: Path | None = None) -> dict:
    """Run independent 3-way physics peer review on a run directory (or fail if synthetic metrics file)."""
    run_path = Path(run_id_or_dir)
    if run_path.is_file():
        synth_data = _load_json(run_path)
        if synth_data.get("synthetic") is True:
            raise ValueError(f"Refusing to peer-review synthetic metrics file '{run_path}' (synthetic=true).")
        run_dir = run_path.parent
    else:
        run_dir = run_path if run_path.exists() else (OUTPUT / str(run_id_or_dir))
    run_id = run_dir.name

    issues: list[str] = []
    info_notes: list[str] = []
    feedback: list[dict] = []

    if run_dir.exists() and not (run_dir / "metrics.json").exists():
        for cand in sorted(run_dir.glob("SYNTHETIC_*_metrics.json")):
            if _load_json(cand).get("synthetic") is True:
                raise ValueError(f"Refusing to peer-review synthetic directory '{run_dir}' ({cand.name}).")

    metrics = _load_json(run_dir / "metrics.json") if (run_dir / "metrics.json").exists() else {}
    if metrics.get("synthetic") is True:
        raise ValueError(f"Refusing to peer-review synthetic metrics in '{run_dir}' (synthetic=true).")

    missing = [name for name in REQUIRED if not (run_dir / name).exists()]
    if missing:
        issues.append(f"missing artifacts: {missing}")
        feedback.append(
            {
                "issue": f"missing required artifacts: {missing}",
                "observed": missing,
                "expected": list(REQUIRED),
                "suggested_action": "Run pipeline.build_and_sim to produce all required artifacts",
            }
        )

    if metrics:
        _check_provenance(metrics, issues, feedback)
        _check_kl_q_errors(metrics, issues, info_notes, feedback)
        if metrics.get("synth_mode") == "tmm_fir" or metrics.get("independent_of_tmm") is False:
            msg = "synth_mode='tmm_fir' uses the TMM transfer function itself (independent_of_tmm=false); refusing circular measured-vs-TMM validation."
            issues.append(msg)
            feedback.append(
                {
                    "issue": msg,
                    "observed": metrics.get("synth_mode"),
                    "expected": "independent time-domain waveguide synthesis ('kl')",
                    "suggested_action": "Use --synth kl (3-port wall junction) instead of --synth tmm_fir for reviewable runs.",
                }
            )
        if metrics.get("clipping_detected") or (metrics.get("audio") or {}).get("clipping_detected"):
            issues.append("metrics reports clipping_detected=true")

    for name in ("source.wav", "filtered.wav"):
        if (run_dir / name).exists():
            issues.extend(_wav_issues(run_dir / name))

    # Measure WAVs independently if source.wav & filtered.wav exist
    measured: dict | None = None
    if (run_dir / "measured.json").exists() and metrics.get("schema_version") in ("2.0", "2.1"):
        measured = _load_json(run_dir / "measured.json")
    elif (run_dir / "source.wav").exists() and (run_dir / "filtered.wav").exists():
        measured = measure_run_dir(run_dir, write_outputs=False)

    # Obtain independent TMM predictions
    tmm_summary: dict | None = None
    if (run_dir / "tmm_lossy.json").exists():
        tmm_summary = _load_json(run_dir / "tmm_lossy.json")
    elif metrics.get("formants_predicted_tmm"):
        tmm_summary = {"peaks": metrics["formants_predicted_tmm"]}
    else:
        csv_fallback = profile_csv or (ROOT / "airway_features.csv")
        if csv_fallback.exists():
            tmm_summary = run_tmm_profile(csv_fallback, run_dir, write_outputs=False)["summary"]

    comparisons: list[dict] = []
    sh_metrics: dict | None = None
    shaping_gain_db = None
    f1_band_frac = None

    if measured is not None and tmm_summary is not None:
        tf_peaks = measured["formants_measured"]["transfer_function"]
        lpc_peaks = measured["formants_measured"]["lpc"]["f_hz"]
        meas_source = "transfer_function" if tf_peaks else "lpc"
        meas_freqs = [float(p["f_hz"]) for p in tf_peaks] if tf_peaks else [float(f) for f in lpc_peaks]

        tmm_peaks = [p for p in tmm_summary["peaks"] if p["f_hz"] <= REVIEW_CONFIG["max_compare_freq_hz"]]
        if not tmm_peaks:
            issues.append("no TMM predicted peaks below 2000 Hz")
        elif not meas_freqs:
            issues.append("no measured formant peaks detected from audio")
            feedback.append(
                {
                    "issue": "no measured formant peaks detected from audio",
                    "observed": [],
                    "expected": [p["f_hz"] for p in tmm_peaks],
                    "suggested_action": "increase source spectral richness and coherence across 100-2000 Hz (M4)",
                }
            )
        else:
            for idx, t_pk in enumerate(tmm_peaks):
                f_pred = float(t_pk["f_hz"])
                nearest_f = min(meas_freqs, key=lambda mf: abs(mf - f_pred))
                rel_err = abs(nearest_f - f_pred) / max(f_pred, 1e-9)
                tol = (
                    REVIEW_CONFIG["f1_measured_vs_tmm_max_rel_err"]
                    if idx == 0
                    else REVIEW_CONFIG["f2_f3_measured_vs_tmm_max_rel_err"]
                )
                ok = rel_err <= tol
                comparisons.append(
                    {
                        "formant": f"F{idx + 1}",
                        "measured_hz": round(nearest_f, 3),
                        "measured_source": meas_source,
                        "predicted_tmm_hz": round(f_pred, 3),
                        "rel_err": round(rel_err, 5),
                        "tolerance": tol,
                        "passed": ok,
                    }
                )
                if not ok:
                    msg = (
                        f"F{idx + 1} measured ({nearest_f:.1f} Hz via {meas_source}) deviates "
                        f"{rel_err * 100:.2f}% from TMM ({f_pred:.1f} Hz, tol={tol * 100:.0f}%)"
                    )
                    issues.append(msg)
                    feedback.append(
                        {
                            "issue": msg,
                            "observed": round(nearest_f, 2),
                            "expected": round(f_pred, 2),
                            "suggested_action": "use coupled Kelly-Lochbaum waveguide built from TMM profile (M4)",
                        }
                    )

            # TMM vs analytic check (closed-open with end correction L_eff, coupled with wall mass if active)
            analytic_block = metrics.get("analytic") or metrics.get("theoretical") or {}
            analytic_eff = tmm_summary.get("analytic_closed_open_eff_hz") or analytic_block.get("closed_open_hz") or []
            if analytic_eff:
                f1_tmm = float(tmm_peaks[0]["f_hz"])
                f1_ana_rigid = float(analytic_eff[0])
                fw_closed = tmm_summary.get("wall_resonance_closed_hz")
                if fw_closed is not None and float(fw_closed) > 0.0:
                    f1_ana = float(np.sqrt(f1_ana_rigid**2 + float(fw_closed) ** 2))
                    ana_label = f"coupled wall-air analytic sqrt((c/4L_eff)^2 + f_w^2) ({f1_ana:.1f} Hz)"
                else:
                    f1_ana = f1_ana_rigid
                    ana_label = f"analytic ({f1_ana:.1f} Hz)"
                ana_err = abs(f1_tmm - f1_ana) / max(f1_ana, 1e-9)
                if ana_err > REVIEW_CONFIG["tmm_vs_analytic_fail_rel_err"]:
                    issues.append(
                        f"TMM F1 ({f1_tmm:.1f} Hz) deviates {ana_err * 100:.1f}% > 30% from {ana_label}"
                    )
                elif ana_err > REVIEW_CONFIG["tmm_vs_analytic_info_rel_err"]:
                    info_notes.append(
                        f"TMM F1 ({f1_tmm:.1f} Hz) differs {ana_err * 100:.1f}% from {ana_label} "
                        "due to cross-sectional area flaring."
                    )

        # Spectral-shaping check on actual WAVs around TMM F1
        if tmm_peaks and (run_dir / "source.wav").exists() and (run_dir / "filtered.wav").exists():
            _, src_x, _ = read_wav_float(run_dir / "source.wav")
            sr_f, flt_x, _ = read_wav_float(run_dir / "filtered.wav")
            f1_tmm = float(tmm_peaks[0]["f_hz"])
            b1_tmm = float(tmm_peaks[0].get("bw_3db_hz", 70.0))
            sh_metrics = compute_shaping_metrics(src_x, flt_x, sr_f, f1_tmm, b1_tmm)
            shaping_gain_db = float(sh_metrics["shaping_gain_db_F1_b3db"])
            f1_band_frac = float(sh_metrics["f1_band_energy_fraction_measured_pm60"])
            if shaping_gain_db < REVIEW_CONFIG["min_shaping_gain_db"]:
                msg = "cavity filter does not shape the spectrum"
                issues.append(f"{msg} (shaping_gain_db_F1_b3db={shaping_gain_db:.2f} dB < {REVIEW_CONFIG['min_shaping_gain_db']} dB)")
                feedback.append(
                    {
                        "issue": msg,
                        "observed": round(shaping_gain_db, 3),
                        "expected": f">= {REVIEW_CONFIG['min_shaping_gain_db']} dB",
                        "suggested_action": "increase source spectral richness (M4)",
                    }
                )
            if f1_band_frac < REVIEW_CONFIG["min_f1_band_energy_fraction"]:
                msg = (
                    f"F1 band energy fraction in filtered.wav too low "
                    f"({f1_band_frac * 100:.3f}% < {REVIEW_CONFIG['min_f1_band_energy_fraction'] * 100:.1f}%)"
                )
                issues.append(msg)
                feedback.append(
                    {
                        "issue": "insufficient energy in F1 band of filtered.wav",
                        "observed": round(f1_band_frac, 6),
                        "expected": f">= {REVIEW_CONFIG['min_f1_band_energy_fraction']}",
                        "suggested_action": "increase source spectral richness (M4) and match KL Q-factor to TMM",
                    }
                )

    passed = len(issues) == 0
    geom_status = (metrics.get("geometry_qc") or {}).get("status")
    terminal_state = resolve_terminal_state(passed, geom_status)

    report = {
        "state": "PEER_REVIEW",
        "terminal_state": terminal_state,
        "contract_status": "CONTRACT_OK" if passed else "REVIEW_FAILED",
        "run_id": run_id,
        "passed": passed,
        "issues": issues,
        "info_notes": info_notes,
        "tolerances": dict(REVIEW_CONFIG),
        "comparisons": comparisons,
        "shaping_check": {
            "shaping_gain_db_F1": None if sh_metrics is None else sh_metrics["shaping_gain_db_F1_pm60"],
            "shaping_gain_db_F1_b3db": None if sh_metrics is None else sh_metrics["shaping_gain_db_F1_b3db"],
            "shaping_gain_db_F1_pm60": None if sh_metrics is None else sh_metrics["shaping_gain_db_F1_pm60"],
            "filtered_f1_band_fraction": None if f1_band_frac is None else round(f1_band_frac, 6),
            "f1_band_energy_fraction_measured_pm60": None if f1_band_frac is None else round(f1_band_frac, 6),
            "b3db_band_hz": None if sh_metrics is None else sh_metrics["b3db"]["band_hz"],
            "pm60_band_hz": None if sh_metrics is None else sh_metrics["pm60"]["band_hz"],
        },
        "retry_feedback": feedback,
    }
    if run_dir.exists():
        (run_dir / "peer_review.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return report

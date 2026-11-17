"""État 2 — BUILD_AND_SIM: geometry QC, lossy TMM, KL/legacy synth, M1 measure, M3 review, report."""

from __future__ import annotations

import argparse
import importlib.metadata
import importlib.util
import json
import logging
import subprocess
import sys
import time
import uuid
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.io import wavfile

from pipeline.geometry_qc import (
    GeometryInputError,
    evaluate_geometry_qc,
    extract_airway_from_stl,
    plot_geometry_qc,
)
from pipeline.hypothesis import load_metadata, write_hypothesis_bundle
from pipeline.measure import compute_shaping_gain_db, measure_run_dir
from pipeline.paths import OUTPUT, ROOT
from pipeline.peer_review import review_run
from pipeline.source_tract import (
    build_kl_waveguide,
    compare_kl_vs_tmm,
    kl_impulse_response,
    synthesize_kl,
    synthesize_legacy_iir,
)
from pipeline.tmm_lossy import Air, run_tmm_profile

logger = logging.getLogger(__name__)


def _new_run_id(taxon_slug: str) -> str:
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    return f"{taxon_slug}_{stamp}_{uuid.uuid4().hex[:6]}"


def _git_commit() -> str:
    for git_bin in ("git", r"C:\Program Files\Git\cmd\git.exe"):
        try:
            proc = subprocess.run(
                [git_bin, "rev-parse", "--short", "HEAD"],
                cwd=ROOT,
                capture_output=True,
                text=True,
                timeout=3,
                check=True,
            )
            return proc.stdout.strip()
        except Exception:
            continue
    return "unknown"


def _package_versions(warnings_out: list[str]) -> dict[str, str | None]:
    pkgs = ("numpy", "scipy", "matplotlib", "soundfile", "trimesh", "jsonschema", "pytest", "vmtk", "gmsh", "dolfinx")
    versions: dict[str, str | None] = {}
    for pkg in pkgs:
        spec = importlib.util.find_spec(pkg)
        if spec is None:
            versions[pkg] = None
            if pkg in ("soundfile", "trimesh"):
                msg = f"Optional dependency '{pkg}' is not installed."
                logger.warning(msg)
                warnings_out.append(msg)
        else:
            try:
                versions[pkg] = importlib.metadata.version(pkg)
            except importlib.metadata.PackageNotFoundError:
                versions[pkg] = "installed"
    return versions


def _write_wav(path: Path, audio: np.ndarray, sr: int) -> None:
    arr = np.asarray(audio, dtype=np.float64)
    assert float(np.max(np.abs(arr))) < 0.999, f"Refusing to write clipped WAV to {path}"
    pcm = np.int16(np.round(np.clip(arr, -1.0, 1.0) * 32767.0))
    wavfile.write(path, sr, pcm)


def _spectrogram(path_png: Path, audio: np.ndarray, sr: int) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.specgram(audio, NFFT=1024, Fs=sr, noverlap=512, cmap="magma")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Frequency (Hz)")
    ax.set_ylim(0, 4000)
    ax.set_title("Filtered call spectrogram")
    fig.tight_layout()
    fig.savefig(path_png, dpi=140)
    plt.close(fig)


def _validate_metrics_schema(metrics: dict, warnings_out: list[str]) -> None:
    schema_path = ROOT / "schemas" / "metrics.schema.json"
    schema = json.loads(schema_path.read_text(encoding="utf-8"))
    try:
        import jsonschema

        jsonschema.validate(instance=metrics, schema=schema)
    except ImportError:
        msg = "jsonschema not installed; fell back to required-key structural check."
        logger.warning(msg)
        warnings_out.append(msg)
        for req in schema["required"]:
            if req not in metrics:
                raise ValueError(f"metrics.json missing required key '{req}'")


def _render_report_md(metrics: dict, kl_cmp: dict | None, review: dict) -> str:
    state = metrics["state"]
    qc = metrics["geometry_qc"]
    tf_peaks = metrics["formants_measured"]["transfer_function"]
    lpc = metrics["formants_measured"]["lpc"]
    tmm_peaks = metrics["formants_predicted_tmm"]
    ana = metrics["analytic"]
    why = "; ".join(review["issues"]) if review["issues"] else "All physics, spectral-shaping, and schema checks passed."

    lines: list[str] = [
        f"# Paleo-Acoustic Run Report — `{metrics['run_id']}`",
        "",
        "## 1. Summary",
        "",
        "| Field | Value |",
        "| :--- | :--- |",
        f"| Taxon / Specimen | {metrics['taxon']} ({metrics['specimen_id']}) |",
        f"| Synth Mode | `{metrics['synth_mode']}` |",
        f"| Terminal State | **`{state}`** |",
        f"| Why | {why} |",
        f"| Measured f0 (YIN) | {metrics['f0']['measured_hz']:.2f} ± {metrics['f0']['measured_std_hz']:.2f} Hz (prior: {metrics['f0']['prior_hz']} Hz) |",
        f"| Shaping Gain (F1 band) | {metrics['shaping_gain_db_F1']:.2f} dB |",
        f"| Source Spectral Tilt | {metrics['source_spectral_tilt_db_per_oct']:.2f} dB/oct (target: -15 to -9 dB/oct) |",
        "",
        "## 2. Geometry QC",
        "",
    ]

    if qc["status"] == "BONE_DERIVED_EXPLORATORY":
        lines.extend(
            [
                "> [!WARNING]",
                "> **BONE-DERIVED GEOMETRY — EXPLORATORY ONLY.** This profile comes from the bone mesh skeleton rather than a segmented nasal airway lumen. Cross-sectional areas (`A_mean = "
                + f"{qc['A_mean_mm2']:.1f} mm²`, `A_max = {qc['A_max_mm2']:.1f} mm²`, `Dh_mean = {qc['Dh_mean_mm']:.2f} mm`) and single-pass length (`L = {qc['L_mm']:.1f} mm`) do not represent a looped soft-tissue airway. This run cannot be `SCIENTIFICALLY_VALIDATED`.",
                "",
            ]
        )

    lines.extend(
        [
            f"- **Status:** `{qc['status']}`",
            f"- **Length:** {qc['L_mm']:.3f} mm | **Mean Area:** {qc['A_mean_mm2']:.2f} mm² | **Max Area:** {qc['A_max_mm2']:.2f} mm² | **Mean Dh:** {qc['Dh_mean_mm']:.2f} mm | **Sections:** {qc['n_sections']}",
            f"- **Source Mesh:** `{qc['source_mesh']}` (`{qc['extraction_method']}`)",
            f"- **Flags:** {qc['flags'] if qc['flags'] else 'None'}",
            "",
            "![Geometry QC](geometry_qc.png)",
            "",
            "## 3. Formants: Measured vs Predicted (Lossy TMM) vs Analytic",
            "",
            "| Mode | Measured H(f) (Hz) | Measured LPC (Hz) | Predicted Lossy TMM (Hz) | Analytic c/4L (Hz) | Analytic c/4L_eff (Hz) |",
            "| :--- | :--- | :--- | :--- | :--- | :--- |",
        ]
    )

    for idx in range(3):
        m_tf = f"{tf_peaks[idx]['f_hz']:.2f} (Q={tf_peaks[idx]['q']:.1f}, γ²={tf_peaks[idx]['coherence']:.3f})" if idx < len(tf_peaks) else "—"
        m_lpc = f"{lpc['f_hz'][idx]:.2f}" if idx < len(lpc["f_hz"]) else "—"
        p_tmm = f"{tmm_peaks[idx]['f_hz']:.2f} (Q={tmm_peaks[idx]['q']:.1f})" if idx < len(tmm_peaks) else "—"
        c_raw = f"{ana['closed_open_hz'][idx]:.2f}" if idx < len(ana["closed_open_hz"]) else "—"
        c_eff = f"{((2*idx+1)*ana['c_m_s']/(4.0*ana['L_eff_m'])):.2f}"
        lines.append(f"| F{idx+1} | {m_tf} | {m_lpc} | {p_tmm} | {c_raw} | {c_eff} |")

    lines.extend(
        [
            "",
            f"- **LPC Reliability:** `{lpc['reliability']}` ({lpc.get('reliability_reasons', [])})",
            f"- **Why TMM F1 differs from c/4L:** Unflanged Levine–Schwinger end correction lengthens the effective tube from $L = {ana['L_m']:.4f}\\text{{ m}}$ to $L_{{\\text{{eff}}}} = {ana['L_eff_m']:.4f}\\text{{ m}}$ (lowering uniform $c/4L_{{\\text{{eff}}}}$ to `{((ana['c_m_s'])/(4.0*ana['L_eff_m'])):.2f}` Hz), whereas strong cross-sectional flaring from `{qc['A_mean_mm2']:.0f}` mm² average (up to `{qc['A_max_mm2']:.0f}` mm²) acts like an expanding horn that pushes the fundamental resonance upward to `{tmm_peaks[0]['f_hz']:.2f}` Hz.",
            "",
            "## 4. Measured Transfer Function H(f)",
            "",
            "![Transfer Function](transfer_function.png)",
            "",
            "## 5. Lossy Kirchhoff TMM Spectrum",
            "",
            "![Lossy TMM](tmm_lossy.png)",
            "",
            "## 6. Kelly–Lochbaum Waveguide vs Lossy TMM",
            "",
        ]
    )

    if kl_cmp and kl_cmp.get("comparisons"):
        lines.extend(
            [
                "| Formant | TMM (Hz) | KL Impulse (Hz) | Rel Error | TMM Q | KL Q |",
                "| :--- | :--- | :--- | :--- | :--- | :--- |",
            ]
        )
        for item in kl_cmp["comparisons"]:
            lines.append(
                f"| — | {item['tmm_hz']:.2f} | {item['kl_hz']:.2f} | {item['rel_err']*100:.3f}% | {item['tmm_q']:.1f} | {item['kl_q']:.1f} |"
            )
        lines.append("")
    lines.extend(
        [
            "![KL vs TMM](kl_vs_tmm.png)",
            "",
            "## 7. Band Energy Fractions (Filtered vs Source)",
            "",
            "| Band (Hz) | Filtered Fraction | Source Fraction |",
            "| :--- | :--- | :--- |",
        ]
    )

    flt_bands = metrics["band_energy_fractions"]["filtered"]
    src_bands = metrics["band_energy_fractions"]["source"]
    for band_k, f_val in flt_bands.items():
        s_val = src_bands.get(band_k, 0.0)
        lines.append(f"| {band_k} | {f_val*100:.3f}% | {s_val*100:.3f}% |")

    lines.extend(
        [
            "",
            "## 8. Source Spectral Tilt",
            "",
            f"- Measured harmonic tilt (100–2000 Hz): **{metrics['source_spectral_tilt_db_per_oct']:.2f} dB/octave**",
            "",
            "## 9. Assumptions",
            "",
            f"```json\n{json.dumps(metrics['assumptions'], indent=2)}\n```",
            "",
            "## 10. Warnings",
            "",
            "\n".join(f"- {w}" for w in metrics["warnings"]) if metrics["warnings"] else "- None",
            "",
            "## 11. Limits",
            "",
            "- **Bone-derived geometry:** Unless a validated CT lumen or endocast profile with `.source.json` is supplied, the cross-sectional profile reflects bone wall geometry rather than the true looped nasal passage.",
            "- **Unknown larynx / syrinx:** No fossilized laryngeal soft tissues are preserved; the two-mass vocal-fold model and subglottal pressure are biomechanical hypotheses.",
            "- **Soft tissue absent:** Mucosal lining, wall compliance, and nasal valve geometry are not directly observable in the fossil STL.",
            f"- **f0 prior:** `f0_prior_hz = {metrics['f0']['prior_hz']} Hz` is a soft-tissue assumption, not an anatomical measurement.",
            "",
            "## 12. Artifacts",
            "",
        ]
    )
    for k, v in metrics.get("artifacts", {}).items():
        lines.append(f"- `{k}`: `{v}`")

    return "\n".join(lines) + "\n"


def run_build_and_sim(
    taxon_slug: str,
    run_id: str | None = None,
    *,
    synth_mode: str = "kl",
    f0_hz: float | None = None,
    temperature_c: float | None = None,
    radiation: str = "unflanged",
    seed: int = 0,
    duration_s: float = 2.5,
    sample_rate: int = 44100,
    fmax: float = 4000.0,
    airway_profile: Path | None = None,
    airway_stl: Path | None = None,
    inlet: tuple[float, float, float] | None = None,
    outlet: tuple[float, float, float] | None = None,
    cli_command: str | None = None,
) -> dict:
    t_start = time.perf_counter()
    meta = load_metadata(taxon_slug)
    run_id = run_id or _new_run_id(taxon_slug)
    run_dir = OUTPUT / run_id
    run_dir.mkdir(parents=True, exist_ok=True)

    warnings_list: list[str] = []
    pkg_versions = _package_versions(warnings_list)

    # Write HYPOTHESIS bundle & state_HYPOTHESIS.json if not already emitted
    if not (run_dir / "antigravity_prompt.md").exists():
        hyp_state = write_hypothesis_bundle(taxon_slug, run_id)
    else:
        hyp_state = {"state": "HYPOTHESIS", "run_id": run_id, "taxon": taxon_slug}
    (run_dir / "state_HYPOTHESIS.json").write_text(json.dumps(hyp_state, indent=2) + "\n", encoding="utf-8")

    # Resolve geometry input (M5)
    if airway_stl is not None:
        extracted_csv = run_dir / "extracted_airway.csv"
        profile_path = extract_airway_from_stl(Path(airway_stl), inlet, outlet, extracted_csv)
        qc, prof = evaluate_geometry_qc(
            profile_path,
            source_mesh=str(airway_stl),
            extraction_method="lumen_stl_centerline",
            require_sidecar=False,
        )
    elif airway_profile is not None:
        cand = Path(airway_profile)
        profile_path = cand if cand.exists() else (ROOT / cand)
        qc, prof = evaluate_geometry_qc(profile_path, require_sidecar=True)
    else:
        rel_csv = meta["geometry"]["airway_features_csv"]
        profile_path = ROOT / rel_csv
        qc, prof = evaluate_geometry_qc(
            profile_path,
            source_mesh=meta["geometry"]["stl_path"],
            extraction_method="bone_mesh_skeleton_axis_slice",
            require_sidecar=False,
        )

    geom_png = run_dir / "geometry_qc.png"
    plot_geometry_qc(qc, prof, geom_png)

    # Air & soft-tissue priors
    priors = meta["acoustic_priors"]
    soft = meta["soft_tissue_assumptions"]
    temp_used = float(soft.get("temperature_c", 20.0)) if temperature_c is None else float(temperature_c)
    f0_prior = float(priors.get("f0_prior_hz", 60.0)) if f0_hz is None else float(f0_hz)
    air = Air.from_temperature(temp_used)

    # M2: Lossy TMM
    tmm_out = run_tmm_profile(
        profile_path,
        run_dir,
        fmax=fmax,
        df=0.5,
        temperature_c=temp_used,
        radiation=radiation,
        write_outputs=True,
    )
    tmm_summary = tmm_out["summary"]
    tmm_raw = tmm_out["raw"]
    warnings_list.extend(tmm_summary.get("warnings", []))

    # M4: Synthesis (--synth kl or --synth legacy_iir)
    kl_png = run_dir / "kl_vs_tmm.png"
    if synth_mode == "kl":
        synth = synthesize_kl(
            prof["s_mm"],
            prof["area_mm2"],
            air,
            tmm_summary,
            duration_s=duration_s,
            sr=sample_rate,
            f0_hz=f0_prior,
            radiation=radiation,
            seed=seed,
            kl_plot_path=kl_png,
            tmm_raw=tmm_raw,
        )
        kl_cmp = synth["kl_check"]
    elif synth_mode == "legacy_iir":
        synth = synthesize_legacy_iir(
            list(priors["expected_formants_hz"]),
            list(priors["default_q_factors"]),
            duration_s=duration_s,
            sr=sample_rate,
            f0_hz=f0_prior,
            seed=seed,
        )
        # Still run KL IR check for diagnostics & plot so artifact contract is complete
        wg = build_kl_waveguide(prof["s_mm"], prof["area_mm2"], air, radiation=radiation, sr=sample_rate)
        ir_wg = kl_impulse_response(wg, ir_duration_s=0.8)
        kl_cmp = compare_kl_vs_tmm(wg, ir_wg, tmm_summary, out_png=kl_png, tmm_raw=tmm_raw)
    else:
        raise ValueError(f"Unsupported synth mode: {synth_mode}")

    source_wav = run_dir / "source.wav"
    filtered_wav = run_dir / "filtered.wav"
    _write_wav(source_wav, synth["source"], sample_rate)
    _write_wav(filtered_wav, synth["filtered"], sample_rate)

    spec_png = run_dir / "spectrogram.png"
    _spectrogram(spec_png, synth["filtered"], sample_rate)

    (run_dir / "state_BUILD_AND_SIM.json").write_text(
        json.dumps({"state": "BUILD_AND_SIM", "run_id": run_id, "synth_mode": synth_mode}, indent=2) + "\n",
        encoding="utf-8",
    )

    # M1: Audio-only measurement
    measured = measure_run_dir(run_dir, write_outputs=True)

    # Compute F1 shaping gain around TMM F1
    f1_tmm_pk = tmm_summary["peaks"][0]
    shaping_db, _, _ = compute_shaping_gain_db(
        synth["source"],
        synth["filtered"],
        sample_rate,
        float(f1_tmm_pk["f_hz"]),
        max(float(f1_tmm_pk["bw_3db_hz"]), 120.0),
    )

    artifacts_map = {
        "source.wav": str(source_wav),
        "filtered.wav": str(filtered_wav),
        "spectrogram.png": str(spec_png),
        "transfer_function.png": str(run_dir / "transfer_function.png"),
        "tmm_lossy.png": str(run_dir / "tmm_lossy.png"),
        "kl_vs_tmm.png": str(kl_png),
        "geometry_qc.png": str(geom_png),
        "measured.json": str(run_dir / "measured.json"),
        "tmm_lossy.json": str(run_dir / "tmm_lossy.json"),
        "report.md": str(run_dir / "report.md"),
        "build_log.txt": str(run_dir / "build_log.txt"),
    }

    metrics = {
        "schema_version": "2.0",
        "taxon": meta["taxon"],
        "specimen_id": meta["specimen_id"],
        "taxon_slug": taxon_slug,
        "run_id": run_id,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "git_commit": _git_commit(),
        "seed": int(seed),
        "sample_rate": int(sample_rate),
        "duration_s": float(duration_s),
        "synth_mode": synth_mode,
        "state": "PEER_REVIEW",
        "f0": {
            "prior_hz": f0_prior,
            "measured_hz": measured["f0"]["measured_hz"],
            "measured_std_hz": measured["f0"]["measured_std_hz"],
            "method": measured["f0"]["method"],
        },
        "formant_targets_hz": [float(x) for x in priors["expected_formants_hz"]],
        "formants_measured": measured["formants_measured"],
        "formants_predicted_tmm": [
            {"f_hz": p["f_hz"], "bw_3db_hz": p["bw_3db_hz"], "q": p["q"]}
            for p in tmm_summary["peaks"][:4]
        ],
        "analytic": {
            "L_m": tmm_summary["L_m"],
            "L_eff_m": tmm_summary["L_eff_m"],
            "c_m_s": air.c,
            "closed_open_hz": tmm_summary["analytic_closed_open_hz"],
            "open_open_hz": tmm_summary["analytic_open_open_hz"],
        },
        "band_energy_fractions": measured["band_energy_fractions"],
        "shaping_gain_db_F1": round(float(shaping_db), 3),
        "source_spectral_tilt_db_per_oct": measured["source_spectral_tilt_db_per_oct"],
        "kl_vs_tmm_rel_err": kl_cmp["rel_errors"] if kl_cmp else [],
        "geometry_qc": qc,
        "review": {
            "passed": False,
            "issues": [],
            "tolerances": {},
            "retry_feedback": [],
        },
        "assumptions": {
            "air": tmm_summary["air"],
            "soft_tissue": {
                **soft,
                "f0_prior_hz": f0_prior,
                "source_model": synth["source_description"],
                "tract_model": synth["filtered_description"],
            },
            "radiation_model": radiation,
        },
        "audio": {
            "peak_dbfs": round(max(synth["audio"]["source_peak_dbfs"], synth["audio"]["filtered_peak_dbfs"]), 3),
            "source_peak_dbfs": synth["audio"]["source_peak_dbfs"],
            "filtered_peak_dbfs": synth["audio"]["filtered_peak_dbfs"],
            "clipping_detected": bool(synth["audio"]["clipping_detected"] or measured["audio"]["clipping_detected"]),
        },
        "warnings": warnings_list,
        "deprecated": {
            "formant_peaks_hz": "removed: was a copy of filter targets",
        },
        "artifacts": artifacts_map,
    }

    metrics_path = run_dir / "metrics.json"
    metrics_path.write_text(json.dumps(metrics, indent=2) + "\n", encoding="utf-8")

    # M3: Non-circular peer review
    review = review_run(run_dir, profile_csv=profile_path)
    metrics["review"] = {
        "passed": review["passed"],
        "issues": review["issues"],
        "tolerances": review["tolerances"],
        "retry_feedback": review["retry_feedback"],
    }
    metrics["state"] = review["terminal_state"]

    _validate_metrics_schema(metrics, warnings_list)
    metrics["warnings"] = warnings_list
    metrics_path.write_text(json.dumps(metrics, indent=2) + "\n", encoding="utf-8")

    # State JSONs & Report & Build log
    (run_dir / f"state_{metrics['state']}.json").write_text(
        json.dumps(
            {
                "state": metrics["state"],
                "passed": review["passed"],
                "contract_status": review["contract_status"],
                "geometry_status": qc["status"],
                "issues": review["issues"],
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    report_md = _render_report_md(metrics, kl_cmp, review)
    (run_dir / "report.md").write_text(report_md, encoding="utf-8")

    elapsed = time.perf_counter() - t_start
    cmd_str = cli_command or f"python -m pipeline.build_and_sim --taxon {taxon_slug} --run-id {run_id} --synth {synth_mode} --seed {seed}"
    build_log = (
        f"command: {cmd_str}\n"
        f"run_id: {run_id}\n"
        f"duration_wall_s: {elapsed:.3f}\n"
        f"final_state: {metrics['state']}\n"
        f"geometry_qc_status: {qc['status']}\n"
        f"python: {sys.version.splitlines()[0]}\n"
        f"packages: {json.dumps(pkg_versions)}\n"
        f"warnings: {json.dumps(warnings_list)}\n"
    )
    (run_dir / "build_log.txt").write_text(build_log, encoding="utf-8")

    return metrics


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Paleo BUILD_AND_SIM artifact producer (v2.0)")
    parser.add_argument("--taxon", required=True, help="Slug under models/")
    parser.add_argument("--run-id", default=None, help="Run ID (auto-generated if omitted)")
    parser.add_argument("--synth", choices=["kl", "legacy_iir"], default="kl")
    parser.add_argument("--f0", type=float, default=None, help="f0 soft-tissue prior in Hz")
    parser.add_argument("--temperature", type=float, default=None, help="Air temperature in C (default 20)")
    parser.add_argument("--radiation", choices=["unflanged", "flanged"], default="unflanged")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--duration", type=float, default=2.5)
    parser.add_argument("--sr", type=int, default=44100)
    parser.add_argument("--fmax", type=float, default=4000.0)
    parser.add_argument("--airway-profile", type=Path, default=None, help="Validated airway CSV with .source.json sidecar")
    parser.add_argument("--airway-stl", type=Path, default=None, help="Lumen STL mesh (requires --inlet and --outlet)")
    parser.add_argument("--inlet", type=float, nargs=3, metavar=("X", "Y", "Z"), default=None)
    parser.add_argument("--outlet", type=float, nargs=3, metavar=("X", "Y", "Z"), default=None)
    args = parser.parse_args(argv)

    cli_cmd = "python -m pipeline.build_and_sim " + " ".join(argv if argv is not None else sys.argv[1:])
    try:
        metrics = run_build_and_sim(
            args.taxon,
            args.run_id,
            synth_mode=args.synth,
            f0_hz=args.f0,
            temperature_c=args.temperature,
            radiation=args.radiation,
            seed=args.seed,
            duration_s=args.duration,
            sample_rate=args.sr,
            fmax=args.fmax,
            airway_profile=args.airway_profile,
            airway_stl=args.airway_stl,
            inlet=tuple(args.inlet) if args.inlet is not None else None,
            outlet=tuple(args.outlet) if args.outlet is not None else None,
            cli_command=cli_cmd,
        )
    except GeometryInputError as exc:
        print(f"ERROR (GeometryInputError): {exc}", file=sys.stderr)
        raise SystemExit(2) from exc

    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()

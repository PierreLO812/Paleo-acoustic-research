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

from pipeline.epb import WITMER_1995_REF, apply_epb_table_to_config
from pipeline.geometry_qc import (
    GeometryInputError,
    apply_end_section_trim,
    evaluate_geometry_qc,
    extract_airway_from_stl,
    load_looped_centerline,
    plot_geometry_qc,
)
from pipeline.hypothesis import load_metadata, write_hypothesis_bundle
from pipeline.measure import compute_shaping_metrics, measure_run_dir, measure_transfer_function
from pipeline.paths import MODELS, OUTPUT, ROOT
from pipeline.peer_review import review_run
from pipeline.sensitivity import run_sensitivity
from pipeline.soft_tissue import (
    SoftTissueError,
    apply_lumen_narrowing,
    load_hypotheses,
    rigid_soft_tissue_config,
)
from pipeline.source_tract import (
    build_kl_waveguide,
    compare_kl_vs_tmm,
    kl_impulse_response,
    synthesize_kl,
    synthesize_legacy_iir,
    synthesize_tmm_fir,
)
from pipeline.tmm_lossy import Air, WallModel, plot_tmm_lossy, run_tmm_profile, solve_tmm

logger = logging.getLogger(__name__)


def _new_run_id(taxon_slug: str) -> str:
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    return f"{taxon_slug}_{stamp}_{uuid.uuid4().hex[:6]}"


_CACHED_GIT_BIN: str | None = None


def _run_git(args: list[str], root: Path = ROOT) -> str | None:
    global _CACHED_GIT_BIN
    candidates = [_CACHED_GIT_BIN] if _CACHED_GIT_BIN else [r"C:\Program Files\Git\cmd\git.exe", "git"]
    for git_bin in candidates:
        try:
            proc = subprocess.run(
                [git_bin, *args],
                cwd=root,
                capture_output=True,
                text=True,
                timeout=5,
                check=True,
            )
            _CACHED_GIT_BIN = git_bin
            return proc.stdout
        except FileNotFoundError:
            continue
        except Exception:
            return None
    return None


def _git_provenance(root: Path = ROOT) -> dict:
    """Capture runtime git commit, full hash, dirty status (ignoring output/), and describe."""
    short = (_run_git(["rev-parse", "--short", "HEAD"], root) or "").strip() or "unknown"
    full = (_run_git(["rev-parse", "HEAD"], root) or "").strip() or "unknown"
    status_out = _run_git(["status", "--porcelain"], root)
    dirty_paths: list[str] = []
    if status_out is not None:
        for raw_line in status_out.splitlines():
            if len(raw_line) < 4:
                continue
            path_part = raw_line[3:].strip().strip('"').replace("\\", "/")
            if " -> " in path_part:
                path_part = path_part.split(" -> ")[-1].strip().strip('"')
            if path_part.startswith("output/") or path_part == "output":
                continue
            dirty_paths.append(path_part)
        git_dirty = len(dirty_paths) > 0
    else:
        git_dirty = False
    desc = (_run_git(["describe", "--tags", "--always", "--dirty"], root) or "").strip()
    if not desc:
        desc = f"{short}-dirty" if git_dirty else short
    return {
        "git_commit": short,
        "git_commit_full": full,
        "git_dirty": git_dirty,
        "git_describe": desc,
        "dirty_paths": dirty_paths,
    }


def _git_commit() -> str:
    return _git_provenance(ROOT)["git_commit"]


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
    fig.savefig(path_png, dpi=96)
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

    prov = metrics.get("provenance") or {}
    dirty_badge = " ⚠ dirty tree" if prov.get("git_dirty") else ""
    shaping = metrics.get("shaping") or {}
    sg_b3db = shaping.get("shaping_gain_db_F1_b3db", metrics["shaping_gain_db_F1"])
    sg_pm60 = shaping.get("shaping_gain_db_F1_pm60", metrics["shaping_gain_db_F1"])
    e_pm60 = shaping.get("f1_band_energy_fraction_measured_pm60")
    kl_block = metrics.get("kl_vs_tmm") or {}

    lines: list[str] = [
        f"# Paleo-Acoustic Run Report — `{metrics['run_id']}`",
        "",
        "## 1. Summary",
        "",
        "| Field | Value |",
        "| :--- | :--- |",
        f"| Taxon / Specimen | {metrics['taxon']} ({metrics['specimen_id']}) |",
        f"| Synth Mode | `{metrics['synth_mode']}` (KL loss: `{kl_block.get('loss_model', 'legacy')}`) |",
        f"| Provenance | `{prov.get('git_commit', metrics['git_commit'])}` (`{prov.get('git_describe', 'unknown')}`){dirty_badge} |",
        f"| Terminal State | **`{state}`** |",
        f"| Why | {why} |",
        f"| Measured f0 (YIN) | {metrics['f0']['measured_hz']:.2f} ± {metrics['f0']['measured_std_hz']:.2f} Hz (prior: {metrics['f0']['prior_hz']} Hz) |",
        f"| Shaping Gain (`F1 ± B1/2`) | {sg_b3db:.2f} dB (legacy `±60 Hz`: {sg_pm60:.2f} dB; floor `±60 Hz`: {e_pm60*100:.2f}% if available) |".replace(
            " if available", ""
        )
        if e_pm60 is not None
        else f"| Shaping Gain (`F1 ± B1/2`) | {sg_b3db:.2f} dB (legacy `±60 Hz`: {sg_pm60:.2f} dB) |",
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
            "![Contrôle qualité géométrique du conduit](geometry_qc.png)",
            "",
            "*Figure 1 — Profil géométrique issu du maillage osseux (section transversale $A(s)$ et diamètre hydraulique $D_h(s)$) : les valeurs dépassent les seuils d'un lumen nasal et classent ce profil en mode exploratoire.*",
            "",
        ]
    )

    end_sec = qc.get("end_section") or {}
    if end_sec.get("end_section_artifact") or end_sec.get("trimmed"):
        lines.extend(
            [
                "### End-Section Artifact QC (With vs Without Trim)",
                "",
                f"- **Flagged end indices:** {end_sec.get('flagged_indices', [])} | **Trim mode:** `{end_sec.get('trim_mode', '0')}` | **Trimmed:** `{end_sec.get('trimmed', False)}`",
                f"- **Length change:** {end_sec.get('original_L_mm')} mm -> {end_sec.get('trimmed_L_mm')} mm",
            ]
        )
        cmp_trim = end_sec.get("formants_untrimmed_vs_trimmed")
        if cmp_trim:
            lines.extend(
                [
                    "",
                    "| Mode | Untrimmed F (Hz) | Untrimmed Q | Trimmed F (Hz) | Trimmed Q |",
                    "| :--- | :--- | :--- | :--- | :--- |",
                ]
            )
            for row in cmp_trim:
                lines.append(
                    f"| {row['mode']} | {row['untrimmed_f_hz']:.2f} | {row['untrimmed_q']:.1f} | {row['trimmed_f_hz']:.2f} | {row['trimmed_q']:.1f} |"
                )
        lines.append("")

    lines.extend(
        [
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
            "![Fonction de transfert acoustique mesurée H(f)](transfer_function.png)",
            "",
            "*Figure 2 — Fonction de transfert acoustique mesurée $H(f) = P_{sf}/P_{ss}$ et cohérence $\\gamma^2(f)$ calculées uniquement à partir des signaux audio `source.wav` et `filtered.wav`.*",
            "",
            "## 5. Lossy Kirchhoff TMM Spectrum",
            "",
            "![Spectre TMM Kirchhoff avec pertes](tmm_lossy.png)",
            "",
            "*Figure 3 — Réponse en fréquence prédite par le solveur indépendant TMM (pertes thermovisqueuses de Kirchhoff, impédance de rayonnement et parois compliantes si activées).*",
            "",
            "## 6. Kelly–Lochbaum Waveguide vs Lossy TMM",
            "",
        ]
    )

    if kl_cmp and kl_cmp.get("comparisons"):
        lines.extend(
            [
                f"- **Loss model:** `{kl_cmp.get('loss_model', 'legacy')}`"
                + (
                    f" (fit RMS dB: median={kl_cmp['fit_rms_db']['median']}, max={kl_cmp['fit_rms_db']['max']})"
                    if kl_cmp.get("fit_rms_db")
                    else ""
                ),
                "",
                "| Formant | TMM (Hz) | KL Impulse (Hz) | F Rel Error | TMM Q | KL Q | Q Rel Error |",
                "| :--- | :--- | :--- | :--- | :--- | :--- | :--- |",
            ]
        )
        for item in kl_cmp["comparisons"]:
            lines.append(
                f"| {item.get('formant', '—')} | {item['tmm_hz']:.2f} | {item['kl_hz']:.2f} | "
                f"{item['rel_err']*100:.3f}% | {item['tmm_q']:.1f} | {item['kl_q']:.1f} | "
                f"{item.get('q_rel_err', 0.0)*100:.2f}% |"
            )
        lines.append("")
    lines.extend(
        [
            "![Comparaison guide d'onde KL vs TMM](kl_vs_tmm.png)",
            "",
            "*Figure 4 — Superposition de la réponse impulsionnelle du guide d'onde numérique Kelly–Lochbaum et de la courbe TMM avec pertes : les filtres biquad par section alignent à la fois les fréquences formantiques et les facteurs de qualité $Q$.*",
            "",
            "## 7. Band Energy Fractions & Spectral Shaping",
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

    if shaping:
        b3 = shaping.get("b3db") or {}
        p60 = shaping.get("pm60") or {}
        lines.extend(
            [
                "",
                f"- **Shaping Gain (`F1 ± B1/2`, {b3.get('band_hz')} Hz):** **{b3.get('shaping_gain_db_F1_b3db')} dB** (used by peer review +3 dB gate; filtered energy fraction: {b3.get('f1_band_energy_fraction_measured_b3db', 0.0)*100:.2f}%)",
                f"- **Shaping Gain (legacy `F1 ± 60 Hz`, {p60.get('band_hz')} Hz):** **{p60.get('shaping_gain_db_F1_pm60')} dB**; energy floor `f1_band_energy_fraction_measured_pm60` = **{p60.get('f1_band_energy_fraction_measured_pm60', 0.0)*100:.3f}%** (threshold >= 2.0%).",
                "- *Why both gates matter:* The legacy `+3 dB` gain check alone did not fail the old `eacbdf` run (which scored `+24.59 dB` because its source had almost zero energy near F1); it is the **>= 2 % F1-band energy floor** (`f1_band_energy_fraction_measured_pm60`) that catches an inaudible filter.",
            ]
        )

    lines.extend(
        [
            "",
            "## 8. Source Spectral Tilt",
            "",
            f"- Measured harmonic tilt (100–2000 Hz): **{metrics['source_spectral_tilt_db_per_oct']:.2f} dB/octave**",
            "",
            "## 9. Hypothèses tissus mous",
            "",
            "**Ces valeurs sont des hypothèses ; elles ne changent pas l'état EXPLORATORY.**",
            "",
        ]
    )

    st_block = (metrics.get("assumptions") or {}).get("soft_tissue") or {}
    st_params = st_block.get("parameters") or {}
    kind_labels = {
        "hypothesis": "hypothèse",
        "cited": "citée",
        "analog_human": "analogue humain",
        "analog_extant": "analogue actuel",
        "epb_bracket": "citée (EPB)",
    }
    if st_params:
        lines.extend(
            [
                "| Paramètre | Valeur nominale | Fourchette | Unité | Type (citée / hypothèse / analogue) | Source | Confiance |",
                "| :--- | :--- | :--- | :--- | :--- | :--- | :--- |",
            ]
        )
        for pname, pinfo in st_params.items():
            rng = pinfo.get("range", [None, None])
            rng_s = f"[{rng[0]:.4g}, {rng[1]:.4g}]"
            src_disp = pinfo.get("source", "hypothesis")
            if pinfo.get("model_source"):
                src_disp = f"{src_disp} (modèle : {pinfo['model_source']})"
            lines.append(
                f"| `{pname}` | {pinfo.get('effective_value')} | {rng_s} | `{pinfo.get('unit')}` | "
                f"{kind_labels.get(pinfo.get('kind', ''), pinfo.get('kind', ''))} | {src_disp} | {pinfo.get('confidence')} |"
            )
        lines.append("")

    # Bone only vs soft tissue (nominal)
    bone_pks = metrics.get("formants_predicted_tmm_bone") or []
    soft_pks = metrics.get("formants_predicted_tmm_soft")
    lines.extend(["## 10. Bone only vs soft tissue (nominal)", ""])
    if bone_pks and soft_pks:
        lines.extend(
            [
                "| Mode | Os seul F (Hz) | Tissus mous F (Hz) | ΔF (Hz / %) | Os seul Q | Tissus mous Q | ΔQ | Os seul B (Hz) | Tissus mous B (Hz) | ΔB (Hz) |",
                "| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |",
            ]
        )
        for i in range(min(3, len(bone_pks), len(soft_pks))):
            bp, sp = bone_pks[i], soft_pks[i]
            df = sp["f_hz"] - bp["f_hz"]
            df_pct = 100.0 * df / max(bp["f_hz"], 1e-9)
            dq = sp["q"] - bp["q"]
            db = sp["bw_3db_hz"] - bp["bw_3db_hz"]
            lines.append(
                f"| F{i+1} | {bp['f_hz']:.2f} | {sp['f_hz']:.2f} | {df:+.2f} ({df_pct:+.1f}%) | "
                f"{bp['q']:.1f} | {sp['q']:.1f} | {dq:+.1f} | {bp['bw_3db_hz']:.2f} | {sp['bw_3db_hz']:.2f} | {db:+.2f} |"
            )
        lines.append("")
        lines.append(
            "*Note physique : le retrait de muqueuse est un exercice de sensibilité sur géométrie osseuse (pas une reconstruction de lumen) ; les parois massiques relèvent $F_1$ et élargissent la bande passante $B_1$.*"
        )
    else:
        lines.append("- Soft-tissue modeling disabled (`--soft-tissue-off`); showing rigid bone baseline only.")
    lines.append("")

    # Sensitivity section
    sens = metrics.get("sensitivity") or {}
    sens_sum = sens.get("summary")
    if sens_sum:
        rel_sens = sens.get("rel_dir", "sensitivity")
        lines.extend(
            [
                "## 11. Sensitivity & Uncertainty Analysis",
                "",
                f"- **Samples:** n={sens_sum.get('n_samples')} (valid={sens_sum.get('n_valid')}, invalid/occluded={sens_sum.get('n_invalid')}), seed={sens_sum.get('seed')}",
                "",
                f"![Diagramme tornade des sensibilités univariées]({rel_sens}/tornado.png)",
                "",
                "*Figure 5 — Diagramme tornade (OAT) montrant l'impact individuel de chaque hypothèse de tissu mou sur $F_1$, $Q_1$ et le gain de filtrage prédit.*",
                "",
                f"![Enveloppe spectrale de sensibilité]({rel_sens}/spaghetti_H.png)",
                "",
                "*Figure 6 — Faisceau de réponses $|H_{TMM}(f)|$ issues de l'échantillonnage LHS (gris), intervalle 5–95 %, scénario nominal (noir) et référence os seul de l'étape 2 (tireté bleu).*",
                "",
                f"![Distributions des formants et facteurs Q]({rel_sens}/formant_distributions.png)",
                "",
                "*Figure 7 — Distributions empiriques de $F_1\\text{–}F_3$ et $Q_1\\text{–}Q_3$ sur l'ensemble des tirages non occlus, comparées aux repères os seul et nominal.*",
                "",
                "### 5–95 % Intervals",
                "",
                "| Output | Median | 5 % | 95 % |",
                "| :--- | :--- | :--- | :--- |",
            ]
        )
        ivals = sens_sum.get("intervals_p05_median_p95") or {}
        for k_out in ("F1_hz", "F2_hz", "F3_hz", "Q1", "Q2", "Q3", "B1_hz", "shaping_gain_db_F1_b3db_predicted"):
            iv = ivals.get(k_out) or {}
            lines.append(f"| `{k_out}` | {iv.get('median')} | {iv.get('p05')} | {iv.get('p95')} |")
        lines.extend(["", "### Parameter Influence Ranking (Spearman / PRCC)", ""])
        ranks = (sens_sum.get("importance") or {}).get("rankings") or {}
        for tgt in ("F1_hz", "Q1", "shaping_gain_db_F1_b3db_predicted"):
            top = ranks.get(tgt) or []
            top_str = ", ".join(f"`{r['parameter']}` (ρ={r['spearman_rho']:+.2f})" for r in top[:4])
            lines.append(f"- **{tgt}:** {top_str}")
        scens = sens_sum.get("scenarios") or []
        if scens:
            lines.extend(["", "### Sensitivity Audio Scenarios (min / median / max damping)", ""])
            for sc in scens:
                p = sc.get("parameters") or {}
                pk1 = (sc.get("predicted_tmm_peaks") or [{}])[0]
                lines.append(
                    f"- [`{rel_sens}/{sc['wav_file']}`]({rel_sens}/{sc['wav_file']}) — Scénario **{sc['scenario']}** : "
                    f"$t={p.get('mucosa_thickness_mm', 0):.2f}$ mm, $R_w={p.get('wall_resistance_per_area', 0):.0f}$ Pa·s/m, "
                    f"$f_0={p.get('f0_prior_hz', 60):.1f}$ Hz → $F_1={pk1.get('f_hz', '—')}$ Hz, $B_1={pk1.get('bw_3db_hz', '—')}$ Hz "
                    f"(pic {sc.get('audio_peak_dbfs')} dBFS, clipping={sc.get('clipping_detected')})."
                )
        lines.append("")

    # What-if synthetic section
    what_if = metrics.get("what_if") or {}
    if what_if.get("enabled"):
        lines.extend(
            [
                "## 12. Synthetic Length-Scale What-If Scenario",
                "",
                "> [!CAUTION]",
                "> **SCÉNARIO SYNTHÉTIQUE ILLUSTRATIF (`synthetic: true`) — illustrative scenario, not a loop, not validated.** Stretching the bone profile along $s$ does not represent the anatomical nasal loop and produces no terminal validation state.",
                "",
                f"- **Scale factor:** {what_if.get('scale')} | **Outputs:** {what_if.get('paths')}",
                "",
            ]
        )

    # EPB & Looped airway scaffolding status
    epb_blk = metrics.get("epb") or {}
    loop_blk = metrics.get("looped_airway") or {}
    lines.extend(
        [
            "## 13. EPB & Looped Airway Status",
            "",
            f"- **EPB comparative data:** `{epb_blk.get('status')}` (concept ref: {WITMER_1995_REF})",
            f"- **Looped nasal airway:** `{loop_blk.get('status')}` (`loop_to_crest_ratio={loop_blk.get('loop_to_crest_ratio')}`)",
            "",
            "## 14. Citations to verify",
            "",
        ]
    )
    cites = list(st_block.get("citations_to_verify") or [])
    cites.append(
        {
            "parameter": "radiation_impedance",
            "role": "model_equation",
            "source": "Silva, F., Guillemain, P., Kergomard, J., Debut, V., & Kemp, A. (2009). Approximation formulae for the acoustic radiation impedance of a cylindrical pipe. J. Sound Vib., 322(1-2), 255-263.",
            "source_locator": "Eq. (unflanged rational fit)",
            "verified_by_human": False,
        }
    )
    cites.append(
        {
            "parameter": "epb_framework",
            "role": "methodology",
            "source": WITMER_1995_REF,
            "source_locator": "Chapter 2",
            "verified_by_human": False,
        }
    )
    for c_item in cites:
        lines.append(
            f"- `{c_item.get('parameter')}` ({c_item.get('role')}): {c_item.get('source')} "
            f"[{c_item.get('source_locator')}] — `verified_by_human: {c_item.get('verified_by_human', False)}`"
        )

    lines.extend(
        [
            "",
            "## 15. Assumptions (Full JSON)",
            "",
            f"```json\n{json.dumps(metrics['assumptions'], indent=2)}\n```",
            "",
            "## 16. Warnings",
            "",
            "\n".join(f"- {w}" for w in metrics["warnings"]) if metrics["warnings"] else "- None",
            "",
            "## 17. Limits",
            "",
            "- **Bone-derived profile:** Unless a validated CT lumen or endocast profile with `.source.json` is supplied, the cross-sectional profile reflects bone wall geometry (`RightSkeleton_RAM14000_printable2.stl`, L = 263.8 mm single pass) rather than the true looped nasal passage.",
            "- **Uniform mucosa sensitivity exercise:** Subtracting a uniform inward offset `t` from bone cross-sections (`A_lumen = A_bone - P·t + π·t²`) is strictly a sensitivity exercise on a bone-derived profile, not a lumen reconstruction.",
            "- **Human-analog wall equations:** Yielding walls use the lumped RLC wall impedance formulation from human vocal-tract acoustics (Ishizaka et al. 1975; Sondhi 1974) with wide unverified `hypothesis` parameter ranges.",
            "- **Unknown body temperature:** Airway temperature is unobservable in fossils; `{20, 37}` °C compares ambient vs endotherm-like scenarios.",
            f"- **f0 prior:** `f0_prior_hz = {metrics['f0']['prior_hz']} Hz` is a soft-tissue assumption, not an anatomical measurement.",
            "- **Loop not represented:** The ascending/descending nasal loop through the crest is not geometrically segmented (`looped airway length unverified`).",
            "- **Internal-consistency-only agreement:** Measured `H(f)`, TMM, and KL all derive from the same bone-derived `A(s)` and cannot externally validate fossil vocalizations.",
            "",
            "## 18. Artifacts",
            "",
        ]
    )
    for k, v in metrics.get("artifacts", {}).items():
        lines.append(f"- `{k}`: `{v}`")

    return "\n".join(lines) + "\n"


def _run_what_if_length_scale(
    prof_active: dict,
    air: Air,
    wall: WallModel | None,
    extra_damping: float,
    scale: float,
    *,
    radiation: str,
    f0_prior: float,
    seed: int,
    duration_s: float,
    sample_rate: int,
    fmax: float,
    kl_loss_model: str,
    what_if_dir: Path,
) -> dict:
    """Stretch A(s) by scale into output/<RUN_ID>/what_if/SYNTHETIC_L_x* with no state file."""
    if scale <= 0.0:
        raise ValueError(f"--what-if-length-scale must be > 0 (got {scale}).")
    what_if_dir.mkdir(parents=True, exist_ok=True)
    tag = str(int(round(scale))) if abs(scale - round(scale)) < 1e-6 else f"{scale:.2f}".replace(".", "p")
    prefix = f"SYNTHETIC_L_x{tag}"

    s0 = float(prof_active["s_mm"][0])
    s_scaled = s0 + float(scale) * (np.asarray(prof_active["s_mm"], dtype=np.float64) - s0)
    area_scaled = np.asarray(prof_active["area_mm2"], dtype=np.float64).copy()
    dh_scaled = np.asarray(prof_active["dh_mm"], dtype=np.float64).copy()
    perim_scaled = (
        np.asarray(prof_active["perimeter_mm"], dtype=np.float64).copy()
        if prof_active.get("perimeter_mm") is not None
        else None
    )

    tmm_wi = solve_tmm(
        s_scaled,
        area_scaled,
        dh_scaled,
        air=air,
        fmax=fmax,
        df=0.5,
        radiation=radiation,
        losses=True,
        wall=wall,
        extra_damping_factor=extra_damping,
    )
    tmm_png = what_if_dir / f"{prefix}_tmm_lossy.png"
    plot_tmm_lossy(tmm_wi, tmm_png)

    synth_wi = synthesize_kl(
        s_scaled,
        area_scaled,
        air,
        {"peaks": tmm_wi["h_peaks"]},
        dh_mm=dh_scaled,
        perimeter_mm=perim_scaled,
        duration_s=duration_s,
        sr=sample_rate,
        f0_hz=f0_prior,
        radiation=radiation,
        seed=seed,
        kl_loss_model=kl_loss_model,
        extra_damping_factor=extra_damping,
        wall=wall,
    )
    src_wav = what_if_dir / f"{prefix}_source.wav"
    flt_wav = what_if_dir / f"{prefix}_filtered.wav"
    _write_wav(src_wav, synth_wi["source"], sample_rate)
    _write_wav(flt_wav, synth_wi["filtered"], sample_rate)

    tf_wi = measure_transfer_function(synth_wi["source"], synth_wi["filtered"], sample_rate)
    from pipeline.measure import plot_transfer_function

    tf_png = what_if_dir / f"{prefix}_transfer_function.png"
    plot_transfer_function(tf_wi, tf_png)

    metrics_wi = {
        "synthetic": True,
        "banner": "illustrative scenario, not a loop, not validated",
        "scale": float(scale),
        "L_mm": round(float(s_scaled[-1] - s_scaled[0]), 3),
        "formants_predicted_tmm": [
            {"f_hz": p["f_hz"], "bw_3db_hz": p["bw_3db_hz"], "q": p["q"]}
            for p in tmm_wi["h_peaks"][:4]
        ],
        "formants_measured_tf": tf_wi["peaks"][:4],
        "audio": synth_wi["audio"],
        "artifacts": {
            "source_wav": str(src_wav),
            "filtered_wav": str(flt_wav),
            "tmm_lossy_png": str(tmm_png),
            "transfer_function_png": str(tf_png),
        },
    }
    metrics_wi_path = what_if_dir / f"{prefix}_metrics.json"
    metrics_wi_path.write_text(json.dumps(metrics_wi, indent=2) + "\n", encoding="utf-8")

    return {
        "enabled": True,
        "synthetic": True,
        "scale": float(scale),
        "banner": "illustrative scenario, not a loop, not validated",
        "paths": [
            str(src_wav),
            str(flt_wav),
            str(tmm_png),
            str(tf_png),
            str(metrics_wi_path),
        ],
        "formants_predicted_tmm": metrics_wi["formants_predicted_tmm"],
    }


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
    kl_loss_model: str = "kirchhoff_fit",
    trim_end_sections: str = "0",
    soft_tissue_path: Path | None = None,
    soft_tissue_off: bool = False,
    sensitivity_n: int | None = None,
    epb_table: Path | None = None,
    airway_centerline: Path | None = None,
    what_if_length_scale: float | None = None,
) -> dict:
    t_start = time.perf_counter()
    meta = load_metadata(taxon_slug)
    run_id = run_id or _new_run_id(taxon_slug)
    run_dir = OUTPUT / run_id
    run_dir.mkdir(parents=True, exist_ok=True)

    warnings_list: list[str] = []
    pkg_versions = _package_versions(warnings_list)
    prov = _git_provenance(ROOT)
    if prov["git_dirty"]:
        dirty_msg = f"Git working tree is dirty (uncommitted changes: {prov['dirty_paths']})."
        logger.warning(dirty_msg)
        warnings_list.append(dirty_msg)

    # Write HYPOTHESIS bundle & state_HYPOTHESIS.json if not already emitted
    if not (run_dir / "antigravity_prompt.md").exists():
        hyp_state = write_hypothesis_bundle(taxon_slug, run_id)
    else:
        hyp_state = {"state": "HYPOTHESIS", "run_id": run_id, "taxon": taxon_slug}
    (run_dir / "state_HYPOTHESIS.json").write_text(json.dumps(hyp_state, indent=2) + "\n", encoding="utf-8")

    # Resolve geometry input (M5 + M4 centerline scaffolding)
    crest_len_mm = float(meta["geometry"].get("conduit_length_mm", 263.804))
    if airway_centerline is not None:
        profile_path = Path(airway_centerline)
        if not profile_path.exists() and (ROOT / profile_path).exists():
            profile_path = ROOT / profile_path
        qc, prof = load_looped_centerline(profile_path, crest_length_mm=crest_len_mm)
        looped_info = qc["looped_airway"]
    elif airway_stl is not None:
        extracted_csv = run_dir / "extracted_airway.csv"
        profile_path = extract_airway_from_stl(Path(airway_stl), inlet, outlet, extracted_csv)
        qc, prof = evaluate_geometry_qc(
            profile_path,
            source_mesh=str(airway_stl),
            extraction_method="lumen_stl_centerline",
            require_sidecar=False,
        )
        looped_info = {"status": "stl_centerline_unverified_loop", "loop_to_crest_ratio": None}
    elif airway_profile is not None:
        cand = Path(airway_profile)
        profile_path = cand if cand.exists() else (ROOT / cand)
        qc, prof = evaluate_geometry_qc(profile_path, require_sidecar=True)
        looped_info = {"status": "validated_airway_profile", "loop_to_crest_ratio": None}
    else:
        rel_csv = meta["geometry"]["airway_features_csv"]
        profile_path = ROOT / rel_csv
        qc, prof = evaluate_geometry_qc(
            profile_path,
            source_mesh=meta["geometry"]["stl_path"],
            extraction_method="bone_mesh_skeleton_axis_slice",
            require_sidecar=False,
        )
        looped_info = {
            "status": "unverified",
            "flag": "looped airway length unverified" if "looped airway length unverified" in qc["flags"] else None,
            "loop_to_crest_ratio": None,
        }

    qc_untrimmed, prof_untrimmed = qc, prof
    qc, prof, trim_msgs = apply_end_section_trim(qc, prof, trim_end_sections)
    warnings_list.extend(trim_msgs)

    geom_png = run_dir / "geometry_qc.png"
    plot_geometry_qc(qc, prof, geom_png)

    # Resolve soft-tissue hypothesis layer (M1) and EPB scaffolding (M3)
    priors = meta["acoustic_priors"]
    soft_meta = meta["soft_tissue_assumptions"]
    if soft_tissue_off:
        st_cfg = rigid_soft_tissue_config()
        soft_enabled = False
    else:
        default_hyp_path = MODELS / taxon_slug / "soft_tissue_hypotheses.json"
        if soft_tissue_path is not None:
            st_cfg = load_hypotheses(Path(soft_tissue_path))
            soft_enabled = True
        elif default_hyp_path.exists():
            st_cfg = load_hypotheses(default_hyp_path)
            soft_enabled = True
        else:
            msg = f"Missing '{default_hyp_path}'; falling back to rigid walls."
            logger.warning(msg)
            warnings_list.append(msg)
            st_cfg = rigid_soft_tissue_config()
            soft_enabled = False

    st_cfg, epb_info, epb_warnings = apply_epb_table_to_config(st_cfg, epb_table)
    for ew in epb_warnings:
        logger.warning(ew)
        warnings_list.append(ew)

    nom = st_cfg.nominal_values()
    if soft_enabled and "lining_thickness_mm" in soft_meta:
        meta_thick = float(soft_meta["lining_thickness_mm"])
        if abs(nom["mucosa_thickness_mm"] - meta_thick) > 1e-6:
            msg = (
                f"metadata.json lining_thickness_mm={meta_thick} mm (superseded_by soft_tissue_hypotheses.json) "
                f"differs from hypothesis nominal mucosa_thickness_mm={nom['mucosa_thickness_mm']:.4g} mm."
            )
            logger.warning(msg)
            warnings_list.append(msg)

    temp_default = float(soft_meta.get("temperature_c", 20.0)) if not soft_enabled else float(nom["air_temperature_c"])
    temp_used = temp_default if temperature_c is None else float(temperature_c)
    f0_default = float(priors.get("f0_prior_hz", 60.0)) if not soft_enabled else float(nom["f0_prior_hz"])
    f0_prior = f0_default if f0_hz is None else float(f0_hz)
    air = Air.from_temperature(temp_used)

    if soft_enabled:
        prof_active = apply_lumen_narrowing(
            prof,
            nom["mucosa_thickness_mm"],
            min_area_fraction=nom["lumen_min_area_fraction"],
        )
        wall_active = WallModel(
            mass_per_area=nom["wall_mass_per_area"],
            stiffness_per_area=nom["wall_stiffness_per_area"],
            resistance_per_area=nom["wall_resistance_per_area"],
        )
        extra_damping = float(nom["extra_damping_factor"])
        lumen_model = prof_active.get("lumen_model", "uniform_inward_offset")
    else:
        prof_active = prof
        wall_active = None
        extra_damping = 1.0
        lumen_model = "rigid_unmodified"

    # Bone-only rigid TMM baseline (for bone vs soft-tissue table; equals active TMM when soft_tissue_off)
    if soft_enabled:
        bone_tmm = solve_tmm(
            prof["s_mm"],
            prof["area_mm2"],
            prof["dh_mm"],
            air=Air.from_temperature(float(soft_meta.get("temperature_c", 20.0))),
            fmax=fmax,
            df=0.5,
            radiation=radiation,
            losses=True,
        )
        bone_peaks = [
            {"f_hz": p["f_hz"], "bw_3db_hz": p["bw_3db_hz"], "q": p["q"]}
            for p in bone_tmm["h_peaks"][:4]
        ]
    else:
        bone_peaks = None

    # M2: Lossy TMM on active geometry & wall configuration
    tmm_out = run_tmm_profile(
        profile_path,
        run_dir,
        fmax=fmax,
        df=0.5,
        temperature_c=temp_used,
        radiation=radiation,
        write_outputs=True,
        prof=prof_active,
        wall=wall_active,
        extra_damping_factor=extra_damping,
    )
    tmm_summary = tmm_out["summary"]
    tmm_raw = tmm_out["raw"]
    warnings_list.extend(tmm_summary.get("warnings", []))

    if bone_peaks is None:
        bone_peaks = [
            {"f_hz": p["f_hz"], "bw_3db_hz": p["bw_3db_hz"], "q": p["q"]}
            for p in tmm_summary["peaks"][:4]
        ]
    soft_peaks = (
        [
            {"f_hz": p["f_hz"], "bw_3db_hz": p["bw_3db_hz"], "q": p["q"]}
            for p in tmm_summary["peaks"][:4]
        ]
        if soft_enabled
        else None
    )

    # Compare F1-F3 with vs without end-section trim whenever end_section_artifact is present or trim applied
    end_sec = qc.get("end_section") or {}
    if end_sec.get("end_section_artifact") or end_sec.get("trimmed"):
        if end_sec.get("trimmed"):
            untrim_tmm = solve_tmm(
                prof_untrimmed["s_mm"],
                prof_untrimmed["area_mm2"],
                prof_untrimmed["dh_mm"],
                air=air,
                fmax=fmax,
                df=0.5,
                radiation=radiation,
            )
            p_untrim = untrim_tmm["h_peaks"]
            p_trim_cmp = solve_tmm(
                prof["s_mm"],
                prof["area_mm2"],
                prof["dh_mm"],
                air=air,
                fmax=fmax,
                df=0.5,
                radiation=radiation,
            )["h_peaks"]
        else:
            _, prof_auto, _ = apply_end_section_trim(qc, prof, "auto")
            trim_tmm = solve_tmm(
                prof_auto["s_mm"],
                prof_auto["area_mm2"],
                prof_auto["dh_mm"],
                air=air,
                fmax=fmax,
                df=0.5,
                radiation=radiation,
            )
            p_untrim = solve_tmm(
                prof["s_mm"],
                prof["area_mm2"],
                prof["dh_mm"],
                air=air,
                fmax=fmax,
                df=0.5,
                radiation=radiation,
            )["h_peaks"]
            p_trim_cmp = trim_tmm["h_peaks"]
        qc["end_section"]["formants_untrimmed_vs_trimmed"] = [
            {
                "mode": f"F{i + 1}",
                "untrimmed_f_hz": p_untrim[i]["f_hz"],
                "untrimmed_q": p_untrim[i]["q"],
                "trimmed_f_hz": p_trim_cmp[i]["f_hz"],
                "trimmed_q": p_trim_cmp[i]["q"],
            }
            for i in range(min(3, len(p_untrim), len(p_trim_cmp)))
        ]

    # M4: Synthesis (--synth kl, legacy_iir, or tmm_fir)
    kl_png = run_dir / "kl_vs_tmm.png"
    if synth_mode == "kl":
        synth = synthesize_kl(
            prof_active["s_mm"],
            prof_active["area_mm2"],
            air,
            tmm_summary,
            dh_mm=prof_active["dh_mm"],
            perimeter_mm=prof_active.get("perimeter_mm"),
            duration_s=duration_s,
            sr=sample_rate,
            f0_hz=f0_prior,
            radiation=radiation,
            seed=seed,
            kl_plot_path=kl_png,
            tmm_raw=tmm_raw,
            kl_loss_model=kl_loss_model,
            extra_damping_factor=extra_damping,
            wall=wall_active,
        )
        kl_cmp = synth["kl_check"]
        independent_of_tmm = True
    elif synth_mode == "legacy_iir":
        synth = synthesize_legacy_iir(
            list(priors["expected_formants_hz"]),
            list(priors["default_q_factors"]),
            duration_s=duration_s,
            sr=sample_rate,
            f0_hz=f0_prior,
            seed=seed,
        )
        wg = build_kl_waveguide(
            prof_active["s_mm"],
            prof_active["area_mm2"],
            air,
            dh_mm=prof_active["dh_mm"],
            perimeter_mm=prof_active.get("perimeter_mm"),
            radiation=radiation,
            sr=sample_rate,
            kl_loss_model=kl_loss_model,
            extra_damping_factor=extra_damping,
            wall=wall_active,
            f1_hint_hz=float(tmm_summary["peaks"][0]["f_hz"]) if tmm_summary["peaks"] else 400.0,
        )
        ir_wg = kl_impulse_response(wg, ir_duration_s=0.8)
        kl_cmp = compare_kl_vs_tmm(wg, ir_wg, tmm_summary, out_png=kl_png, tmm_raw=tmm_raw)
        independent_of_tmm = True
    elif synth_mode == "tmm_fir":
        synth = synthesize_tmm_fir(
            tmm_raw,
            air,
            duration_s=duration_s,
            sr=sample_rate,
            f0_hz=f0_prior,
            seed=seed,
        )
        wg = build_kl_waveguide(
            prof_active["s_mm"],
            prof_active["area_mm2"],
            air,
            dh_mm=prof_active["dh_mm"],
            perimeter_mm=prof_active.get("perimeter_mm"),
            radiation=radiation,
            sr=sample_rate,
            kl_loss_model=kl_loss_model,
            extra_damping_factor=extra_damping,
            wall=wall_active,
            f1_hint_hz=float(tmm_summary["peaks"][0]["f_hz"]) if tmm_summary["peaks"] else 400.0,
        )
        ir_wg = kl_impulse_response(wg, ir_duration_s=0.8)
        kl_cmp = compare_kl_vs_tmm(wg, ir_wg, tmm_summary, out_png=kl_png, tmm_raw=tmm_raw)
        independent_of_tmm = False
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

    # Compute both F1 +/- B1/2 and legacy F1 +/- 60 Hz shaping metrics
    f1_tmm_pk = tmm_summary["peaks"][0]
    shaping = compute_shaping_metrics(
        synth["source"],
        synth["filtered"],
        sample_rate,
        float(f1_tmm_pk["f_hz"]),
        float(f1_tmm_pk["bw_3db_hz"]),
    )

    # Optional chained sensitivity sweep (--sensitivity N)
    sensitivity_block: dict = {"status": "not_run", "path": None, "summary": None}
    if sensitivity_n is not None and int(sensitivity_n) > 0:
        sens_dir = run_dir / "sensitivity"
        sens_sum = run_sensitivity(
            taxon_slug,
            n_samples=int(sensitivity_n),
            seed=seed,
            out_dir=sens_dir,
            soft_tissue_path=soft_tissue_path,
            epb_table=epb_table,
            trim_end_sections=trim_end_sections,
            duration_s=duration_s,
            render_audio=True,
        )
        sensitivity_block = {
            "status": "completed",
            "path": str(sens_dir),
            "rel_dir": "sensitivity",
            "summary": {
                "n_samples": sens_sum["n_samples"],
                "n_valid": sens_sum["n_valid"],
                "n_invalid": sens_sum["n_invalid"],
                "seed": sens_sum["seed"],
                "intervals_p05_median_p95": sens_sum["intervals_p05_median_p95"],
                "importance": sens_sum["importance"],
                "scenarios": sens_sum["scenarios"],
            },
        }

    # Optional synthetic length-scale what-if (--what-if-length-scale)
    what_if_block: dict = {"enabled": False, "synthetic": False, "paths": []}
    if what_if_length_scale is not None:
        what_if_block = _run_what_if_length_scale(
            prof_active,
            air,
            wall_active,
            extra_damping,
            float(what_if_length_scale),
            radiation=radiation,
            f0_prior=f0_prior,
            seed=seed,
            duration_s=duration_s,
            sample_rate=sample_rate,
            fmax=fmax,
            kl_loss_model=kl_loss_model,
            what_if_dir=run_dir / "what_if",
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

    soft_metrics_block = st_cfg.to_metrics_dict(
        lumen_model=lumen_model,
        wall_model_active=wall_active is not None,
        soft_tissue_enabled=soft_enabled,
    )
    soft_metrics_block.update(
        {
            "f0_prior_hz": f0_prior,
            "source_model": synth["source_description"],
            "tract_model": synth["filtered_description"],
            "legacy_metadata_lining_thickness_mm": soft_meta.get("lining_thickness_mm"),
        }
    )

    metrics = {
        "schema_version": "2.1",
        "taxon": meta["taxon"],
        "specimen_id": meta["specimen_id"],
        "taxon_slug": taxon_slug,
        "run_id": run_id,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "git_commit": prov["git_commit"],
        "provenance": {
            "git_commit": prov["git_commit"],
            "git_commit_full": prov["git_commit_full"],
            "git_dirty": prov["git_dirty"],
            "git_describe": prov["git_describe"],
        },
        "seed": int(seed),
        "sample_rate": int(sample_rate),
        "duration_s": float(duration_s),
        "synth_mode": synth_mode,
        "independent_of_tmm": bool(independent_of_tmm),
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
        "formants_predicted_tmm_bone": bone_peaks,
        "formants_predicted_tmm_soft": soft_peaks,
        "analytic": {
            "L_m": tmm_summary["L_m"],
            "L_eff_m": tmm_summary["L_eff_m"],
            "c_m_s": air.c,
            "closed_open_hz": tmm_summary["analytic_closed_open_hz"],
            "open_open_hz": tmm_summary["analytic_open_open_hz"],
        },
        "band_energy_fractions": measured["band_energy_fractions"],
        "shaping": shaping,
        "shaping_gain_db_F1": shaping["shaping_gain_db_F1_pm60"],
        "source_spectral_tilt_db_per_oct": measured["source_spectral_tilt_db_per_oct"],
        "kl_vs_tmm_rel_err": kl_cmp["rel_errors"] if kl_cmp else [],
        "kl_vs_tmm": {
            "freq_rel_err": kl_cmp["rel_errors"] if kl_cmp else [],
            "q_rel_err": kl_cmp["q_rel_errors"] if kl_cmp else [],
            "loss_model": kl_cmp["loss_model"] if kl_cmp else kl_loss_model,
            "fit_rms_db": kl_cmp["fit_rms_db"] if kl_cmp else None,
        },
        "geometry_qc": qc,
        "review": {
            "passed": False,
            "issues": [],
            "tolerances": {},
            "retry_feedback": [],
        },
        "assumptions": {
            "air": tmm_summary["air"],
            "soft_tissue": soft_metrics_block,
            "radiation_model": radiation,
        },
        "sensitivity": sensitivity_block,
        "epb": epb_info,
        "looped_airway": looped_info,
        "what_if": what_if_block,
        "audio": {
            "peak_dbfs": round(max(synth["audio"]["source_peak_dbfs"], synth["audio"]["filtered_peak_dbfs"]), 3),
            "source_peak_dbfs": synth["audio"]["source_peak_dbfs"],
            "filtered_peak_dbfs": synth["audio"]["filtered_peak_dbfs"],
            "clipping_detected": bool(synth["audio"]["clipping_detected"] or measured["audio"]["clipping_detected"]),
        },
        "warnings": warnings_list,
        "deprecated": {
            "formant_peaks_hz": "removed: was a copy of filter targets",
            "shaping_gain_db_F1": "deprecated alias of shaping_gain_db_F1_pm60; review gate uses shaping.shaping_gain_db_F1_b3db",
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
    parser = argparse.ArgumentParser(description="Paleo BUILD_AND_SIM artifact producer (v2.1)")
    parser.add_argument("--taxon", required=True, help="Slug under models/")
    parser.add_argument("--run-id", default=None, help="Run ID (auto-generated if omitted)")
    parser.add_argument("--synth", choices=["kl", "legacy_iir", "tmm_fir"], default="kl")
    parser.add_argument("--f0", type=float, default=None, help="f0 soft-tissue prior in Hz")
    parser.add_argument("--temperature", type=float, default=None, help="Air temperature in C")
    parser.add_argument("--radiation", choices=["unflanged", "flanged"], default="unflanged")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--duration", type=float, default=2.5)
    parser.add_argument("--sr", type=int, default=44100)
    parser.add_argument("--fmax", type=float, default=4000.0)
    parser.add_argument("--airway-profile", type=Path, default=None, help="Validated airway CSV with .source.json sidecar")
    parser.add_argument("--airway-stl", type=Path, default=None, help="Lumen STL mesh (requires --inlet and --outlet)")
    parser.add_argument("--inlet", type=float, nargs=3, metavar=("X", "Y", "Z"), default=None)
    parser.add_argument("--outlet", type=float, nargs=3, metavar=("X", "Y", "Z"), default=None)
    parser.add_argument("--kl-loss", choices=["kirchhoff_fit", "legacy"], default="kirchhoff_fit")
    parser.add_argument("--trim-end-sections", default="0", help="Trim flagged end artifacts: 0 (default), auto, or N")
    parser.add_argument("--soft-tissue", type=Path, default=None, help="Path to soft_tissue_hypotheses.json")
    parser.add_argument("--soft-tissue-off", action="store_true", help="Disable soft tissue (rigid walls, t=0)")
    parser.add_argument("--sensitivity", type=int, default=None, metavar="N", help="Chain N-sample sensitivity sweep")
    parser.add_argument("--epb-table", type=Path, default=None, help="EPB comparative CSV table")
    parser.add_argument("--airway-centerline", type=Path, default=None, help="Looped airway centerline CSV with .source.json")
    parser.add_argument("--what-if-length-scale", type=float, default=None, help="Generate synthetic stretched-length what-if outputs")
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
            kl_loss_model=args.kl_loss,
            trim_end_sections=args.trim_end_sections,
            soft_tissue_path=args.soft_tissue,
            soft_tissue_off=args.soft_tissue_off,
            sensitivity_n=args.sensitivity,
            epb_table=args.epb_table,
            airway_centerline=args.airway_centerline,
            what_if_length_scale=args.what_if_length_scale,
        )
    except (GeometryInputError, SoftTissueError) as exc:
        print(f"ERROR ({type(exc).__name__}): {exc}", file=sys.stderr)
        raise SystemExit(2) from exc

    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()

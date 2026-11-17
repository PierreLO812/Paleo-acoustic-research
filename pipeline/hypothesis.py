"""État 1 — HYPOTHESIS: read taxon metadata, emit Antigravity engineering prompt."""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

from pipeline.paths import MODELS, OUTPUT, ROOT


def load_metadata(taxon_slug: str) -> dict:
    path = MODELS / taxon_slug / "metadata.json"
    if not path.exists():
        raise FileNotFoundError(f"Missing taxon metadata: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def list_taxa() -> list[str]:
    if not MODELS.exists():
        return []
    return sorted(
        p.name for p in MODELS.iterdir() if p.is_dir() and (p / "metadata.json").exists()
    )


def theoretical_fundamentals(meta: dict) -> dict:
    geom = meta["geometry"]
    soft = meta["soft_tissue_assumptions"]
    L = float(geom["conduit_length_mm"]) / 1000.0
    c = float(soft["speed_of_sound_m_s"])
    closed_open = [((2 * n - 1) * c) / (4.0 * L) for n in range(1, 4)]
    open_open = [(n * c) / (2.0 * L) for n in range(1, 4)]
    return {
        "L_m": L,
        "c_m_s": c,
        "f_closed_open_hz": closed_open,
        "f_open_open_hz": open_open,
        "preferred_formula": meta["acoustic_priors"]["theoretical_formula"],
    }


def build_antigravity_prompt(
    taxon_slug: str,
    *,
    run_id: str,
    retry_feedback: str | list | None = None,
) -> str:
    meta = load_metadata(taxon_slug)
    theory = theoretical_fundamentals(meta)
    out_dir = OUTPUT / run_id
    stl = meta["geometry"]["stl_path"]
    csv_path = meta["geometry"]["airway_features_csv"]
    formants = meta["acoustic_priors"]["expected_formants_hz"]
    q = meta["acoustic_priors"]["default_q_factors"]

    retry_block = ""
    if retry_feedback:
        fb_text = (
            json.dumps(retry_feedback, indent=2)
            if not isinstance(retry_feedback, str)
            else retry_feedback.strip()
        )
        retry_block = f"""
## RETRY FEEDBACK FROM PEER_REVIEW (fix this without altering measured values)
{fb_text}
"""

    return f"""# Antigravity — BUILD_AND_SIM for `{taxon_slug}`

You are the engineering agent for Pierre's Paleo-acoustic research repo.
Repo root: `{ROOT}`
Run id: `{run_id}`
Output directory (create it): `{out_dir}`

## Scientific goal
Produce a physically constrained vocalization sketch for **{meta['taxon']}** ({meta['specimen_id']}),
not a cinematic roar. Soft tissue is hypothesized; geometry comes from fossil STL / extracted airway.

## Raw anatomical parameters (from models/{taxon_slug}/metadata.json)
```json
{json.dumps(meta, indent=2)}
```

## Theoretical resonance anchors
- Conduit length L = {theory['L_m']:.6f} m
- c = {theory['c_m_s']} m/s
- Closed-open quarters: {[round(x, 2) for x in theory['f_closed_open_hz']]} Hz  (f = (2n-1)c / 4L)
- Open-open halves: {[round(x, 2) for x in theory['f_open_open_hz']]} Hz  (f = n c / 2L)
- Preferred model: `{theory['preferred_formula']}`
- Prior formants from metadata: {formants} Hz, Q≈{q} (targets only, never copy as measured peaks)

## Mandatory engineering steps
1. Use the independent solvers:
   - Geometry QC: `pipeline/geometry_qc.py`
   - Lossy Kirchhoff TMM: `pipeline/tmm_lossy.py`
   - Coupled KL waveguide + two-mass source: `pipeline/source_tract.py`
   - Audio-only measurement: `pipeline/measure.py`
2. Call `python -m pipeline.build_and_sim --taxon {taxon_slug} --run-id {run_id}`.
3. STL reference: `{stl}` ; airway table: `{csv_path}`
4. Do **not** invent fossil measurements. If a file is missing, fail clearly.

## Mandatory artifacts in `{out_dir}/`
- `source.wav`, `filtered.wav` (-1 dBFS peak normalized, no clipping)
- `spectrogram.png`, `transfer_function.png`, `tmm_lossy.png`, `kl_vs_tmm.png`, `geometry_qc.png`
- `metrics.json` (schema v2.0), `measured.json`, `tmm_lossy.json`, `state_*.json`, `report.md`, `build_log.txt`

{retry_block}
## Done criteria
All artifacts exist, `metrics.json` validates against `schemas/metrics.schema.json`, and `report.md` / `build_log.txt` are written.
"""


def write_hypothesis_bundle(
    taxon_slug: str,
    run_id: str,
    retry_feedback: str | list | None = None,
) -> dict:
    meta = load_metadata(taxon_slug)
    prompt = build_antigravity_prompt(taxon_slug, run_id=run_id, retry_feedback=retry_feedback)
    run_dir = OUTPUT / run_id
    run_dir.mkdir(parents=True, exist_ok=True)
    prompt_path = run_dir / "antigravity_prompt.md"
    prompt_path.write_text(prompt, encoding="utf-8")
    state = {
        "state": "HYPOTHESIS",
        "taxon": taxon_slug,
        "run_id": run_id,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "metadata_path": str(MODELS / taxon_slug / "metadata.json"),
        "prompt_path": str(prompt_path),
        "theory": theoretical_fundamentals(meta),
        "retry_feedback": retry_feedback,
    }
    (run_dir / "hypothesis.json").write_text(json.dumps(state, indent=2) + "\n", encoding="utf-8")
    return state

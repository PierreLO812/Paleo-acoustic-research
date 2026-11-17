# Step 2 Plan: Honest Acoustics, Lossy TMM, Non-Circular Review & Geometry QC

## Files to Modify
- `pipeline/build_and_sim.py` — Integrate M1–M7 (`--synth kl|legacy_iir`, `--f0`, `--temperature`, `--radiation`, `--seed`, `--airway-profile`, `--airway-stl`), schema v2.0, `report.md`, `build_log.txt`; remove feedback target blending.
- `pipeline/peer_review.py` — Replace circular target comparison with 3-way review (measured H(f)/LPC vs TMM vs analytic), spectral-shaping gate (+3 dB, >=2% F1 energy), provenance check, structured `retry_feedback`, and terminal states (`CONTRACT_OK`, `EXPLORATORY`, `SCIENTIFICALLY_VALIDATED`, `REVIEW_FAILED`).
- `pipeline/orchestrator.py` — Support new terminal states and `VALIDATED` -> `CONTRACT_OK` alias.
- `pipeline/hypothesis.py` — Update prompt contract notes for schema v2.0.
- `.gitignore` & `README.md` — Ignore `output/*` except one example run; document new CLI and modules.

## New Files to Create
- `pipeline/measure.py` — Audio-only YIN f0, LPC formants, Welch cross-spectral transfer function H(f), band energy fractions, source tilt, and CLI.
- `pipeline/tmm_lossy.py` — Independent lossy Kirchhoff TMM solver with Levine-Schwinger radiation, grid-independence check, and CLI.
- `pipeline/geometry_qc.py` — Geometry sanity gate (`BONE_DERIVED_EXPLORATORY` vs `AIRWAY_OK`), sidecar validation, `--airway-stl` handling, and `geometry_qc.png`.
- `pipeline/source_tract.py` — Two-mass/LF glottal source (-12 dB/oct tilt) and coupled Kelly-Lochbaum waveguide synthesis with KL-vs-TMM check.
- `schemas/metrics.schema.json` — JSON Schema (draft 2020-12) for `metrics.json` v2.0.
- `requirements.txt`, `tests/test_pipeline.py`, `tests/fixtures/old_metrics_eacbdf.json`.

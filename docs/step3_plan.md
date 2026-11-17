# Step 3 Plan: Physics Calibration, Soft-Tissue Hypotheses, Sensitivity & EPB/Loop Scaffolding

## Files Touched
- `pipeline/source_tract.py` — Per-section Kirchhoff frequency-dependent loss biquad fit, refined radiation reflection filter, 3-port yielding-wall junction, Q error in `compare_kl_vs_tmm`, and `tmm_fir` fallback.
- `pipeline/tmm_lossy.py` — Optional `WallModel` shunt admittance (`Y_w = P*l/Z_w`), `extra_damping_factor`, closed-end support, and `wall_resonance_closed_hz` (unchanged path when rigid).
- `pipeline/measure.py` — Both `F1 ± B1/2` (min 2 Welch bins) and legacy `F1 ± 60 Hz` shaping gains + band energy floor metrics.
- `pipeline/geometry_qc.py` — End-section jump QC (`end_section_artifact`), explicit `--trim-end-sections` helper, and looped centerline CSV+sidecar loader.
- `pipeline/peer_review.py` — `kl_vs_tmm_max_q_rel_err` (0.20), `_b3db` shaping gate, `_pm60` energy floor, `tmm_fir` non-independence guard, `synthetic` refusal, and strict `EXPLORATORY` cap.
- `pipeline/build_and_sim.py` — Runtime git provenance (`git_dirty`), new CLI flags (`--kl-loss`, `--soft-tissue*`, `--trim-end-sections`, `--sensitivity`, `--epb-table`, `--airway-centerline`, `--what-if-length-scale`), schema v2.1, and updated `report.md`.
- `models/parasaurolophus_ram14000/metadata.json` & `schemas/metrics.schema.json` — Mark `lining_thickness_mm` superseded; support metrics v2.0 and v2.1.

## New Files
- `pipeline/soft_tissue.py`, `schemas/soft_tissue.schema.json`, `models/parasaurolophus_ram14000/soft_tissue_hypotheses.json` — Strict hypothesis schema/loader and Steiner inward lumen offset.
- `pipeline/sensitivity.py` — Seeded LHS/OAT uncertainty sweep, Spearman/PRCC ranking, tornado/spaghetti/distribution plots, and min/median/max damping WAV scenarios.
- `pipeline/epb.py` & `data/epb/epb_comparative_template.csv` — Empty Extant Phylogenetic Bracket comparative CSV template and bracket loader/validator.
- `tests/fixtures/step2_reference_peaks.json` & `tests/test_step3.py` — Step-2 regression fixture and 12 step-3 verification tests.

## Milestone Order
1. **M0:** KL Kirchhoff losses + Q check, runtime git provenance, dual shaping bands, end-section QC/trim.
2. **M1:** Soft-tissue hypothesis config, Steiner lumen narrowing, yielding walls in TMM and KL.
3. **M2:** Deterministic LHS + OAT sensitivity sweep, plots, and 3 damping-scenario WAVs.
4. **M3–M4:** EPB comparative scaffolding and looped centerline / synthetic length-scale what-if.
5. **M5–M6:** Step-3 tests (`tests/test_step3.py`), schema v2.1 + report, clean-tree committed run, and push.

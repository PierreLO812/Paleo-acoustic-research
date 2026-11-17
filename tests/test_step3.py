"""M5 (Step 3) — Verification suite for yielding walls, soft-tissue hypotheses, sensitivity, EPB, and schema v2.1."""

from __future__ import annotations

import csv
import json
import shutil
import subprocess
from pathlib import Path

import jsonschema
import numpy as np
import pytest

from pipeline.build_and_sim import _git_provenance, run_build_and_sim
from pipeline.epb import EPBDataMissing, EPBError, derive_bracket, load_epb_table
from pipeline.geometry_qc import (
    apply_end_section_trim,
    evaluate_geometry_qc,
    load_area_profile,
    load_looped_centerline,
)
from pipeline.paths import MODELS, OUTPUT, ROOT
from pipeline.peer_review import resolve_terminal_state, review_metrics_only, review_run
from pipeline.sensitivity import run_sensitivity
from pipeline.soft_tissue import (
    SoftTissueError,
    apply_lumen_narrowing,
    load_hypotheses,
)
from pipeline.source_tract import (
    build_kl_waveguide,
    compare_kl_vs_tmm,
    extract_kl_peaks,
    kl_impulse_response,
)
from pipeline.tmm_lossy import Air, WallModel, check_grid_independence, solve_tmm


def _uniform_tube(length_mm: float = 250.0, radius_mm: float = 10.0, n_pts: int = 51):
    s_mm = np.linspace(0.0, length_mm, n_pts)
    area_mm2 = np.full(n_pts, np.pi * radius_mm**2, dtype=np.float64)
    dh_mm = np.full(n_pts, 2.0 * radius_mm, dtype=np.float64)
    return s_mm, area_mm2, dh_mm


def test_wall_tmm_closed_closed_wall_resonance():
    length_mm, radius_mm = 250.0, 10.0
    s_mm, area_mm2, dh_mm = _uniform_tube(length_mm, radius_mm)
    air = Air.from_temperature(20.0)
    mw = 10.0
    r_m = radius_mm * 1e-3
    p_m = 2.0 * np.pi * r_m
    a_m2 = np.pi * r_m**2
    fw_analytic = (1.0 / (2.0 * np.pi)) * np.sqrt(air.rho * (air.c**2) * p_m / (a_m2 * mw))

    res = solve_tmm(
        s_mm,
        area_mm2,
        dh_mm,
        air=air,
        fmin=10.0,
        fmax=2000.0,
        df=0.5,
        radiation="closed",
        losses=False,
        wall=WallModel(mass_per_area=mw, stiffness_per_area=0.0, resistance_per_area=0.0),
    )
    assert res["zin_peaks"], "Expected Z_in resonance peak for closed-closed mass-wall tube."
    f_lowest = res["zin_peaks"][0]["f_hz"]
    # No zero-frequency mode below f_w
    assert f_lowest > 0.5 * fw_analytic
    assert abs(f_lowest - fw_analytic) / fw_analytic < 0.03


def test_wall_tmm_closed_open_shift():
    length_mm, radius_mm = 250.0, 10.0
    s_mm, area_mm2, dh_mm = _uniform_tube(length_mm, radius_mm)
    air = Air.from_temperature(20.0)
    mw = 10.0
    r_m = radius_mm * 1e-3
    fw = (1.0 / (2.0 * np.pi)) * np.sqrt(air.rho * (air.c**2) * (2.0 * np.pi * r_m) / (np.pi * r_m**2 * mw))

    res_rigid = solve_tmm(s_mm, area_mm2, dh_mm, air=air, fmax=2500.0, df=0.5, radiation="unflanged", losses=True)
    res_mass = solve_tmm(
        s_mm,
        area_mm2,
        dh_mm,
        air=air,
        fmax=2500.0,
        df=0.5,
        radiation="unflanged",
        losses=True,
        wall=WallModel(mass_per_area=mw, stiffness_per_area=0.0, resistance_per_area=0.0),
    )
    res_damp = solve_tmm(
        s_mm,
        area_mm2,
        dh_mm,
        air=air,
        fmax=2500.0,
        df=0.5,
        radiation="unflanged",
        losses=True,
        wall=WallModel(mass_per_area=mw, stiffness_per_area=0.0, resistance_per_area=5000.0),
    )

    f1_r, f1_m = res_rigid["h_peaks"][0]["f_hz"], res_mass["h_peaks"][0]["f_hz"]
    f3_r, f3_m = res_rigid["h_peaks"][2]["f_hz"], res_mass["h_peaks"][2]["f_hz"]
    assert f1_m > f1_r
    f1_sq_pred = f1_r**2 + fw**2
    assert abs(f1_m**2 - f1_sq_pred) / f1_sq_pred < 0.05
    assert abs(f3_m - f3_r) / f3_r < abs(f1_m - f1_r) / f1_r

    b1_r, b1_d = res_rigid["h_peaks"][0]["bw_3db_hz"], res_damp["h_peaks"][0]["bw_3db_hz"]
    b3_r, b3_d = res_rigid["h_peaks"][2]["bw_3db_hz"], res_damp["h_peaks"][2]["bw_3db_hz"]
    assert (b1_d - b1_r) / b1_r > (b3_d - b3_r) / b3_r


def test_soft_tissue_off_regression():
    ref = json.loads((ROOT / "tests" / "fixtures" / "step2_reference_peaks.json").read_text(encoding="utf-8"))
    prof = load_area_profile(ROOT / "airway_features.csv")
    prof_zero = apply_lumen_narrowing(prof, 0.0)
    assert np.allclose(prof_zero["area_mm2"], prof["area_mm2"], rtol=1e-12)

    air = Air.from_temperature(20.0)
    tmm = check_grid_independence(
        prof_zero["s_mm"],
        prof_zero["area_mm2"],
        prof_zero["dh_mm"],
        air=air,
        fmin=10.0,
        fmax=4000.0,
        df=0.5,
        radiation="unflanged",
        losses=True,
        wall=None,
        extra_damping_factor=1.0,
    )
    for got, exp in zip(tmm["h_peaks"], ref["tmm_peaks_df025"]):
        assert np.isclose(got["f_hz"], exp["f_hz"], rtol=1e-9)
        assert np.isclose(got["q"], exp["q"], rtol=1e-9)

    wg = build_kl_waveguide(
        prof_zero["s_mm"],
        prof_zero["area_mm2"],
        air,
        radiation="unflanged",
        sr=44100,
        oversample=2,
        kl_loss_model="legacy",
        extra_damping_factor=1.0,
    )
    ir = kl_impulse_response(wg, ir_duration_s=0.8)
    _, _, kl_pks = extract_kl_peaks(ir, wg["fs_wg"], fmax_hz=3500.0)
    for got, exp in zip(kl_pks, ref["kl_peaks"]):
        assert np.isclose(got["f_hz"], exp["f_hz"], rtol=1e-9)
        assert np.isclose(got["q"], exp["q"], rtol=1e-9)


def test_lumen_offset_circle():
    r_mm = 10.0
    s_mm = np.linspace(0.0, 100.0, 11)
    area_mm2 = np.full(11, np.pi * r_mm**2)
    perim_mm = np.full(11, 2.0 * np.pi * r_mm)
    dh_mm = 4.0 * area_mm2 / perim_mm
    prof = {"s_mm": s_mm, "area_mm2": area_mm2, "perimeter_mm": perim_mm, "dh_mm": dh_mm}

    narrowed = apply_lumen_narrowing(prof, 2.0, min_area_fraction=0.05)
    assert narrowed["lumen_model"] == "uniform_inward_offset"
    assert np.allclose(narrowed["area_mm2"], np.pi * 8.0**2, rtol=1e-9)
    assert np.allclose(narrowed["perimeter_mm"], 2.0 * np.pi * 8.0, rtol=1e-9)
    assert np.allclose(narrowed["dh_mm"], 16.0, rtol=1e-9)

    with pytest.raises(SoftTissueError, match="occluded lumen"):
        apply_lumen_narrowing(prof, 9.9, min_area_fraction=0.05)
    with pytest.raises(SoftTissueError, match="occluded lumen"):
        apply_lumen_narrowing(prof, 12.0, min_area_fraction=0.05)


def test_hypothesis_validation(tmp_path: Path):
    base_path = MODELS / "parasaurolophus_ram14000" / "soft_tissue_hypotheses.json"
    cfg = load_hypotheses(base_path)
    assert abs(cfg.mucosa_thickness_mm.nominal_from_range - 2.25) < 1e-9
    assert abs(cfg.wall_mass_per_area.nominal_from_range - 10.0) < 1e-9
    base = json.loads(base_path.read_text(encoding="utf-8"))

    def _write_and_expect_err(mutator):
        data = json.loads(json.dumps(base))
        mutator(data)
        p = tmp_path / "bad_hyp.json"
        p.write_text(json.dumps(data), encoding="utf-8")
        with pytest.raises(SoftTissueError):
            load_hypotheses(p)

    # Missing source
    _write_and_expect_err(lambda d: d["parameters"]["mucosa_thickness_mm"].__setitem__("source", ""))
    # kind="hypothesis" with non-null value
    _write_and_expect_err(lambda d: d["parameters"]["mucosa_thickness_mm"].__setitem__("value", 2.0))
    # Missing range
    _write_and_expect_err(lambda d: d["parameters"]["mucosa_thickness_mm"].pop("range"))
    # Cited value without source_locator
    _write_and_expect_err(
        lambda d: d["parameters"]["mucosa_thickness_mm"].update(
            {"kind": "cited", "source": "Doe 2020", "source_locator": None, "value": 2.0}
        )
    )
    # lo > hi
    _write_and_expect_err(lambda d: d["parameters"]["mucosa_thickness_mm"].__setitem__("range", [4.0, 0.5]))
    # Wrong unit
    _write_and_expect_err(lambda d: d["parameters"]["mucosa_thickness_mm"].__setitem__("unit", "kg/m^2"))
    # Unknown key
    _write_and_expect_err(lambda d: d["parameters"]["mucosa_thickness_mm"].__setitem__("bogus_key", 123))


def test_sensitivity_deterministic(tmp_path: Path):
    d0a = tmp_path / "s0a"
    d0b = tmp_path / "s0b"
    d1 = tmp_path / "s1"

    sum0a = run_sensitivity("parasaurolophus_ram14000", n_samples=8, seed=0, out_dir=d0a, duration_s=0.18, render_audio=True, render_plots=True)
    sum0b = run_sensitivity("parasaurolophus_ram14000", n_samples=8, seed=0, out_dir=d0b, duration_s=0.18, render_audio=False, render_plots=False)
    sum1 = run_sensitivity("parasaurolophus_ram14000", n_samples=8, seed=1, out_dir=d1, duration_s=0.18, render_audio=False, render_plots=False)

    assert (d0a / "samples.csv").read_bytes() == (d0b / "samples.csv").read_bytes()
    for k in ("n_samples", "n_valid", "n_invalid", "invalid_reasons", "intervals_p05_median_p95", "oat_sensitivity", "importance"):
        assert sum0a[k] == sum0b[k]
    assert (d0a / "samples.csv").read_bytes() != (d1 / "samples.csv").read_bytes()
    assert sum0a["intervals_p05_median_p95"] != sum1["intervals_p05_median_p95"]

    for tag in ("min", "median", "max"):
        assert (d0a / f"scenario_damping_{tag}.wav").exists()
        assert (d0a / f"scenario_damping_{tag}.json").exists()
    for sc in sum0a["scenarios"]:
        assert sc["clipping_detected"] is False
        assert sc["audio_peak_dbfs"] <= -0.95


def test_kl_q_within_20pct():
    prof = load_area_profile(ROOT / "airway_features.csv")
    air = Air.from_temperature(20.0)
    tmm_rigid = solve_tmm(prof["s_mm"], prof["area_mm2"], prof["dh_mm"], air=air, fmax=3500.0, df=0.5, radiation="unflanged")
    wg_rigid = build_kl_waveguide(
        prof["s_mm"], prof["area_mm2"], air, dh_mm=prof["dh_mm"], radiation="unflanged", sr=44100, kl_loss_model="kirchhoff_fit"
    )
    ir_rigid = kl_impulse_response(wg_rigid, ir_duration_s=0.45)
    cmp_rigid = compare_kl_vs_tmm(wg_rigid, ir_rigid, {"peaks": tmm_rigid["h_peaks"]}, max_freq_hz=2000.0, max_rel_err=0.02)
    assert len(cmp_rigid["q_rel_errors"]) >= 3
    assert all(err <= 0.02 for err in cmp_rigid["rel_errors"][:3])
    assert all(q_err <= 0.20 for q_err in cmp_rigid["q_rel_errors"][:3])


@pytest.fixture(scope="module")
def _step3_what_if_run():
    run_id = "pytest_step3_shared"
    try:
        metrics = run_build_and_sim(
            "parasaurolophus_ram14000",
            run_id=run_id,
            synth_mode="kl",
            duration_s=0.35,
            seed=0,
            what_if_length_scale=2.0,
        )
        yield metrics, OUTPUT / run_id
    finally:
        shutil.rmtree(OUTPUT / run_id, ignore_errors=True)


def test_schema_v21(_step3_what_if_run):
    schema = json.loads((ROOT / "schemas" / "metrics.schema.json").read_text(encoding="utf-8"))
    step2_metrics = json.loads(
        (ROOT / "output" / "parasaurolophus_ram14000_20261117T104254Z_19c2be" / "metrics.json").read_text(encoding="utf-8")
    )
    jsonschema.validate(instance=step2_metrics, schema=schema)

    metrics, _ = _step3_what_if_run
    assert metrics["schema_version"] == "2.1"
    jsonschema.validate(instance=metrics, schema=schema)


def test_git_provenance(monkeypatch):
    def _mock_clean(cmd, **_kwargs):
        sub = " ".join(cmd)
        if "rev-parse --short" in sub:
            return subprocess.CompletedProcess(cmd, 0, stdout="abc1234\n", stderr="")
        if "rev-parse HEAD" in sub:
            return subprocess.CompletedProcess(cmd, 0, stdout="abc1234567890abcdef\n", stderr="")
        if "status --porcelain" in sub:
            # Only output/ untracked changes must NOT mark tree dirty
            return subprocess.CompletedProcess(cmd, 0, stdout="?? output/tmp_run/metrics.json\n", stderr="")
        return subprocess.CompletedProcess(cmd, 0, stdout="abc1234\n", stderr="")

    monkeypatch.setattr(subprocess, "run", _mock_clean)
    prov_clean = _git_provenance(ROOT)
    assert prov_clean["git_commit"] == "abc1234"
    assert prov_clean["git_dirty"] is False

    def _mock_dirty(cmd, **_kwargs):
        sub = " ".join(cmd)
        if "status --porcelain" in sub:
            return subprocess.CompletedProcess(cmd, 0, stdout=" M pipeline/measure.py\n", stderr="")
        if "describe" in sub:
            return subprocess.CompletedProcess(cmd, 0, stdout="abc1234-dirty\n", stderr="")
        return _mock_clean(cmd)

    monkeypatch.setattr(subprocess, "run", _mock_dirty)
    prov_dirty = _git_provenance(ROOT)
    assert prov_dirty["git_dirty"] is True
    assert prov_dirty["git_describe"] == "abc1234-dirty"


def test_end_section_flag_and_trim(tmp_path: Path):
    csv_path = tmp_path / "end_jump.csv"
    s = np.linspace(0.0, 200.0, 20)
    area = np.full(20, 500.0)
    perim = np.full(20, 80.0)
    area[-1] = 2000.0
    perim[-1] = 260.0
    with csv_path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["s_mm", "area_mm2", "perimeter_mm"])
        for si, ai, pi in zip(s, area, perim):
            w.writerow([f"{si:.3f}", f"{ai:.3f}", f"{pi:.3f}"])

    qc, prof = evaluate_geometry_qc(csv_path, source_mesh="mesh.stl", extraction_method="bone_slice")
    assert "end_section_artifact" in qc["flags"]
    assert qc["end_section"]["end_section_artifact"] is True
    assert qc["end_section"]["flagged_indices"] == [19]

    qc_0, prof_0, msgs_0 = apply_end_section_trim(qc, prof, "0")
    assert len(prof_0["s_mm"]) == 20
    assert qc_0["end_section"]["trimmed"] is False
    assert msgs_0 == []

    qc_t, prof_t, msgs_t = apply_end_section_trim(qc, prof, "auto")
    assert len(prof_t["s_mm"]) == 19
    assert qc_t["end_section"]["trimmed"] is True
    assert qc_t["end_section"]["trimmed_indices"] == [19]
    assert qc_t["L_mm"] < qc["L_mm"]
    assert msgs_t and "19" in msgs_t[0]


def test_epb_loader(tmp_path: Path):
    with pytest.raises(EPBDataMissing):
        load_epb_table(tmp_path / "nonexistent_epb.csv")

    with pytest.raises(EPBDataMissing, match="template only"):
        load_epb_table(ROOT / "data" / "epb" / "epb_comparative_template.csv")

    # Missing citation must raise EPBError
    no_cite = tmp_path / "no_cite.csv"
    with no_cite.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["taxon", "clade", "specimen_id", "region", "measure", "value", "unit", "n", "method", "citation", "doi_or_url", "page_or_figure", "notes"])
        w.writerow(["TestCroc", "Crocodylia", "TC-1", "nasal", "mucosa_thickness_mm", "1.5", "mm", "1", "CT", "", "10.0/fake", "p.1", ""])
    with pytest.raises(EPBError, match="citation"):
        load_epb_table(no_cite)

    # Single clade only -> bracket incomplete
    one_clade = tmp_path / "one_clade.csv"
    with one_clade.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["taxon", "clade", "specimen_id", "region", "measure", "value", "unit", "n", "method", "citation", "doi_or_url", "page_or_figure", "notes"])
        w.writerow(["TestCroc", "Crocodylia", "TC-1", "nasal", "mucosa_thickness_mm", "1.5", "mm", "1", "CT", "TEST FIXTURE – not data", "", "Fig. 1", ""])
    rows_one = load_epb_table(one_clade)
    with pytest.raises(EPBError, match="bracket incomplete"):
        derive_bracket(rows_one, "mucosa_thickness_mm")

    # Both Crocodylia and Aves -> complete bracket
    both = tmp_path / "both_clades.csv"
    with both.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["taxon", "clade", "specimen_id", "region", "measure", "value", "unit", "n", "method", "citation", "doi_or_url", "page_or_figure", "notes"])
        w.writerow(["TestCroc", "Crocodylia", "TC-1", "nasal", "mucosa_thickness_mm", "1.2", "mm", "1", "CT", "TEST FIXTURE – not data", "", "Fig. 1", ""])
        w.writerow(["TestBird", "Aves", "TB-1", "nasal", "mucosa_thickness_mm", "0.8", "mm", "2", "CT", "TEST FIXTURE – not data", "", "Table 2", ""])
    rows_both = load_epb_table(both)
    br = derive_bracket(rows_both, "mucosa_thickness_mm")
    assert br["range"] == [0.8, 1.2]
    assert br["kind"] == "epb_bracket"
    assert br["verified_by_human"] is False


def test_state_cap(tmp_path: Path, _step3_what_if_run):
    assert resolve_terminal_state(True, "BONE_DERIVED_EXPLORATORY") == "EXPLORATORY"
    assert resolve_terminal_state(False, "AIRWAY_OK") == "REVIEW_FAILED"

    metrics, run_dir = _step3_what_if_run
    assert metrics["state"] == "EXPLORATORY"
    assert metrics["what_if"]["synthetic"] is True
    wi_dir = run_dir / "what_if"
    assert not list(wi_dir.glob("state_*.json"))
    synth_json = wi_dir / "SYNTHETIC_L_x2_metrics.json"
    assert synth_json.exists()

    with pytest.raises(ValueError, match="synthetic"):
        review_run(synth_json)
    with pytest.raises(ValueError, match="synthetic"):
        review_run(wi_dir)
    with pytest.raises(ValueError, match="synthetic"):
        review_metrics_only({"synthetic": True})

    # Also test looped centerline loader on synthetic tmp CSV
    cl_csv = tmp_path / "loop_centerline.csv"
    with cl_csv.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["s_mm", "x_mm", "y_mm", "z_mm", "area_mm2", "perimeter_mm", "limb"])
        for i, si in enumerate(np.linspace(0.0, 520.0, 20)):
            limb = "ascending" if i < 10 else "descending"
            w.writerow([f"{si:.2f}", f"{si*0.5:.2f}", "10.0", "0.0", "80.0", "32.0", limb])
    cl_csv.with_suffix(".source.json").write_text(
        json.dumps({"origin": "CT segmentation", "reference": "Synthetic loop test", "author": "pytest"}),
        encoding="utf-8",
    )
    qc_loop, _ = load_looped_centerline(cl_csv, crest_length_mm=263.8)
    assert qc_loop["status"] == "AIRWAY_OK"
    assert qc_loop["looped_airway"]["loop_to_crest_ratio"] > 1.5

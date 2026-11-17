"""M6 — Fast (<60s), self-contained verification suite for M1–M7."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import jsonschema
import numpy as np
import pytest
from scipy import signal

from pipeline.build_and_sim import run_build_and_sim
from pipeline.geometry_qc import (
    GeometryInputError,
    evaluate_geometry_qc,
    extract_airway_from_stl,
)
from pipeline.measure import (
    band_energy_fractions,
    measure_f0,
    measure_formants_lpc,
    measure_transfer_function,
)
from pipeline.paths import ROOT
from pipeline.peer_review import normalize_legacy_state, review_metrics_only
from pipeline.source_tract import (
    build_kl_waveguide,
    check_clipping_raw,
    compare_kl_vs_tmm,
    kl_impulse_response,
    normalize_audio,
    synthesize_kl,
)
from pipeline.tmm_lossy import Air, check_grid_independence, solve_tmm


def _uniform_tube_arrays(length_mm: float = 250.0, radius_mm: float = 10.0, n_pts: int = 51):
    s_mm = np.linspace(0.0, length_mm, n_pts)
    area_mm2 = np.full(n_pts, np.pi * radius_mm**2, dtype=np.float64)
    dh_mm = np.full(n_pts, 2.0 * radius_mm, dtype=np.float64)
    return s_mm, area_mm2, dh_mm


def test_tmm_uniform_closed_open():
    length_m = 0.25
    radius_m = 0.010
    air = Air.from_temperature(20.0)
    s_mm, area_mm2, dh_mm = _uniform_tube_arrays(length_m * 1e3, radius_m * 1e3)

    # Case 1: ideal open end (Zr = 0)
    res_ideal = solve_tmm(s_mm, area_mm2, dh_mm, air=air, fmax=3000.0, df=0.5, radiation="ideal_open", losses=True)
    assert len(res_ideal["h_peaks"]) >= 4
    for n in range(1, 5):
        f_analytic = ((2 * n - 1) * air.c) / (4.0 * length_m)
        f_tmm = res_ideal["h_peaks"][n - 1]["f_hz"]
        assert abs(f_tmm - f_analytic) / f_analytic < 0.01

    # Case 2: unflanged Levine-Schwinger radiation (effective length L + 0.6133*r)
    res_unfl = solve_tmm(s_mm, area_mm2, dh_mm, air=air, fmax=3000.0, df=0.5, radiation="unflanged", losses=True)
    assert len(res_unfl["h_peaks"]) >= 4
    l_eff = length_m + 0.6133 * radius_m
    for n in range(1, 5):
        f_analytic = ((2 * n - 1) * air.c) / (4.0 * l_eff)
        f_tmm = res_unfl["h_peaks"][n - 1]["f_hz"]
        assert abs(f_tmm - f_analytic) / f_analytic < 0.01


def test_tmm_grid_independence():
    s_mm, area_mm2, dh_mm = _uniform_tube_arrays(260.0, 12.0)
    res = check_grid_independence(s_mm, area_mm2, dh_mm, air=Air(), fmax=3500.0, df=0.5, radiation="unflanged")
    assert res["grid_independence"]["passed"] is True
    assert res["grid_independence"]["max_freq_shift_pct"] < 0.2
    assert res["grid_independence"]["max_mag_shift_db"] < 0.5


def test_tmm_finite():
    s_mm, area_mm2, dh_mm = _uniform_tube_arrays(260.0, 18.0)
    res = solve_tmm(s_mm, area_mm2, dh_mm, air=Air(), fmax=4000.0, df=0.5, radiation="unflanged")
    assert np.all(np.isfinite(np.abs(res["h_transfer"])))
    assert np.all(np.isfinite(np.abs(res["z_in"])))


def _synth_allpole(f0: float, formants: list[float], bws: list[float], fs: int = 44100, dur: float = 2.0):
    n = int(dur * fs)
    src = np.zeros(n, dtype=np.float64)
    period = int(round(fs / f0))
    src[::period] = 1.0
    # Add tiny low-level white noise so H(f) between harmonics retains coherence > 0.8
    rng = np.random.default_rng(42)
    src_noisy = src + 0.01 * rng.normal(0.0, 1.0, size=n)
    y = src_noisy.copy()
    for f, bw in zip(formants, bws):
        r = np.exp(-np.pi * bw / fs)
        w = 2.0 * np.pi * f / fs
        a = np.array([1.0, -2.0 * r * np.cos(w), r**2])
        b = np.array([1.0 - 2.0 * r * np.cos(w) + r**2])
        y = signal.lfilter(b, a, y)
    return src_noisy, y


def test_lpc_synthetic():
    targets = [400.0, 1200.0, 2200.0]
    bws = [80.0, 120.0, 160.0]
    sr = 44100

    # f0 = 120 Hz -> LPC recovers within 5%, H(f) within 3%
    src120, y120 = _synth_allpole(120.0, targets, bws, fs=sr)
    f0_est = measure_f0(y120, sr)["measured_hz"]
    assert abs(f0_est - 120.0) / 120.0 < 0.03

    lpc120 = measure_formants_lpc(y120, sr, f0_hz=120.0)
    assert len(lpc120["f_hz"]) >= 3
    for f_meas, f_true in zip(lpc120["f_hz"][:3], targets):
        assert abs(f_meas - f_true) / f_true < 0.05

    tf120 = measure_transfer_function(src120, y120, sr)
    assert len(tf120["peaks"]) >= 3
    for f_true in targets:
        nearest = min(tf120["peaks"], key=lambda p: abs(p["f_hz"] - f_true))
        assert abs(nearest["f_hz"] - f_true) / f_true < 0.03

    # f0 = 60 Hz -> flagged low reliability
    _, y60 = _synth_allpole(60.0, targets, bws, fs=sr)
    lpc60 = measure_formants_lpc(y60, sr, f0_hz=60.0)
    assert lpc60["reliability"] == "low"


def test_transfer_function_white_noise():
    sr = 44100
    rng = np.random.default_rng(123)
    noise = rng.normal(0.0, 1.0, size=sr * 4)
    f_true = 650.0
    q_true = 12.0
    b, a = signal.iirpeak(f_true, Q=q_true, fs=sr)
    y = signal.lfilter(b, a, noise)

    tf = measure_transfer_function(noise, y, sr)
    assert len(tf["peaks"]) >= 1
    pk = min(tf["peaks"], key=lambda p: abs(p["f_hz"] - f_true))
    assert abs(pk["f_hz"] - f_true) / f_true < 0.02
    assert abs(pk["q"] - q_true) / q_true < 0.15


def test_band_fractions_sum_to_one():
    sr = 44100
    rng = np.random.default_rng(7)
    x = rng.normal(0.0, 1.0, size=sr * 2)
    fracs = band_energy_fractions(x, sr)
    assert abs(sum(fracs.values()) - 1.0) < 1e-6


def test_no_clipping(tmp_path: Path):
    s_mm, area_mm2, dh_mm = _uniform_tube_arrays(250.0, 12.0)
    air = Air()
    tmm_res = solve_tmm(s_mm, area_mm2, dh_mm, air=air, fmax=4000.0, df=0.5)
    tmm_summary = {"peaks": tmm_res["h_peaks"]}
    synth = synthesize_kl(s_mm, area_mm2, air, tmm_summary, duration_s=0.6, sr=44100, f0_hz=100.0, seed=0)

    for key in ("source", "filtered"):
        assert np.max(np.abs(synth[key])) <= 0.8913 + 1e-4
        assert check_clipping_raw(synth[key]) is False
    assert synth["audio"]["clipping_detected"] is False

    # Deliberately clipped signal must be detected as True
    clipped = np.array([0.0, 0.5, 1.0, -1.0], dtype=np.float64)
    assert check_clipping_raw(clipped) is True
    normed, _, norm_clip = normalize_audio(clipped, target_dbfs=-1.0)
    assert norm_clip is False and np.max(np.abs(normed)) <= 0.8913 + 1e-4


def test_review_fails_on_circular_metrics():
    fixture = ROOT / "tests" / "fixtures" / "old_metrics_eacbdf.json"
    old_metrics = json.loads(fixture.read_text(encoding="utf-8"))
    rev = review_metrics_only(old_metrics)
    assert rev["passed"] is False
    assert rev["terminal_state"] == "REVIEW_FAILED"
    assert len(rev["retry_feedback"]) >= 1
    assert any("circular metric" in iss for iss in rev["issues"])


def test_geometry_qc_flags_bone(tmp_path: Path):
    # 1. Bone-scale profile fixture (mean ~1371 mm2, max ~3000 mm2, Dh ~20.5 mm)
    bone_csv = tmp_path / "bone_profile.csv"
    s = np.linspace(0.0, 263.8, 50)
    area = np.linspace(400.0, 2342.0, 50)
    area[-1] = 2996.0
    perim = 4.0 * area / 20.5
    with bone_csv.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["distance_mm", "area_mm2", "perimeter_mm"])
        for si, ai, pi in zip(s, area, perim):
            w.writerow([f"{si:.3f}", f"{ai:.3f}", f"{pi:.3f}"])

    qc_bone, _ = evaluate_geometry_qc(
        bone_csv, source_mesh="RightSkeleton.stl", extraction_method="bone_mesh_skeleton_axis_slice"
    )
    assert qc_bone["status"] == "BONE_DERIVED_EXPLORATORY"
    assert len(qc_bone["flags"]) >= 2

    # 2. Valid synthetic airway CSV + sidecar passes as AIRWAY_OK
    air_csv = tmp_path / "lumen_airway.csv"
    s_air = np.linspace(0.0, 480.0, 40)
    r_air = np.full(40, 5.0)  # Dh = 10 mm <= 15 mm, area = 78.5 mm2 <= 1500 mm2
    with air_csv.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["s_mm", "area_mm2", "perimeter_mm"])
        for si, ri in zip(s_air, r_air):
            w.writerow([f"{si:.3f}", f"{np.pi * ri**2:.3f}", f"{2 * np.pi * ri:.3f}"])
    air_csv.with_suffix(".source.json").write_text(
        json.dumps({"origin": "CT segmentation", "reference": "Synthetic lumen test", "author": "pytest"}),
        encoding="utf-8",
    )
    qc_ok, _ = evaluate_geometry_qc(air_csv, require_sidecar=True)
    assert qc_ok["status"] == "AIRWAY_OK"

    # 3. Calling --airway-stl without inlet/outlet fails with GeometryInputError
    dummy_stl = tmp_path / "dummy.stl"
    dummy_stl.write_text("solid dummy\nendsolid dummy\n", encoding="utf-8")
    with pytest.raises(GeometryInputError, match="inlet"):
        extract_airway_from_stl(dummy_stl, inlet=None, outlet=None, out_csv=tmp_path / "out.csv")


def test_kl_matches_tmm():
    air = Air.from_temperature(20.0)

    # Case 1: uniform tube (L=0.25 m, r=10 mm)
    s_u, a_u, dh_u = _uniform_tube_arrays(250.0, 10.0)
    tmm_u = solve_tmm(s_u, a_u, dh_u, air=air, fmax=3000.0, df=0.5, radiation="unflanged")
    wg_u = build_kl_waveguide(s_u, a_u, air, radiation="unflanged", sr=44100)
    ir_u = kl_impulse_response(wg_u, ir_duration_s=0.8)
    cmp_u = compare_kl_vs_tmm(wg_u, ir_u, {"peaks": tmm_u["h_peaks"]}, max_freq_hz=2000.0, max_rel_err=0.02)
    assert cmp_u["passed"] is True

    # Case 2: two-section stepped tube (10 mm -> 20 mm radius)
    s_2 = np.linspace(0.0, 260.0, 52)
    r_2 = np.where(s_2 < 130.0, 10.0, 20.0)
    a_2 = np.pi * r_2**2
    dh_2 = 2.0 * r_2
    tmm_2 = solve_tmm(s_2, a_2, dh_2, air=air, fmax=3000.0, df=0.5, radiation="unflanged")
    wg_2 = build_kl_waveguide(s_2, a_2, air, radiation="unflanged", sr=44100)
    ir_2 = kl_impulse_response(wg_2, ir_duration_s=0.8)
    cmp_2 = compare_kl_vs_tmm(wg_2, ir_2, {"peaks": tmm_2["h_peaks"]}, max_freq_hz=2000.0, max_rel_err=0.02)
    assert cmp_2["passed"] is True


def test_metrics_schema(tmp_path: Path):
    # Create a validated synthetic airway CSV + sidecar and run full build_and_sim end-to-end
    air_csv = tmp_path / "valid_airway.csv"
    s_air = np.linspace(0.0, 263.8, 40)
    r_air = np.full(40, 6.0)  # Dh = 12 mm <= 15 mm
    with air_csv.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["s_mm", "area_mm2", "perimeter_mm"])
        for si, ri in zip(s_air, r_air):
            w.writerow([f"{si:.3f}", f"{np.pi * ri**2:.3f}", f"{2 * np.pi * ri:.3f}"])
    air_csv.with_suffix(".source.json").write_text(
        json.dumps({"origin": "CT segmentation", "reference": "pytest synthetic tube", "author": "pytest"}),
        encoding="utf-8",
    )

    run_id = "pytest_schema_check"
    metrics = run_build_and_sim(
        "parasaurolophus_ram14000",
        run_id=run_id,
        synth_mode="kl",
        duration_s=1.5,
        f0_hz=90.0,
        seed=0,
        airway_profile=air_csv,
    )
    schema = json.loads((ROOT / "schemas" / "metrics.schema.json").read_text(encoding="utf-8"))
    jsonschema.validate(instance=metrics, schema=schema)
    assert "formant_peaks_hz" not in metrics
    assert metrics["state"] == "SCIENTIFICALLY_VALIDATED"
    assert normalize_legacy_state("VALIDATED") == "CONTRACT_OK"

    # Clean up test output directory
    import shutil

    shutil.rmtree(ROOT / "output" / run_id, ignore_errors=True)

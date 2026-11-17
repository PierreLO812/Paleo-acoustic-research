"""M1 (Step 3) — Explicit soft-tissue hypothesis layer and Steiner inward lumen offset.

All unmeasured soft-tissue and physiological parameters live in a validated
``soft_tissue_hypotheses.json`` file with provenance, ranges, and uncertainty
scales. Subtracting mucosa thickness from a bone cross-section is strictly a
geometric sensitivity exercise on a bone-derived profile, never a true lumen
reconstruction, and cannot upgrade any terminal state above EXPLORATORY.
"""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np


class SoftTissueError(ValueError):
    """Raised on invalid hypothesis schema/provenance or physically occluded lumen."""


VALID_KINDS = {"cited", "analog_human", "analog_extant", "epb_bracket", "hypothesis"}
VALID_CONFIDENCE = {"high", "medium", "low", "unknown"}
VALID_SCALES = {"linear", "log"}
CITED_KINDS = {"cited", "analog_human", "analog_extant", "epb_bracket"}

TOP_LEVEL_KEYS = {"schema_version", "taxon_slug", "notes", "parameters"}
ENTRY_ALLOWED_KEYS = {
    "value",
    "range",
    "unit",
    "scale",
    "source",
    "source_locator",
    "model_source",
    "kind",
    "confidence",
    "verified_by_human",
    "used_in",
    "discrete_values",
    "notes",
}
ENTRY_REQUIRED_KEYS = {
    "value",
    "range",
    "unit",
    "source",
    "source_locator",
    "kind",
    "confidence",
    "notes",
}

EXPECTED_UNITS: dict[str, set[str]] = {
    "mucosa_thickness_mm": {"mm"},
    "wall_mass_per_area": {"kg/m^2"},
    "wall_stiffness_per_area": {"Pa/m", "N/m^3"},
    "wall_resistance_per_area": {"Pa*s/m"},
    "extra_damping_factor": {"dimensionless"},
    "f0_prior_hz": {"Hz"},
    "air_temperature_c": {"degC"},
    "lumen_min_area_fraction": {"dimensionless"},
}
REQUIRED_PARAMETERS = tuple(EXPECTED_UNITS.keys())


@dataclass(frozen=True)
class Hypothesis:
    """Single soft-tissue or physiological parameter with full provenance."""

    name: str
    value: float | None
    range: tuple[float, float]
    unit: str
    source: str
    source_locator: str | None
    kind: str
    confidence: str
    notes: str
    scale: str = "linear"
    model_source: str | None = None
    verified_by_human: bool = False
    used_in: tuple[str, ...] = ()
    discrete_values: tuple[float, ...] | None = None

    @property
    def nominal_from_range(self) -> float:
        lo, hi = float(self.range[0]), float(self.range[1])
        if self.scale == "log":
            return float(np.sqrt(lo * hi))
        return 0.5 * (lo + hi)

    @property
    def effective_value(self) -> float:
        if self.value is not None:
            return float(self.value)
        return self.nominal_from_range

    def to_dict(self) -> dict:
        return {
            "name": self.name,
            "value": None if self.value is None else float(self.value),
            "nominal_from_range": round(self.nominal_from_range, 6),
            "effective_value": round(self.effective_value, 6),
            "range": [float(self.range[0]), float(self.range[1])],
            "unit": self.unit,
            "scale": self.scale,
            "source": self.source,
            "source_locator": self.source_locator,
            "model_source": self.model_source,
            "kind": self.kind,
            "confidence": self.confidence,
            "verified_by_human": bool(self.verified_by_human),
            "used_in": list(self.used_in),
            "discrete_values": None if self.discrete_values is None else list(self.discrete_values),
            "notes": self.notes,
        }


@dataclass(frozen=True)
class SoftTissueConfig:
    """Validated bundle of all soft-tissue hypotheses for a run."""

    mucosa_thickness_mm: Hypothesis
    wall_mass_per_area: Hypothesis
    wall_stiffness_per_area: Hypothesis
    wall_resistance_per_area: Hypothesis
    extra_damping_factor: Hypothesis
    f0_prior_hz: Hypothesis
    air_temperature_c: Hypothesis
    lumen_min_area_fraction: Hypothesis
    schema_version: str = "1.0"
    taxon_slug: str = "parasaurolophus_ram14000"
    notes: str = ""
    source_path: str | None = None

    def parameters(self) -> dict[str, Hypothesis]:
        return {
            "mucosa_thickness_mm": self.mucosa_thickness_mm,
            "wall_mass_per_area": self.wall_mass_per_area,
            "wall_stiffness_per_area": self.wall_stiffness_per_area,
            "wall_resistance_per_area": self.wall_resistance_per_area,
            "extra_damping_factor": self.extra_damping_factor,
            "f0_prior_hz": self.f0_prior_hz,
            "air_temperature_c": self.air_temperature_c,
            "lumen_min_area_fraction": self.lumen_min_area_fraction,
        }

    def nominal_values(self) -> dict[str, float]:
        return {name: hyp.effective_value for name, hyp in self.parameters().items()}

    def citations_to_verify(self) -> list[dict]:
        """Collect all parameter citations and model equation citations pending human check."""
        out: list[dict] = []
        seen: set[tuple[str, str, str]] = set()
        for name, hyp in self.parameters().items():
            if hyp.kind in CITED_KINDS:
                key = (name, "parameter_value", hyp.source)
                if key not in seen:
                    seen.add(key)
                    out.append(
                        {
                            "parameter": name,
                            "role": "parameter_value",
                            "kind": hyp.kind,
                            "source": hyp.source,
                            "source_locator": hyp.source_locator,
                            "verified_by_human": bool(hyp.verified_by_human),
                        }
                    )
            if hyp.model_source:
                key = (name, "model_equation", hyp.model_source)
                if key not in seen:
                    seen.add(key)
                    out.append(
                        {
                            "parameter": name,
                            "role": "model_equation",
                            "kind": "model_source",
                            "source": hyp.model_source,
                            "source_locator": "equation form",
                            "verified_by_human": False,
                        }
                    )
        return out

    def with_parameter(self, name: str, hyp: Hypothesis) -> "SoftTissueConfig":
        if name not in EXPECTED_UNITS:
            raise SoftTissueError(f"Unknown soft-tissue parameter '{name}'.")
        return replace(self, **{name: hyp})

    def to_metrics_dict(
        self,
        *,
        lumen_model: str = "uniform_inward_offset",
        wall_model_active: bool = True,
        soft_tissue_enabled: bool = True,
    ) -> dict:
        params = self.parameters()
        return {
            "enabled": bool(soft_tissue_enabled),
            "source_file": self.source_path,
            "lumen_model": lumen_model,
            "wall_model": "yielding_lumped_shunt_RLC" if wall_model_active else "rigid",
            "state_effect": "none (hypotheses cannot upgrade state)",
            "nominal_values": {k: round(v.effective_value, 6) for k, v in params.items()},
            "parameters": {k: v.to_dict() for k, v in params.items()},
            "citations_to_verify": self.citations_to_verify(),
        }


def _validate_hypothesis_entry(name: str, entry: object) -> Hypothesis:
    if not isinstance(entry, dict):
        raise SoftTissueError(f"Parameter '{name}' must be a JSON object.")

    unknown_keys = set(entry.keys()) - ENTRY_ALLOWED_KEYS
    if unknown_keys:
        raise SoftTissueError(f"Unknown keys in parameter '{name}': {sorted(unknown_keys)}")

    missing_keys = ENTRY_REQUIRED_KEYS - set(entry.keys())
    if missing_keys:
        raise SoftTissueError(f"Missing required keys in parameter '{name}': {sorted(missing_keys)}")

    kind = entry["kind"]
    if not isinstance(kind, str) or kind not in VALID_KINDS:
        raise SoftTissueError(f"Invalid kind '{kind}' for '{name}'. Allowed: {sorted(VALID_KINDS)}")

    confidence = entry["confidence"]
    if not isinstance(confidence, str) or confidence not in VALID_CONFIDENCE:
        raise SoftTissueError(f"Invalid confidence '{confidence}' for '{name}'. Allowed: {sorted(VALID_CONFIDENCE)}")

    source = entry["source"]
    if not isinstance(source, str) or not source.strip():
        raise SoftTissueError(f"Parameter '{name}' must have a non-empty 'source'.")

    source_locator = entry["source_locator"]
    if source_locator is not None and not isinstance(source_locator, str):
        raise SoftTissueError(f"Parameter '{name}' source_locator must be a string or null.")

    val = entry["value"]
    if val is not None:
        if isinstance(val, bool) or not isinstance(val, (int, float)) or not np.isfinite(val):
            raise SoftTissueError(f"Parameter '{name}' value must be a finite number or null.")
        val = float(val)

    if kind == "hypothesis":
        if val is not None:
            raise SoftTissueError(f"Parameter '{name}' with kind='hypothesis' must have value=null (got {val}).")
        if source.strip().lower() != "hypothesis":
            raise SoftTissueError(f"Parameter '{name}' with kind='hypothesis' must have source='hypothesis'.")
    else:
        if source.strip().lower() == "hypothesis":
            raise SoftTissueError(f"Parameter '{name}' with kind='{kind}' cannot have source='hypothesis'.")
        if not source_locator or not source_locator.strip():
            raise SoftTissueError(f"Parameter '{name}' with kind='{kind}' requires a non-empty 'source_locator'.")
        if "verified_by_human" not in entry or not isinstance(entry["verified_by_human"], bool):
            raise SoftTissueError(f"Parameter '{name}' with kind='{kind}' requires boolean 'verified_by_human'.")

    rng_raw = entry["range"]
    if not isinstance(rng_raw, list) or len(rng_raw) != 2:
        raise SoftTissueError(f"Parameter '{name}' requires a 2-element 'range' [lo, hi].")
    if any(isinstance(v, bool) or not isinstance(v, (int, float)) or not np.isfinite(v) for v in rng_raw):
        raise SoftTissueError(f"Parameter '{name}' range values must be finite numbers.")
    lo, hi = float(rng_raw[0]), float(rng_raw[1])
    if lo > hi:
        raise SoftTissueError(f"Parameter '{name}' range has lo > hi ({lo} > {hi}).")

    scale = entry.get("scale", "linear")
    if not isinstance(scale, str) or scale not in VALID_SCALES:
        raise SoftTissueError(f"Parameter '{name}' scale must be one of {sorted(VALID_SCALES)}.")
    if scale == "log" and lo <= 0.0:
        raise SoftTissueError(f"Parameter '{name}' with scale='log' requires lo > 0 (got {lo}).")

    unit = entry["unit"]
    if not isinstance(unit, str) or unit not in EXPECTED_UNITS[name]:
        raise SoftTissueError(
            f"Wrong unit '{unit}' for parameter '{name}'. Expected one of {sorted(EXPECTED_UNITS[name])}."
        )

    if name == "extra_damping_factor" and lo < 1.0:
        raise SoftTissueError(f"extra_damping_factor must have range lo >= 1.0 (got {lo}).")
    if name == "mucosa_thickness_mm" and lo < 0.0:
        raise SoftTissueError(f"mucosa_thickness_mm must be >= 0 (got {lo}).")
    if name == "lumen_min_area_fraction" and (lo <= 0.0 or hi >= 1.0):
        raise SoftTissueError(f"lumen_min_area_fraction must lie strictly inside (0, 1) (got [{lo}, {hi}]).")

    notes = entry["notes"]
    if not isinstance(notes, str) or not notes.strip():
        raise SoftTissueError(f"Parameter '{name}' requires non-empty 'notes'.")

    model_source = entry.get("model_source")
    if model_source is not None and not isinstance(model_source, str):
        raise SoftTissueError(f"Parameter '{name}' model_source must be string or null.")

    verified_by_human = bool(entry.get("verified_by_human", False))
    used_in_raw = entry.get("used_in") or []
    if not isinstance(used_in_raw, list) or not all(isinstance(u, str) for u in used_in_raw):
        raise SoftTissueError(f"Parameter '{name}' used_in must be a list of strings.")

    disc_raw = entry.get("discrete_values")
    disc: tuple[float, ...] | None = None
    if disc_raw is not None:
        if not isinstance(disc_raw, list) or not disc_raw or not all(
            not isinstance(d, bool) and isinstance(d, (int, float)) and np.isfinite(d) for d in disc_raw
        ):
            raise SoftTissueError(f"Parameter '{name}' discrete_values must be a non-empty list of numbers.")
        disc = tuple(float(d) for d in disc_raw)

    return Hypothesis(
        name=name,
        value=val,
        range=(lo, hi),
        unit=unit,
        source=source,
        source_locator=source_locator,
        kind=kind,
        confidence=confidence,
        notes=notes,
        scale=scale,
        model_source=model_source,
        verified_by_human=verified_by_human,
        used_in=tuple(used_in_raw),
        discrete_values=disc,
    )


def load_hypotheses(path: Path) -> SoftTissueConfig:
    """Load and strictly validate a soft_tissue_hypotheses.json file."""
    path = Path(path)
    if not path.exists():
        raise SoftTissueError(f"Missing soft-tissue hypotheses file: '{path}'.")
    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise SoftTissueError(f"Invalid JSON in '{path}': {exc}") from exc

    if not isinstance(raw, dict):
        raise SoftTissueError(f"Top-level JSON in '{path}' must be an object.")

    unknown_top = set(raw.keys()) - TOP_LEVEL_KEYS
    if unknown_top:
        raise SoftTissueError(f"Unknown top-level keys in '{path}': {sorted(unknown_top)}")

    params_raw = raw.get("parameters")
    if not isinstance(params_raw, dict):
        raise SoftTissueError(f"Missing or invalid 'parameters' object in '{path}'.")

    unknown_params = set(params_raw.keys()) - set(REQUIRED_PARAMETERS)
    if unknown_params:
        raise SoftTissueError(f"Unknown parameter names in '{path}': {sorted(unknown_params)}")
    missing_params = set(REQUIRED_PARAMETERS) - set(params_raw.keys())
    if missing_params:
        raise SoftTissueError(f"Missing required parameters in '{path}': {sorted(missing_params)}")

    parsed = {name: _validate_hypothesis_entry(name, params_raw[name]) for name in REQUIRED_PARAMETERS}
    return SoftTissueConfig(
        **parsed,
        schema_version=str(raw.get("schema_version", "1.0")),
        taxon_slug=str(raw.get("taxon_slug", "unknown")),
        notes=str(raw.get("notes", "")),
        source_path=str(path),
    )


def rigid_soft_tissue_config() -> SoftTissueConfig:
    """Return an explicit rigid-wall / zero-thickness baseline config (no silent assumptions)."""
    def _rigid_h(name: str, val: float, unit: str, notes: str) -> Hypothesis:
        return Hypothesis(
            name=name,
            value=val,
            range=(val, val),
            unit=unit,
            source="rigid_baseline_off",
            source_locator="soft-tissue disabled",
            kind="cited",
            confidence="high",
            notes=notes,
            scale="linear",
            verified_by_human=True,
        )

    return SoftTissueConfig(
        mucosa_thickness_mm=_rigid_h("mucosa_thickness_mm", 0.0, "mm", "Soft tissue disabled (t=0)."),
        wall_mass_per_area=_rigid_h("wall_mass_per_area", 0.0, "kg/m^2", "Rigid walls (no yielding wall mass)."),
        wall_stiffness_per_area=_rigid_h("wall_stiffness_per_area", 0.0, "Pa/m", "Rigid walls (no wall stiffness)."),
        wall_resistance_per_area=_rigid_h("wall_resistance_per_area", 0.0, "Pa*s/m", "Rigid walls (no wall resistance)."),
        extra_damping_factor=_rigid_h("extra_damping_factor", 1.0, "dimensionless", "Smooth-wall Kirchhoff damping (1.0x)."),
        f0_prior_hz=_rigid_h("f0_prior_hz", 60.0, "Hz", "Metadata f0 prior."),
        air_temperature_c=_rigid_h("air_temperature_c", 20.0, "degC", "Baseline ambient air temperature (20 C)."),
        lumen_min_area_fraction=_rigid_h("lumen_min_area_fraction", 0.05, "dimensionless", "Occlusion guard."),
        taxon_slug="rigid",
        notes="Rigid-wall baseline with soft-tissue modeling disabled.",
    )


def apply_lumen_narrowing(
    prof: dict,
    thickness_mm: float | np.ndarray,
    *,
    min_area_fraction: float = 0.05,
) -> dict:
    """Apply Steiner inward parallel offset to bone cross-sections (sensitivity exercise, not a lumen).

    A_lumen = A_bone - P_bone * t + pi * t^2
    P_lumen = max(P_bone - 2 * pi * t, 0)
    Dh_lumen = 4 * A_lumen / P_lumen
    """
    area_bone = np.asarray(prof["area_mm2"], dtype=np.float64)
    dh_bone = np.asarray(prof["dh_mm"], dtype=np.float64)
    perim_bone = (
        np.asarray(prof["perimeter_mm"], dtype=np.float64)
        if prof.get("perimeter_mm") is not None
        else (4.0 * area_bone / dh_bone)
    )

    t_arr = np.asarray(thickness_mm, dtype=np.float64)
    if t_arr.ndim == 0:
        t_val = float(t_arr)
        is_uniform = True
        if t_val < 0.0:
            raise SoftTissueError(f"Mucosa thickness must be >= 0 mm (got {t_val}).")
        if t_val == 0.0:
            out = dict(prof)
            out["s_mm"] = np.asarray(prof["s_mm"], dtype=np.float64).copy()
            out["area_mm2"] = area_bone.copy()
            out["dh_mm"] = dh_bone.copy()
            if prof.get("perimeter_mm") is not None:
                out["perimeter_mm"] = perim_bone.copy()
            out["lumen_model"] = "rigid_unmodified"
            return out
        t_vec = np.full_like(area_bone, t_val, dtype=np.float64)
    else:
        is_uniform = False
        if t_arr.shape != area_bone.shape:
            raise SoftTissueError(f"Thickness profile shape {t_arr.shape} != area profile shape {area_bone.shape}.")
        if np.any(t_arr < 0.0):
            raise SoftTissueError("Spatially varying mucosa thickness must be >= 0 mm.")
        t_vec = t_arr

    perim_raw = perim_bone - 2.0 * np.pi * t_vec
    area_lumen = area_bone - perim_bone * t_vec + np.pi * (t_vec ** 2)
    perim_lumen = np.maximum(perim_raw, 0.0)
    frac = area_lumen / np.maximum(area_bone, 1e-12)

    bad_p = np.where(perim_raw <= 0.0)[0]
    if bad_p.size > 0:
        idx = int(bad_p[0])
        raise SoftTissueError(
            f"occluded lumen at section {idx} (s={float(prof['s_mm'][idx]):.2f} mm): "
            f"P_lumen={float(perim_raw[idx]):.3f} <= 0 mm at t={float(t_vec[idx]):.3f} mm."
        )
    bad_a = np.where(area_lumen <= 0.0)[0]
    if bad_a.size > 0:
        idx = int(bad_a[0])
        raise SoftTissueError(
            f"occluded lumen at section {idx} (s={float(prof['s_mm'][idx]):.2f} mm): "
            f"A_lumen={float(area_lumen[idx]):.3f} <= 0 mm^2 at t={float(t_vec[idx]):.3f} mm."
        )
    bad_frac = np.where(frac < float(min_area_fraction))[0]
    if bad_frac.size > 0:
        idx = int(bad_frac[0])
        raise SoftTissueError(
            f"occluded lumen at section {idx} (s={float(prof['s_mm'][idx]):.2f} mm): "
            f"A_lumen/A_bone={float(frac[idx]):.4f} < min_area_fraction={min_area_fraction:.4f}."
        )

    dh_lumen = 4.0 * area_lumen / perim_lumen
    out = dict(prof)
    out["s_mm"] = np.asarray(prof["s_mm"], dtype=np.float64).copy()
    out["area_mm2"] = area_lumen
    out["perimeter_mm"] = perim_lumen
    out["dh_mm"] = dh_lumen
    out["lumen_model"] = "uniform_inward_offset" if is_uniform else "spatially_varying_inward_offset"
    out["mucosa_thickness_mm_applied"] = float(t_vec[0]) if is_uniform else float(np.mean(t_vec))
    return out


def load_thickness_profile_csv(csv_path: Path, s_mm: np.ndarray) -> np.ndarray:
    """Load spatially varying mucosa thickness t(s) from CSV with mandatory 'source' column."""
    csv_path = Path(csv_path)
    if not csv_path.exists():
        raise SoftTissueError(f"Missing thickness profile CSV: '{csv_path}'.")
    with csv_path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None:
            raise SoftTissueError(f"Empty CSV '{csv_path}'.")
        fields = [f.strip() for f in reader.fieldnames if f is not None]
        rows = list(reader)
    pos_col = "s_mm" if "s_mm" in fields else ("distance_mm" if "distance_mm" in fields else None)
    if pos_col is None or "thickness_mm" not in fields or "source" not in fields:
        raise SoftTissueError(
            f"Thickness profile '{csv_path}' must contain columns: s_mm (or distance_mm), thickness_mm, source."
        )
    if len(rows) < 2:
        raise SoftTissueError(f"Thickness profile '{csv_path}' must contain at least 2 rows.")
    if any(not (r.get("source") or "").strip() for r in rows):
        raise SoftTissueError(f"Every row in thickness profile '{csv_path}' must specify a non-empty 'source'.")
    s_src = np.array([float(r[pos_col]) for r in rows], dtype=np.float64)
    t_src = np.array([float(r["thickness_mm"]) for r in rows], dtype=np.float64)
    if np.any(np.diff(s_src) <= 0.0) or np.any(t_src < 0.0):
        raise SoftTissueError(f"Invalid s or thickness values in '{csv_path}'.")
    return np.interp(np.asarray(s_mm, dtype=np.float64), s_src, t_src)

"""M5 — Geometry sanity gate and validated airway profile loader.

Detects whether a cross-sectional profile comes from a genuine nasal lumen
segmentation/endocast or from a raw bone skeleton mesh (which captures bone walls
rather than the looped internal airway).
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np


# Documented heuristic thresholds (never silent facts; emitted into metrics.json)
MAX_AIRWAY_DH_MEAN_MM = 15.0
MAX_AIRWAY_AREA_MM2 = 1500.0
MIN_LOOP_CREST_RATIO = 1.5
AREA_JUMP_FACTOR = 3.0
VALID_SIDECAR_ORIGINS = {"CT segmentation", "endocast", "literature"}


class GeometryInputError(FileNotFoundError, ValueError):
    """Raised when required airway geometry, sidecar provenance, or STL inlet/outlet is missing/invalid."""


def load_area_profile(csv_path: Path) -> dict:
    """Load and validate an airway cross-sectional area profile CSV.

    Required columns:
      - arc-length / position in mm: 's_mm' or 'distance_mm'
      - cross-sectional area in mm^2: 'area_mm2'
    Optional columns:
      - 'perimeter_mm'
      - 'hydraulic_diameter_mm'
    """
    csv_path = Path(csv_path)
    if not csv_path.exists():
        raise GeometryInputError(
            f"Missing airway profile CSV: '{csv_path}'. "
            "Expected a CSV with columns: s_mm (or distance_mm), area_mm2[, perimeter_mm]."
        )

    with csv_path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None:
            raise GeometryInputError(f"Empty CSV file or missing header in '{csv_path}'.")
        fields = [f.strip() for f in reader.fieldnames if f is not None]
        rows = list(reader)

    if len(rows) < 2:
        raise GeometryInputError(f"Airway profile '{csv_path}' must contain at least 2 cross-sections.")

    pos_candidates = [c for c in ("s_mm", "distance_mm") if c in fields]
    if len(pos_candidates) != 1:
        raise GeometryInputError(
            f"Cannot unambiguously identify arc-length column in '{csv_path}' (found fields: {fields}). "
            "Expected exactly one of 's_mm' or 'distance_mm' (units: mm)."
        )
    pos_col = pos_candidates[0]

    if "area_mm2" not in fields:
        raise GeometryInputError(
            f"Missing required area column 'area_mm2' (units: mm^2) in '{csv_path}' (found fields: {fields})."
        )

    try:
        s_mm = np.array([float(r[pos_col]) for r in rows], dtype=np.float64)
        area_mm2 = np.array([float(r["area_mm2"]) for r in rows], dtype=np.float64)
        perimeter_mm = (
            np.array([float(r["perimeter_mm"]) for r in rows], dtype=np.float64)
            if "perimeter_mm" in fields
            else None
        )
        dh_col = (
            np.array([float(r["hydraulic_diameter_mm"]) for r in rows], dtype=np.float64)
            if "hydraulic_diameter_mm" in fields
            else None
        )
    except (TypeError, ValueError) as exc:
        raise GeometryInputError(f"Non-numeric geometry values in '{csv_path}': {exc}") from exc

    if not (np.all(np.isfinite(s_mm)) and np.all(np.isfinite(area_mm2))):
        raise GeometryInputError(f"Non-finite (NaN/inf) values in '{csv_path}'.")
    if np.any(np.diff(s_mm) <= 0.0):
        raise GeometryInputError(
            f"Position column '{pos_col}' in '{csv_path}' must be strictly monotonically increasing (units: mm)."
        )
    if np.any(area_mm2 <= 0.0):
        raise GeometryInputError(f"All 'area_mm2' values in '{csv_path}' must be strictly positive (> 0 mm^2).")

    if perimeter_mm is not None:
        if not np.all(np.isfinite(perimeter_mm)) or np.any(perimeter_mm <= 0.0):
            raise GeometryInputError(f"All 'perimeter_mm' values in '{csv_path}' must be finite and > 0.")
        dh_mm = 4.0 * area_mm2 / perimeter_mm
    elif dh_col is not None:
        if not np.all(np.isfinite(dh_col)) or np.any(dh_col <= 0.0):
            raise GeometryInputError(f"All 'hydraulic_diameter_mm' values in '{csv_path}' must be finite and > 0.")
        dh_mm = dh_col
    else:
        dh_mm = 2.0 * np.sqrt(area_mm2 / np.pi)

    return {
        "s_mm": s_mm,
        "area_mm2": area_mm2,
        "perimeter_mm": perimeter_mm,
        "dh_mm": dh_mm,
        "pos_col": pos_col,
        "columns": fields,
    }


def _find_and_validate_sidecar(profile_csv: Path, sidecar_path: Path | None = None) -> tuple[Path, dict]:
    candidates = (
        [sidecar_path]
        if sidecar_path is not None
        else [profile_csv.with_suffix(".source.json"), Path(str(profile_csv) + ".source.json")]
    )
    existing = next((c for c in candidates if c.exists()), None)
    if existing is None:
        expected = candidates[0]
        raise GeometryInputError(
            f"Missing mandatory provenance sidecar '{expected}' for '{profile_csv}'. "
            'Expected JSON format: {"origin": "CT segmentation|endocast|literature", '
            '"reference": "...", "author": "..."}.'
        )
    try:
        data = json.loads(existing.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise GeometryInputError(f"Invalid JSON in sidecar '{existing}': {exc}") from exc

    missing = [k for k in ("origin", "reference", "author") if not data.get(k)]
    if missing:
        raise GeometryInputError(
            f"Sidecar '{existing}' missing required fields {missing}. "
            'Expected: {"origin": "CT segmentation|endocast|literature", "reference": "...", "author": "..."}.'
        )
    if data["origin"] not in VALID_SIDECAR_ORIGINS:
        raise GeometryInputError(
            f"Invalid sidecar origin '{data['origin']}' in '{existing}'. "
            f"Allowed: {sorted(VALID_SIDECAR_ORIGINS)}."
        )
    return existing, data


def evaluate_geometry_qc(
    profile_csv: Path,
    *,
    source_mesh: str = "unknown",
    extraction_method: str = "unknown",
    sidecar_path: Path | None = None,
    require_sidecar: bool = False,
    crest_length_mm: float | None = None,
) -> tuple[dict, dict]:
    """Compute geometry QC metrics, heuristic flags, and status for an airway profile."""
    profile_csv = Path(profile_csv)
    prof = load_area_profile(profile_csv)
    sidecar_data = None
    if require_sidecar:
        sc_path, sidecar_data = _find_and_validate_sidecar(profile_csv, sidecar_path)
        source_mesh = sidecar_data.get("reference", str(sc_path))
        extraction_method = sidecar_data["origin"]

    s_mm = prof["s_mm"]
    area_mm2 = prof["area_mm2"]
    dh_mm = prof["dh_mm"]

    l_mm = float(s_mm[-1] - s_mm[0])
    a_mean = float(np.mean(area_mm2))
    a_max = float(np.max(area_mm2))
    a_p95 = float(np.percentile(area_mm2, 95))
    dh_mean = float(np.mean(dh_mm))

    n_jumps = 0
    for i in range(len(area_mm2) - 1):
        lo = max(0, i - 2)
        hi = min(len(area_mm2), i + 3)
        local_med = float(np.median(area_mm2[lo:hi]))
        if abs(area_mm2[i + 1] - area_mm2[i]) > AREA_JUMP_FACTOR * max(local_med, 1e-6):
            n_jumps += 1

    flags: list[str] = []
    if dh_mean > MAX_AIRWAY_DH_MEAN_MM:
        flags.append(f"suspiciously wide for a nasal airway (Dh_mean={dh_mean:.2f} mm > {MAX_AIRWAY_DH_MEAN_MM} mm)")
    if a_max > MAX_AIRWAY_AREA_MM2:
        flags.append(
            f"area values compatible with bone/crest cross-sections, not lumen "
            f"(A_max={a_max:.1f} mm^2 > {MAX_AIRWAY_AREA_MM2} mm^2)"
        )
    is_bone_method = "bone" in extraction_method.lower() or (
        not require_sidecar and "skeleton" in source_mesh.lower()
    )
    if is_bone_method:
        flags.append("profile derived from bone mesh")

    if crest_length_mm is not None and crest_length_mm > 0.0:
        if l_mm < MIN_LOOP_CREST_RATIO * crest_length_mm:
            flags.append(
                f"L_mm ({l_mm:.1f} mm) < {MIN_LOOP_CREST_RATIO}x crest_length ({crest_length_mm:.1f} mm)"
            )
    elif not require_sidecar:
        flags.append("looped airway length unverified")

    if is_bone_method or dh_mean > MAX_AIRWAY_DH_MEAN_MM or a_max > MAX_AIRWAY_AREA_MM2:
        status = "BONE_DERIVED_EXPLORATORY"
    elif require_sidecar and sidecar_data is not None and not flags:
        status = "AIRWAY_OK"
    else:
        status = "UNKNOWN"

    qc = {
        "status": status,
        "flags": flags,
        "L_mm": round(l_mm, 3),
        "A_mean_mm2": round(a_mean, 3),
        "A_max_mm2": round(a_max, 3),
        "A_p95_mm2": round(a_p95, 3),
        "Dh_mean_mm": round(dh_mean, 3),
        "n_sections": int(len(s_mm)),
        "n_area_jumps": int(n_jumps),
        "source_mesh": str(source_mesh),
        "extraction_method": str(extraction_method),
        "sidecar": sidecar_data,
        "heuristic_thresholds": {
            "MAX_AIRWAY_DH_MEAN_MM": MAX_AIRWAY_DH_MEAN_MM,
            "MAX_AIRWAY_AREA_MM2": MAX_AIRWAY_AREA_MM2,
            "MIN_LOOP_CREST_RATIO": MIN_LOOP_CREST_RATIO,
            "AREA_JUMP_FACTOR": AREA_JUMP_FACTOR,
        },
    }
    return qc, prof


def extract_airway_from_stl(
    stl_path: Path,
    inlet: tuple[float, float, float] | None,
    outlet: tuple[float, float, float] | None,
    out_csv: Path,
) -> Path:
    """Extract an airway centerline profile from a lumen STL mesh using VMTK or shortest-path skeletonization."""
    stl_path = Path(stl_path)
    if not stl_path.exists():
        raise GeometryInputError(
            f"Missing airway STL file: '{stl_path}'. Provide a watertight lumen mesh STL."
        )
    if inlet is None or outlet is None or len(inlet) != 3 or len(outlet) != 3:
        raise GeometryInputError(
            "Both --inlet X Y Z and --outlet X Y Z coordinates (in mesh mm) are required with --airway-stl."
        )

    # Prefer VMTK if available
    try:
        import vmtk  # type: ignore  # noqa: F401
        has_vmtk = True
    except ImportError:
        has_vmtk = False

    if not has_vmtk:
        missing = []
        for mod in ("trimesh", "skimage", "networkx"):
            try:
                __import__(mod)
            except ImportError:
                missing.append(mod)
        if missing:
            raise GeometryInputError(
                f"Cannot extract lumen centerline from '{stl_path}': neither 'vmtk' nor fallback mesh stack "
                f"(missing packages: {missing}) is installed."
            )

    import networkx as nx  # type: ignore
    import trimesh  # type: ignore
    from scipy.ndimage import distance_transform_edt
    from skimage.morphology import skeletonize  # type: ignore

    mesh = trimesh.load(str(stl_path), force="mesh")
    pitch = float(np.max(mesh.extents)) / 120.0
    vox = mesh.voxelized(pitch=pitch).fill()
    mat = vox.matrix.astype(bool)
    _ = distance_transform_edt(mat)
    skel = skeletonize(mat)
    pts_idx = np.argwhere(skel)
    if len(pts_idx) < 2:
        raise GeometryInputError(f"Skeletonization of '{stl_path}' yielded fewer than 2 centerline voxels.")

    idx_to_node = {tuple(p): i for i, p in enumerate(pts_idx)}
    graph = nx.Graph()
    for i, p in enumerate(pts_idx):
        graph.add_node(i, xyz=vox.indices_to_points([p])[0])
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    if dx == 0 and dy == 0 and dz == 0:
                        continue
                    nb = (int(p[0] + dx), int(p[1] + dy), int(p[2] + dz))
                    if nb in idx_to_node:
                        j = idx_to_node[nb]
                        if i < j:
                            dist = float(np.linalg.norm(np.array([dx, dy, dz]) * pitch))
                            graph.add_edge(i, j, weight=dist)

    xyz_all = np.array([graph.nodes[i]["xyz"] for i in graph.nodes])
    start_node = int(np.argmin(np.linalg.norm(xyz_all - np.asarray(inlet), axis=1)))
    end_node = int(np.argmin(np.linalg.norm(xyz_all - np.asarray(outlet), axis=1)))
    if not nx.has_path(graph, start_node, end_node):
        raise GeometryInputError(
            f"No connected skeleton path between inlet {inlet} and outlet {outlet} inside '{stl_path}'."
        )
    path_nodes = nx.shortest_path(graph, start_node, end_node, weight="weight")
    centerline = np.array([graph.nodes[n]["xyz"] for n in path_nodes], dtype=np.float64)

    rows = []
    cum_s = 0.0
    for i in range(len(centerline)):
        if i > 0:
            cum_s += float(np.linalg.norm(centerline[i] - centerline[i - 1]))
        if i < len(centerline) - 1:
            tangent = centerline[i + 1] - centerline[i]
        else:
            tangent = centerline[i] - centerline[i - 1]
        norm = float(np.linalg.norm(tangent)) or 1.0
        tangent = tangent / norm
        sec = mesh.section(plane_origin=centerline[i], plane_normal=tangent)
        if sec is None:
            continue
        planar, _ = sec.to_2D()
        area = float(planar.area)
        perim = float(planar.length)
        if area > 0 and perim > 0:
            rows.append((cum_s, area, perim, 4.0 * area / perim))

    if len(rows) < 2:
        raise GeometryInputError(f"Failed to slice >=2 valid cross-sections along centerline of '{stl_path}'.")

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["s_mm", "area_mm2", "perimeter_mm", "hydraulic_diameter_mm"])
        for r in rows:
            w.writerow([f"{r[0]:.4f}", f"{r[1]:.4f}", f"{r[2]:.4f}", f"{r[3]:.4f}"])
    return out_csv


def plot_geometry_qc(qc: dict, prof: dict, out_png: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    s = prof["s_mm"]
    area = prof["area_mm2"]
    dh = prof["dh_mm"]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
    ax1.plot(s, area, color="steelblue", lw=1.8, label="Cross-sectional area A(s)")
    ax1.axhline(MAX_AIRWAY_AREA_MM2, color="crimson", ls="--", lw=1.2, label=f"Lumen max heuristic ({MAX_AIRWAY_AREA_MM2:.0f} mm²)")
    ax1.set_ylabel("Area (mm²)")
    ax1.set_title(f"Geometry Sanity QC — Status: {qc['status']} (L={qc['L_mm']:.1f} mm, A_mean={qc['A_mean_mm2']:.1f} mm²)")
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc="upper right", fontsize=8)

    ax2.plot(s, dh, color="darkorange", lw=1.8, label="Hydraulic diameter Dh(s)")
    ax2.axhline(MAX_AIRWAY_DH_MEAN_MM, color="crimson", ls="--", lw=1.2, label=f"Airway Dh_mean threshold ({MAX_AIRWAY_DH_MEAN_MM:.0f} mm)")
    ax2.set_xlabel("Arc length s (mm)")
    ax2.set_ylabel("Dh (mm)")
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc="upper right", fontsize=8)

    fig.tight_layout()
    fig.savefig(out_png, dpi=140)
    plt.close(fig)

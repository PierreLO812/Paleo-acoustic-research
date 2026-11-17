"""M3 (Step 3) — Extant Phylogenetic Bracket (EPB) comparative soft-tissue scaffolding.

Concept reference:
  Witmer, L. M. (1995). The Extant Phylogenetic Bracket and the importance of
  reconstructing soft tissues in fossils. In J. J. Thomason (ed.), Functional
  Morphology in Vertebrate Paleontology, pp. 19-33. Cambridge University Press.

Zero fake comparative data is included in the repository; ``data/epb/epb_comparative_template.csv``
contains headers only. Missing or template-only tables warn and skip without upgrading state.
"""

from __future__ import annotations

import csv
from pathlib import Path

from pipeline.soft_tissue import Hypothesis, SoftTissueConfig


EPB_COLUMNS = [
    "taxon",
    "clade",
    "specimen_id",
    "region",
    "measure",
    "value",
    "unit",
    "n",
    "method",
    "citation",
    "doi_or_url",
    "page_or_figure",
    "notes",
]
VALID_CLADES = {"Crocodylia", "Aves", "other"}
VALID_MEASURES = {
    "lumen_to_bone_area_ratio",
    "mucosa_thickness_mm",
    "wall_mass_per_area",
    "wall_stiffness_per_area",
    "wall_resistance_per_area",
}
WITMER_1995_REF = (
    "Witmer, L. M. (1995). The Extant Phylogenetic Bracket and the importance of reconstructing "
    "soft tissues in fossils. In Thomason, J. J. (ed.), Functional Morphology in Vertebrate "
    "Paleontology, pp. 19-33. Cambridge University Press."
)


class EPBError(ValueError):
    """Raised when an EPB table has invalid schema, missing citations, or an incomplete bracket."""


class EPBDataMissing(FileNotFoundError, EPBError):
    """Raised when the EPB comparative table file does not exist or contains only header rows."""


def load_epb_table(path: Path) -> list[dict]:
    """Load and validate an Extant Phylogenetic Bracket comparative CSV table."""
    path = Path(path)
    if not path.exists():
        raise EPBDataMissing(
            f"Missing EPB comparative table '{path}'. "
            f"Expected CSV with columns: {', '.join(EPB_COLUMNS)}."
        )

    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None:
            raise EPBError(f"Empty EPB CSV '{path}'. Expected columns: {', '.join(EPB_COLUMNS)}.")
        fields = [f.strip() for f in reader.fieldnames if f is not None]
        missing_cols = [c for c in EPB_COLUMNS if c not in fields]
        if missing_cols:
            raise EPBError(
                f"EPB CSV '{path}' missing required columns {missing_cols}. "
                f"Expected: {', '.join(EPB_COLUMNS)}."
            )
        rows = list(reader)

    if len(rows) == 0:
        raise EPBDataMissing("template only, no data")

    parsed: list[dict] = []
    for idx, row in enumerate(rows, start=2):
        clade = (row.get("clade") or "").strip()
        measure = (row.get("measure") or "").strip()
        citation = (row.get("citation") or "").strip()
        doi_or_url = (row.get("doi_or_url") or "").strip()
        page_or_figure = (row.get("page_or_figure") or "").strip()
        if clade not in VALID_CLADES:
            raise EPBError(f"Row {idx} in '{path}': invalid clade '{clade}'. Allowed: {sorted(VALID_CLADES)}.")
        if measure not in VALID_MEASURES:
            raise EPBError(f"Row {idx} in '{path}': invalid measure '{measure}'. Allowed: {sorted(VALID_MEASURES)}.")
        if not citation:
            raise EPBError(f"Row {idx} in '{path}': missing mandatory 'citation'.")
        if not (doi_or_url or page_or_figure):
            raise EPBError(
                f"Row {idx} in '{path}': must provide at least one of 'doi_or_url' or 'page_or_figure'."
            )
        try:
            value = float((row.get("value") or "").strip())
        except ValueError as exc:
            raise EPBError(f"Row {idx} in '{path}': non-numeric 'value'.") from exc
        n_raw = (row.get("n") or "").strip()
        n_val = int(n_raw) if n_raw else 1
        parsed.append(
            {
                "taxon": (row.get("taxon") or "").strip(),
                "clade": clade,
                "specimen_id": (row.get("specimen_id") or "").strip(),
                "region": (row.get("region") or "").strip(),
                "measure": measure,
                "value": value,
                "unit": (row.get("unit") or "").strip(),
                "n": n_val,
                "method": (row.get("method") or "").strip(),
                "citation": citation,
                "doi_or_url": doi_or_url,
                "page_or_figure": page_or_figure,
                "notes": (row.get("notes") or "").strip(),
            }
        )
    return parsed


def derive_bracket(table: list[dict], measure: str) -> dict:
    """Derive [min, max] EPB bracket requiring both Crocodylia and Aves observations."""
    if measure not in VALID_MEASURES:
        raise EPBError(f"Unsupported EPB measure '{measure}'. Allowed: {sorted(VALID_MEASURES)}.")

    sub = [r for r in table if r["measure"] == measure]
    croc_rows = [r for r in sub if r["clade"] == "Crocodylia"]
    aves_rows = [r for r in sub if r["clade"] == "Aves"]
    if not croc_rows or not aves_rows:
        missing = []
        if not croc_rows:
            missing.append("Crocodylia")
        if not aves_rows:
            missing.append("Aves")
        raise EPBError(f"bracket incomplete for '{measure}': missing clade(s) {missing}.")

    bracket_rows = croc_rows + aves_rows
    vals = [float(r["value"]) for r in bracket_rows]
    citations = list(dict.fromkeys(r["citation"] for r in bracket_rows))
    locators = list(
        dict.fromkeys(
            (r["page_or_figure"] or r["doi_or_url"]) for r in bracket_rows if (r["page_or_figure"] or r["doi_or_url"])
        )
    )
    return {
        "measure": measure,
        "range": [min(vals), max(vals)],
        "kind": "epb_bracket",
        "citations": citations,
        "source": "; ".join(citations),
        "source_locator": "; ".join(locators),
        "verified_by_human": False,
        "n_crocodylia": len(croc_rows),
        "n_aves": len(aves_rows),
    }


def apply_epb_table_to_config(
    cfg: SoftTissueConfig,
    epb_path: Path | None,
) -> tuple[SoftTissueConfig, dict, list[str]]:
    """Optionally load EPB table and apply complete brackets to SoftTissueConfig; missing/empty warns & skips."""
    warnings_out: list[str] = []
    if epb_path is None:
        return cfg, {
            "status": "not_provided",
            "path": None,
            "concept_reference": WITMER_1995_REF,
            "applied_measures": [],
        }, warnings_out

    epb_path = Path(epb_path)
    try:
        rows = load_epb_table(epb_path)
    except EPBDataMissing as exc:
        msg = f"EPB table skipped ({epb_path}): {exc}"
        warnings_out.append(msg)
        return cfg, {
            "status": "skipped_no_data",
            "reason": str(exc),
            "path": str(epb_path),
            "concept_reference": WITMER_1995_REF,
            "applied_measures": [],
        }, warnings_out

    updated = cfg
    applied: list[str] = []
    brackets: dict[str, dict] = {}
    params = cfg.parameters()
    for measure in sorted({r["measure"] for r in rows}):
        try:
            br = derive_bracket(rows, measure)
        except EPBError as exc:
            msg = f"EPB bracket skipped for '{measure}': {exc}"
            warnings_out.append(msg)
            continue
        brackets[measure] = br
        if measure in params:
            old_h = params[measure]
            new_h = Hypothesis(
                name=old_h.name,
                value=None,
                range=(float(br["range"][0]), float(br["range"][1])),
                unit=old_h.unit,
                source=br["source"],
                source_locator=br["source_locator"],
                kind="epb_bracket",
                confidence="low",
                notes=f"EPB bracket derived across Crocodylia (n={br['n_crocodylia']}) and Aves (n={br['n_aves']}). {old_h.notes}",
                scale=old_h.scale,
                model_source=old_h.model_source,
                verified_by_human=False,
                used_in=old_h.used_in,
                discrete_values=old_h.discrete_values,
            )
            updated = updated.with_parameter(measure, new_h)
            applied.append(measure)

    return updated, {
        "status": "applied" if applied else "loaded_no_complete_bracket",
        "path": str(epb_path),
        "concept_reference": WITMER_1995_REF,
        "n_rows": len(rows),
        "applied_measures": applied,
        "brackets": brackets,
    }, warnings_out

"""Rendering: ASCII canopy maps and inspection reports.

Reports are produced in three formats — Markdown for humans, JSON for
machines, and a compact text digest for terminals and CI logs.
"""

from __future__ import annotations

import json
import math
from collections.abc import Iterable, Sequence

from .detect import Crown
from .grid import Grid
from .inspection import Finding, InspectionReport, Severity

__all__ = [
    "HEIGHT_RAMP",
    "ascii_map",
    "findings_table",
    "render_json",
    "render_markdown",
    "render_text",
]

#: Characters ordered from lowest to highest canopy.
HEIGHT_RAMP = " .:-=+*o0#@"

SEVERITY_MARK = {
    Severity.INFO: "i",
    Severity.NOTICE: "-",
    Severity.WARNING: "!",
    Severity.CRITICAL: "X",
}


def ascii_map(
    grid: Grid,
    *,
    width: int = 72,
    ramp: str = HEIGHT_RAMP,
    markers: Sequence[tuple[float, float, str]] = (),
    aspect: float = 0.5,
    minimum: float | None = None,
    maximum: float | None = None,
) -> str:
    """Render a raster as text using a height ramp.

    The grid is reduced by block maximum, which keeps treetops visible instead
    of averaging them away.  ``markers`` are drawn on top in map coordinates.

    Args:
        grid: Raster to draw.
        width: Target character width.
        ramp: Characters from lowest to highest value.
        markers: ``(x, y, char)`` overlays in map coordinates.
        aspect: Character height/width ratio correction.
        minimum: Value mapped to the first ramp character (default: data min).
        maximum: Value mapped to the last ramp character (default: data max).
    """
    if width < 1:
        raise ValueError("width must be positive")
    if len(ramp) < 2:
        raise ValueError("ramp needs at least two characters")

    cols = min(width, grid.ncols)
    col_step = grid.ncols / cols
    row_step = col_step / aspect
    rows = max(1, round(grid.nrows / row_step))
    row_step = grid.nrows / rows

    stats = grid.stats()
    lo = stats.minimum if minimum is None else minimum
    hi = stats.maximum if maximum is None else maximum
    if lo is None or hi is None:
        return "(empty raster)"
    if math.isclose(hi, lo):
        hi = lo + 1.0

    canvas: list[list[str]] = []
    for r in range(rows):
        line: list[str] = []
        for c in range(cols):
            r0, r1 = int(r * row_step), max(int(r * row_step) + 1, int((r + 1) * row_step))
            c0, c1 = int(c * col_step), max(int(c * col_step) + 1, int((c + 1) * col_step))
            block = [
                grid.values[rr * grid.ncols + cc]
                for rr in range(r0, min(r1, grid.nrows))
                for cc in range(c0, min(c1, grid.ncols))
            ]
            valid = [v for v in block if v is not None]
            if not valid:
                line.append(" ")
                continue
            value = max(valid)
            position = (value - lo) / (hi - lo)
            index = min(len(ramp) - 1, max(0, int(position * (len(ramp) - 1) + 0.5)))
            line.append(ramp[index])
        canvas.append(line)

    for x, y, char in markers:
        cell = grid.cell_of(x, y)
        if cell is None:
            continue
        row, col = cell
        r = min(rows - 1, int(row / row_step))
        c = min(cols - 1, int(col / col_step))
        canvas[r][c] = char[0]

    body = "\n".join("".join(line) for line in canvas)
    legend = f"  low {lo:.1f} m [{ramp}] {hi:.1f} m high"
    return f"{body}\n{legend}"


def crown_markers(crowns: Iterable[Crown], char: str = "^") -> list[tuple[float, float, str]]:
    """Build :func:`ascii_map` markers for a set of crowns."""
    return [(crown.x, crown.y, char) for crown in crowns]


def findings_table(findings: Sequence[Finding]) -> str:
    """Render findings as a Markdown table."""
    if not findings:
        return "_No findings — the stand passed every rule._"
    lines = [
        "| Severity | Code | Subject | Finding |",
        "| --- | --- | --- | --- |",
    ]
    for finding in findings:
        detail = finding.detail.replace("|", "\\|")
        lines.append(
            f"| {finding.severity.label} | {finding.code} | `{finding.subject}` | "
            f"**{finding.title}.** {detail} |"
        )
    return "\n".join(lines)


def _kv_table(title: str, rows: Sequence[tuple[str, object]]) -> str:
    lines = [f"### {title}", "", "| Metric | Value |", "| --- | --- |"]
    for key, value in rows:
        lines.append(f"| {key} | {_fmt(value)} |")
    lines.append("")
    return "\n".join(lines)


def render_markdown(
    report: InspectionReport,
    *,
    title: str = "Silvispect stand inspection",
    chm: Grid | None = None,
    map_width: int = 72,
) -> str:
    """Render a full inspection report as Markdown."""
    metrics = report.metrics
    counts = report.counts()
    parts: list[str] = [f"# {title}", ""]
    parts.append(
        f"Inspected **{metrics.tree_count} trees** over **{report.area_ha:.3f} ha** "
        f"using {report.metrics_source} data. "
        f"Highest severity: **{report.max_severity.label}** "
        f"({len(report.findings)} findings: "
        + ", ".join(f"{count} {name}" for name, count in counts.items() if count)
        + ")."
        if report.findings
        else f"Inspected **{metrics.tree_count} trees** over "
        f"**{report.area_ha:.3f} ha** using {report.metrics_source} data. "
        "No findings."
    )
    parts.append("")

    parts.append("## Findings")
    parts.append("")
    parts.append(findings_table(report.findings))
    parts.append("")

    parts.append("## Stand metrics")
    parts.append("")
    parts.append(
        _kv_table(
            "Density and size",
            [
                ("Stems per hectare", metrics.stems_per_ha),
                ("Basal area (m2/ha)", metrics.basal_area_per_ha),
                ("Quadratic mean diameter (cm)", metrics.quadratic_mean_diameter_cm),
                ("Mean diameter (cm)", metrics.mean_diameter_cm),
                ("Stand density index", metrics.sdi),
                ("Gini of basal area", metrics.gini_basal_area),
            ],
        )
    )
    parts.append(
        _kv_table(
            "Height and yield",
            [
                ("Mean height (m)", metrics.mean_height_m),
                ("Lorey's mean height (m)", metrics.lorey_height_m),
                ("Dominant height (m)", metrics.dominant_height_m),
                ("Maximum height (m)", metrics.max_height_m),
                ("Volume (m3/ha)", metrics.volume_per_ha_m3),
                ("Above-ground biomass (t/ha)", metrics.biomass_per_ha_t),
            ],
        )
    )

    if metrics.species_shares:
        parts.append("### Species composition")
        parts.append("")
        parts.append("| Species | Share of stems |")
        parts.append("| --- | --- |")
        for species, share in metrics.species_shares.items():
            parts.append(f"| {species} | {share:.1%} |")
        parts.append("")

    if metrics.diameter_classes:
        parts.append("### Diameter distribution")
        parts.append("")
        parts.append("| Class (cm) | Stems |")
        parts.append("| --- | --- |")
        peak = max(metrics.diameter_classes.values())
        for label, count in metrics.diameter_classes.items():
            bar = "#" * max(1, round(32 * count / peak)) if count else ""
            parts.append(f"| {label} | {count} {bar} |")
        parts.append("")

    if report.canopy:
        parts.append("## Canopy")
        parts.append("")
        parts.append(
            _kv_table(
                "Cover and structure",
                [
                    ("Canopy cover", _pct(report.canopy.get("cover"))),
                    ("Open fraction", _pct(report.canopy.get("gap_fraction"))),
                    ("Rugosity (m)", report.canopy.get("rugosity_m")),
                    ("Mean canopy height (m)", report.canopy.get("mean_height_m")),
                    ("Maximum canopy height (m)", report.canopy.get("max_height_m")),
                    ("Mapped gaps", report.canopy.get("gap_count")),
                ],
            )
        )
        strata = report.canopy.get("strata")
        if isinstance(strata, dict) and strata:
            parts.append("### Vertical strata")
            parts.append("")
            parts.append("| Stratum | Share of cells |")
            parts.append("| --- | --- |")
            for name, share in strata.items():
                parts.append(f"| {name} | {share:.1%} |")
            parts.append("")

    if report.gaps:
        parts.append("### Canopy gaps")
        parts.append("")
        parts.append("| Gap | Area (m2) | Width (m) | Centre | Edge |")
        parts.append("| --- | --- | --- | --- | --- |")
        for gap in report.gaps[:20]:
            parts.append(
                f"| {gap.gap_id} | {gap.area:.0f} | {gap.width:.1f} | "
                f"({gap.centroid_x:.1f}, {gap.centroid_y:.1f}) | "
                f"{'yes' if gap.touches_edge else 'no'} |"
            )
        parts.append("")

    if report.detection:
        detection = report.detection
        parts.append("## Detection")
        parts.append("")
        parts.append(
            _kv_table(
                "Crowns from the canopy model",
                [
                    ("Detected trees", detection.get("tree_count")),
                    ("Detected stems per hectare", detection.get("density_per_ha")),
                ],
            )
        )
        if report.match:
            match = report.match
            parts.append(
                _kv_table(
                    "Agreement with the field reference",
                    [
                        ("Matched", match.get("matched")),
                        ("Omissions", match.get("omissions")),
                        ("Commissions", match.get("commissions")),
                        ("Recall", _pct(match.get("recall"))),
                        ("Precision", _pct(match.get("precision"))),
                        ("F1", match.get("f1")),
                        ("Height bias (m)", match.get("height_bias_m")),
                        ("Height RMSE (m)", match.get("height_rmse_m")),
                        ("Mean offset (m)", match.get("mean_offset_m")),
                    ],
                )
            )

    if report.allometry:
        allometry = report.allometry
        params = allometry.get("params")
        pretty = (
            ", ".join(f"{k} = {v}" for k, v in params.items()) if isinstance(params, dict) else "-"
        )
        parts.append("## Height-diameter model")
        parts.append("")
        parts.append(
            _kv_table(
                "Fitted curve",
                [
                    ("Family", allometry.get("model")),
                    ("Equation", f"`{allometry.get('equation')}`"),
                    ("Parameters", pretty),
                    ("Observations", allometry.get("n")),
                    ("RMSE (m)", allometry.get("rmse_m")),
                    ("R-squared", allometry.get("r_squared")),
                    ("Fitted to this stand", allometry.get("fitted")),
                ],
            )
        )
        by_species = allometry.get("by_species")
        if isinstance(by_species, dict) and by_species:
            parts.append("### Per-species curves")
            parts.append("")
            parts.append(
                "Fitted separately so that a minority species is not flagged as an "
                "outlier against the majority's curve."
            )
            parts.append("")
            parts.append("| Species | Parameters | n | RMSE (m) | R-squared |")
            parts.append("| --- | --- | --- | --- | --- |")
            for code, fit in by_species.items():
                fit_params = fit.get("params")
                shown = (
                    ", ".join(f"{k} = {v}" for k, v in fit_params.items())
                    if isinstance(fit_params, dict)
                    else "-"
                )
                parts.append(
                    f"| {code} | {shown} | {fit.get('n')} | "
                    f"{_fmt(fit.get('rmse_m'))} | {_fmt(fit.get('r_squared'))} |"
                )
            parts.append("")

    if chm is not None:
        parts.append("## Canopy map")
        parts.append("")
        parts.append("```")
        parts.append(ascii_map(chm, width=map_width))
        parts.append("```")
        parts.append("")

    return "\n".join(parts).rstrip() + "\n"


def render_text(report: InspectionReport) -> str:
    """Render a compact digest suitable for a terminal or a CI log."""
    metrics = report.metrics
    lines = [
        "Silvispect inspection",
        "=" * 21,
        f"area              {report.area_ha:.3f} ha ({report.metrics_source} data)",
        f"trees             {metrics.tree_count}",
        f"stems/ha          {metrics.stems_per_ha:.0f}",
        f"basal area/ha     {_fmt(metrics.basal_area_per_ha)} m2",
        f"qmd               {_fmt(metrics.quadratic_mean_diameter_cm)} cm",
        f"lorey height      {_fmt(metrics.lorey_height_m)} m",
        f"dominant height   {_fmt(metrics.dominant_height_m)} m",
        f"volume/ha         {_fmt(metrics.volume_per_ha_m3)} m3",
    ]
    if report.canopy:
        lines.append(f"canopy cover      {_pct(report.canopy.get('cover'))}")
        lines.append(f"gaps mapped       {report.canopy.get('gap_count')}")
    if report.match:
        lines.append(
            f"detection         recall {_pct(report.match.get('recall'))}, "
            f"precision {_pct(report.match.get('precision'))}"
        )
    lines.append("")
    if not report.findings:
        lines.append("no findings")
    else:
        counts = ", ".join(f"{count} {name}" for name, count in report.counts().items() if count)
        lines.append(f"findings ({counts}):")
        for finding in report.findings:
            mark = SEVERITY_MARK[finding.severity]
            lines.append(
                f"  [{mark}] {finding.code} {finding.subject:<16} {finding.title}: {finding.detail}"
            )
    return "\n".join(lines) + "\n"


def render_json(report: InspectionReport, *, indent: int = 2) -> str:
    """Render the report as JSON.

    ``allow_nan`` is off deliberately: ``NaN`` and ``Infinity`` are Python
    extensions that most JSON parsers reject, so a leak fails here rather than
    producing a document a consumer cannot read.
    """
    return json.dumps(report.as_dict(), indent=indent, sort_keys=False, allow_nan=False) + "\n"


def _fmt(value: object) -> str:
    if value is None:
        return "-"
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, float):
        if math.isnan(value) or math.isinf(value):
            return "-"
        return f"{value:,.2f}"
    return str(value)


def _pct(value: object) -> str:
    if not isinstance(value, (int, float)):
        return "-"
    return f"{float(value):.1%}"

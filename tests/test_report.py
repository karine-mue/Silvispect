"""Tests for ASCII maps and report rendering."""

from __future__ import annotations

import json

import pytest

from silvispect.grid import Grid
from silvispect.inspection import InspectionConfig, inspect_stand
from silvispect.inventory import Tree
from silvispect.report import (
    HEIGHT_RAMP,
    ascii_map,
    crown_markers,
    findings_table,
    render_json,
    render_markdown,
    render_text,
)


@pytest.fixture(scope="module")
def report(stand):
    return inspect_stand(chm=stand.chm, trees=stand.trees, area_ha=stand.area_ha)


def test_ascii_map_shape_and_legend():
    grid = Grid.filled(20, 40, 0.0, cellsize=1.0)
    grid.set(10, 20, 30.0)
    text = ascii_map(grid, width=40)
    lines = text.splitlines()
    assert all(len(line) == 40 for line in lines[:-1])
    assert len(lines) - 1 == 10  # aspect correction halves the row count
    assert lines[-1].startswith("  low 0.0 m")
    assert HEIGHT_RAMP[-1] in text  # the peak reaches the top of the ramp


def test_ascii_map_downsamples_by_maximum():
    grid = Grid.filled(4, 4, 0.0, cellsize=1.0)
    grid.set(0, 0, 25.0)
    assert ascii_map(grid, width=2).splitlines()[0][0] == HEIGHT_RAMP[-1]


def test_ascii_map_draws_nodata_as_blank():
    grid = Grid.filled(8, 8, 5.0, cellsize=1.0)
    for row in range(4):
        for col in range(4):
            grid.set(row, col, None)
    assert " " in ascii_map(grid, width=8).splitlines()[0]


def test_ascii_map_of_an_empty_raster():
    assert ascii_map(Grid.filled(4, 4, None)) == "(empty raster)"


def test_ascii_map_with_a_flat_raster_does_not_divide_by_zero():
    assert ascii_map(Grid.filled(4, 4, 7.0)).splitlines()[-1].count("7.0") == 1


def test_ascii_map_markers_are_drawn_and_clipped():
    grid = Grid.filled(20, 20, 5.0, cellsize=1.0)
    text = ascii_map(grid, width=20, markers=[(10.0, 10.0, "^"), (900.0, 900.0, "X")])
    assert "^" in text
    assert "X" not in text


def test_ascii_map_validates_arguments():
    grid = Grid.filled(4, 4, 1.0)
    with pytest.raises(ValueError, match="width"):
        ascii_map(grid, width=0)
    with pytest.raises(ValueError, match="ramp"):
        ascii_map(grid, ramp="#")


def test_ascii_map_honours_an_explicit_range():
    grid = Grid.filled(4, 4, 5.0)
    assert "10.0 m high" in ascii_map(grid, minimum=0.0, maximum=10.0)


def test_crown_markers(stand):
    from silvispect.detect import detect_trees

    crowns = detect_trees(stand.chm).crowns
    markers = crown_markers(crowns, "*")
    assert len(markers) == len(crowns)
    assert markers[0][2] == "*"


def test_findings_table_empty():
    assert "No findings" in findings_table([])


def test_findings_table_escapes_pipes():
    from silvispect.inspection import Finding, Severity

    finding = Finding("SV999", Severity.WARNING, "Odd", "a | b", subject="stand")
    assert "a \\| b" in findings_table([finding])


def test_render_text_digest(report):
    text = render_text(report)
    assert text.startswith("Silvispect inspection")
    assert "stems/ha" in text
    assert "canopy cover" in text
    assert "detection" in text
    assert text.endswith("\n")


def test_render_text_without_findings():
    clean = inspect_stand(
        trees=[Tree(f"t{i}", i, i, species="PIAB", dbh_cm=30.0, height_m=25.0) for i in range(20)],
        area_ha=0.05,
        config=InspectionConfig(
            max_stems_per_ha=100_000.0,
            max_sdi=100_000.0,
            min_gini=0.0,
            min_shannon=0.0,
            max_species_share=1.0,
        ),
    )
    assert "no findings" in render_text(clean)


def test_render_markdown_sections(report, stand):
    text = render_markdown(report, chm=stand.chm)
    for heading in (
        "# Silvispect stand inspection",
        "## Findings",
        "## Stand metrics",
        "### Species composition",
        "### Diameter distribution",
        "## Canopy",
        "## Detection",
        "## Height-diameter model",
        "## Canopy map",
    ):
        assert heading in text, f"missing section {heading}"
    assert "```" in text  # the map is fenced
    assert text.endswith("\n")


def test_render_markdown_lists_per_species_curves(report):
    text = render_markdown(report)
    assert "### Per-species curves" in text
    assert "| Species | Parameters | n | RMSE (m) | R-squared |" in text
    for code in ("PIAB", "FASY", "BEPE"):
        assert code in text


def test_render_markdown_without_a_map(report):
    assert "## Canopy map" not in render_markdown(report)


def test_render_markdown_custom_title(report):
    assert render_markdown(report, title="Compartment 7").startswith("# Compartment 7")


def test_render_markdown_of_a_clean_report():
    clean = inspect_stand(
        trees=[Tree(f"t{i}", i, i, species="PIAB", dbh_cm=30.0, height_m=25.0) for i in range(20)],
        area_ha=0.05,
        config=InspectionConfig(
            max_stems_per_ha=100_000.0,
            max_sdi=100_000.0,
            min_gini=0.0,
            min_shannon=0.0,
            max_species_share=1.0,
        ),
    )
    text = render_markdown(clean)
    assert "No findings" in text


def test_render_json_is_valid(report):
    payload = json.loads(render_json(report))
    assert payload["summary"]["max_severity"] in {"info", "notice", "warning", "critical"}
    assert payload["metrics"]["tree_count"] == report.metrics.tree_count
    assert isinstance(payload["findings"], list)


def test_render_json_has_no_nan(report):
    """NaN is not valid JSON, so the renderer must never emit it."""
    text = render_json(report)
    assert "NaN" not in text
    assert "Infinity" not in text


def test_render_json_without_field_data_has_no_nan(stand):
    """A default (unfitted) curve has NaN diagnostics that must not leak."""
    text = render_json(inspect_stand(chm=stand.chm))
    assert "NaN" not in text
    payload = json.loads(text)
    assert payload["allometry"]["rmse_m"] is None
    assert payload["allometry"]["fitted"] is False

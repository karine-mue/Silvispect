"""Regression tests — Lanes C & D (metrics/report/inventory/cli/optimize).

Issues #8, #9, #10, #12, #18, #21, #22, #23, #26, #29.
"""

import io
import json
import math

import pytest

from silvispect.canopy import CanopyGap
from silvispect.cli import main
from silvispect.detect import detect_trees
from silvispect.grid import Grid
from silvispect.inspection import InspectionConfig, InspectionReport, inspect_stand
from silvispect.inventory import (
    InventoryError,
    Tree,
    match_trees,
    parse_trees,
)
from silvispect.metrics import (
    MetricsError,
    StandMetrics,
    diameter_distribution,
    dominant_height,
    stand_metrics,
)
from silvispect.optimize import nelder_mead
from silvispect.report import ascii_map, render_markdown


def tree(tid, dbh=None, height=None, species=""):
    return Tree(tid, 0.0, 0.0, species=species, dbh_cm=dbh, height_m=height)


# --- #8: SDI population --------------------------------------------------------

def test_sdi_uses_full_stem_count():
    trees = [tree(f"M{i}", dbh=25.0, height=20.0) for i in range(100)]
    trees += [tree(f"U{i}", dbh=None, height=20.0) for i in range(100)]
    metrics = stand_metrics(trees, 1.0)
    assert metrics.stems_per_ha == 200.0
    assert metrics.sdi == pytest.approx(200.0)


# --- #12: unlabelled species -----------------------------------------------------

def test_species_composition_ignores_unlabelled_stems():
    trees = [tree(f"A{i}", species="") for i in range(50)]
    trees += [tree(f"B{i}", species="PIAB") for i in range(50)]
    metrics = stand_metrics(trees, 1.0)
    assert metrics.species_shares == {"PIAB": 1.0}


def test_diversity_is_unknown_when_no_species_data():
    trees = [tree(f"A{i}", dbh=20.0, height=18.0) for i in range(20)]
    metrics = stand_metrics(trees, 1.0)
    assert metrics.species_shares == {}
    assert metrics.shannon is None
    assert metrics.simpson is None


# --- #23: as_dict non-finite handling -----------------------------------------------

def test_metrics_as_dict_maps_nonfinite_to_none():
    metrics = stand_metrics([tree("BIG", dbh=250.0, height=1e308)], 1.0)
    payload = metrics.as_dict()
    assert payload["volume_per_ha_m3"] is None
    json.dumps(payload, allow_nan=False)


# --- #29: dominant height rounding / finite checks -------------------------------------

def test_dominant_height_rounds_half_up():
    trees = []
    for i in range(20):
        dbh = 40.0 - i
        trees.append(tree(f"T{i}", dbh=dbh, height=100.0 - dbh))
    # 100 * 0.125 = 12.5 stems: half-up keeps 13 (heights 60..72, mean 66).
    assert dominant_height(trees, 0.125) == pytest.approx(66.0)


def test_diameter_distribution_rejects_nonfinite_width():
    with pytest.raises(MetricsError):
        diameter_distribution([tree("A", dbh=20.0)], class_width=float("nan"))
    with pytest.raises(MetricsError):
        diameter_distribution([tree("A", dbh=20.0)], class_width=float("inf"))


def test_stand_metrics_rejects_infinite_area():
    with pytest.raises(MetricsError):
        stand_metrics([tree("A", dbh=20.0, height=18.0)], float("inf"))


def test_report_notes_truncated_gap_list():
    gaps = [
        CanopyGap(
            gap_id=i + 1,
            cells=((0, 0),),
            area=100.0,
            centroid_x=1.0,
            centroid_y=2.0,
            mean_height=1.0,
            max_height=2.0,
            touches_edge=False,
            inscribed_radius=5.0,
        )
        for i in range(25)
    ]
    report = InspectionReport(
        findings=[],
        metrics=StandMetrics.unknown(1.0),
        metrics_source="field",
        config=InspectionConfig(),
        area_ha=1.0,
        gaps=gaps,
    )
    markdown = render_markdown(report)
    assert "of 25 gaps" in markdown


# --- #18: detection density uses the inspection's plot area --------------------------------

def test_inspection_detection_density_uses_plot_area():
    chm = Grid.from_rows([[5.0, 5.0, 5.0]] * 3, cellsize=1.0)
    crowns = detect_trees(chm).crowns
    assert crowns
    report = inspect_stand(chm=chm, area_ha=0.5)
    assert report.detection is not None
    assert report.detection["density_per_ha"] == round(len(crowns) / 0.5, 2)


# --- #23: ascii_map robustness ---------------------------------------------------------------

def test_ascii_map_skips_empty_marker_char():
    grid = Grid.from_rows([[5.0]])
    rendered = ascii_map(grid, markers=[(0.5, 0.5, "")])
    assert "low" in rendered


def test_ascii_map_rejects_nonpositive_aspect():
    grid = Grid.from_rows([[5.0]])
    with pytest.raises(ValueError):
        ascii_map(grid, aspect=0.0)


def test_ascii_map_rejects_inverted_range():
    grid = Grid.from_rows([[5.0]])
    with pytest.raises(ValueError):
        ascii_map(grid, minimum=30.0, maximum=10.0)


# --- #22: empty CSV / Tree coercion -------------------------------------------------------

def test_parse_empty_csv_raises():
    with pytest.raises(InventoryError):
        parse_trees("")


def test_tree_coerces_numeric_strings():
    coerced = Tree("T", "10", "20")
    assert coerced.x == 10.0
    assert isinstance(coerced.x, float)


def test_tree_rejects_bool_dbh():
    with pytest.raises(InventoryError):
        Tree("T", 0.0, 0.0, dbh_cm=True)


def test_tree_rejects_nonnumeric_coordinate():
    with pytest.raises(InventoryError):
        Tree("T", "abc", 0.0)


# --- #9: generated ids -----------------------------------------------------------------------

def test_generated_ids_never_collide_with_explicit_ids():
    trees = parse_trees("tree_id,x,y\n2,10,10\n,20,20\n,30,30\n")
    ids = [t.tree_id for t in trees]
    assert len(ids) == len(set(ids)) == 3


def test_generated_ids_skip_later_explicit_ids():
    trees = parse_trees("tree_id,x,y\n,10,10\n1,20,20\n")
    ids = [t.tree_id for t in trees]
    assert len(ids) == len(set(ids)) == 2


# --- #26: match_trees tolerance -----------------------------------------------------------------

def test_match_trees_rejects_nan_tolerance():
    with pytest.raises(InventoryError):
        match_trees([], [], tolerance=float("nan"))


# --- #23: nelder_mead -inf -------------------------------------------------------------------------

def test_nelder_mead_treats_negative_inf_as_worst():
    result = nelder_mead(lambda p: float("-inf"), [1.0], max_iterations=50)
    assert result.fun == math.inf


# --- #10 / #21: CLI ---------------------------------------------------------------------------------

def test_cli_match_counts_only_observed_stems(tmp_path):
    chm = Grid.from_rows([
        [5.0, 5.0, None],
        [5.0, 5.0, None],
        [5.0, 5.0, None],
    ], cellsize=1.0)
    chm_path = tmp_path / "chm.asc"
    chm.write(chm_path)
    inventory = tmp_path / "trees.csv"
    inventory.write_text("tree_id,x,y\nA,0.5,2.5\nB,2.5,2.5\n", encoding="utf-8")
    out = io.StringIO()
    code = main(["match", str(inventory), "--chm", str(chm_path), "--json"], out=out)
    assert code == 0
    payload = json.loads(out.getvalue())
    assert payload["recall"] == 1.0


def test_cli_directory_input_exits_2(tmp_path):
    code = main(["inspect", "--chm", str(tmp_path)], out=io.StringIO())
    assert code == 2


def test_cli_negative_precision_is_a_clean_usage_error(tmp_path):
    dsm = Grid.from_rows([[10.0]])
    dtm = Grid.from_rows([[0.0]])
    dsm_path = tmp_path / "dsm.asc"
    dtm_path = tmp_path / "dtm.asc"
    dsm.write(dsm_path)
    dtm.write(dtm_path)
    with pytest.raises(SystemExit) as excinfo:
        main(["chm", "--dsm", str(dsm_path), "--dtm", str(dtm_path),
              "-o", str(tmp_path / "out.asc"), "--precision", "-1"],
             out=io.StringIO())
    assert excinfo.value.code == 2

"""Regression tests — Lane B (inspection/allometry/canopy).

Issues #7, #15, #16, #17, #19, #20, #24, #26, #27.
"""

import math

import pytest

from silvispect.allometry import (
    MODELS,
    AllometryError,
    default_model,
    fit_by_species,
    fit_height_diameter,
    residual_scores,
)
from silvispect.canopy import find_gaps, rugosity
from silvispect.grid import Grid
from silvispect.inspection import (
    InspectionConfig,
    InspectionContext,
    check_height_anomalies,
    inspect_stand,
)
from silvispect.inventory import Tree
from silvispect.metrics import stand_metrics


def on_curve(tree_id, dbh, offset=0.0, species="", status="live"):
    model = default_model("")
    return Tree(tree_id, 0.0, 0.0, species=species, dbh_cm=dbh,
                height_m=model.predict(dbh) + offset, status=status)


# --- #7: invert range -------------------------------------------------------

def test_invert_reaches_tall_heights():
    model = default_model("")
    dbh = model.invert(32.5)
    assert dbh is not None
    assert abs(model.predict(dbh) - 32.5) < 0.01


def test_invert_above_asymptote_stays_none():
    assert default_model("").invert(35.0) is None


def test_invert_existing_range_unchanged():
    assert abs(default_model("").invert(30.0) - 167.37) < 1.0


# --- #15/#16: residual_scores ------------------------------------------------

def test_residual_scores_exclude_dead_trees():
    model = default_model("")
    trees = [
        on_curve("L1", 20.0),
        on_curve("L2", 22.0),
        on_curve("D1", 20.0, offset=-10.0, status="dead"),
    ]
    scored = residual_scores(model, trees)
    assert [t.tree_id for t, _, _ in scored] == ["L1", "L2"]


def test_residual_scores_loo_breaks_small_group_ceiling():
    model = default_model("")
    offsets = [0.1, -0.1, 0.05, -0.05, 0.0, 0.08, -0.08, 19.34]
    trees = [on_curve(f"T{i}", 20.0 + i, offset=o) for i, o in enumerate(offsets)]
    scored = residual_scores(model, trees)
    z_by_id = {t.tree_id: z for t, _, z in scored}
    assert abs(z_by_id["T7"]) > 2.5


def test_sv030_direction_matches_z_sign():
    trees = [on_curve(f"N{i}", 15.0 + i, offset=3.0) for i in range(9)]
    trees.append(on_curve("LOW", 25.0, offset=0.2))
    model = default_model("")
    ctx = InspectionContext(
        config=InspectionConfig(),
        area_ha=1.0,
        trees=trees,
        metrics=stand_metrics(trees, 1.0),
        metrics_source="field",
        field_trees=trees,
        model=model,
        species_models={"": model},
    )
    findings = [f for f in check_height_anomalies(ctx) if f.code == "SV030"]
    assert len(findings) == 1
    assert findings[0].subject == "tree:LOW"
    assert "shorter" in findings[0].detail


# --- #15/#16: SV031 -----------------------------------------------------------

def test_sv031_fires_for_lone_outlier_in_small_group():
    trees = [on_curve(f"S{i}", 10.0) for i in range(9)]
    trees.append(on_curve("BIG", 200.0))
    report = inspect_stand(trees=trees, area_ha=1.0)
    assert "SV031" in [f.code for f in report.findings]


def test_sv031_ignores_dead_trees():
    trees = [on_curve(f"S{i}", 10.0 + i) for i in range(10)]
    trees.append(Tree("DEAD", 0.0, 0.0, dbh_cm=200.0, height_m=30.0, status="dead"))
    report = inspect_stand(trees=trees, area_ha=1.0)
    assert "SV031" not in [f.code for f in report.findings]


# --- #16: fit_by_species dominant ---------------------------------------------

def test_fit_by_species_dominant_counts_paired_trees():
    trees = []
    for i in range(10):
        trees.append(Tree(f"P{i}", 0.0, 0.0, species="PIAB", dbh_cm=20.0,
                          height_m=20.0, status="dead"))
    for i in range(5):
        trees.append(Tree(f"Q{i}", 0.0, 0.0, species="QURO", dbh_cm=10.0 + i,
                          height_m=12.0 + i))
    models = fit_by_species(trees)
    assert models[""].params == (1.700, 0.195)  # QURO defaults


# --- #17: rugosity --------------------------------------------------------------

def test_rugosity_includes_cells_at_threshold():
    grid = Grid.from_rows([[2.0, 2.0], [2.0, 10.0]])
    assert rugosity(grid, threshold=2.0) == pytest.approx(4.0)


def test_rugosity_default_counts_zero_height_cells():
    grid = Grid.from_rows([[0.0, 0.0], [0.0, 10.0]])
    assert rugosity(grid) == pytest.approx(5.0)


# --- #20: gaps key ----------------------------------------------------------------

def _few_trees():
    return [on_curve(f"T{i}", 15.0 + i) for i in range(6)]


def test_gaps_key_present_when_chm_supplied_and_no_gaps():
    chm = Grid.from_rows([[10.0] * 4] * 4)
    report = inspect_stand(trees=_few_trees(), area_ha=1.0, chm=chm,
                           run_detection=False)
    assert report.as_dict().get("gaps") == []


def test_gaps_key_absent_without_chm():
    report = inspect_stand(trees=_few_trees(), area_ha=1.0)
    assert "gaps" not in report.as_dict()


# --- #24: config validation ---------------------------------------------------------

def test_from_dict_rejects_numeric_strings():
    with pytest.raises(ValueError):
        InspectionConfig.from_dict({"min_recall": "0.8"})


def test_from_dict_accepts_numbers():
    assert InspectionConfig.from_dict({"min_recall": 0.8}).min_recall == 0.8


def test_match_tolerance_must_be_positive():
    with pytest.raises(ValueError):
        InspectionConfig(match_tolerance_m=0)
    with pytest.raises(ValueError):
        InspectionConfig.from_dict({"match_tolerance_m": 0})


# --- #26: allometry input validation --------------------------------------------------

def test_fit_ignores_nonfinite_pairs():
    pairs = [(10.0, 20.0), (12.0, 21.0), (14.0, 22.0), (16.0, 23.0),
             (18.0, 24.0), (20.0, 25.0), (float("inf"), float("inf"))]
    model = fit_height_diameter(pairs, min_points=6)
    assert math.isfinite(model.rmse)


def test_fit_with_only_nonfinite_pairs_raises():
    with pytest.raises(AllometryError):
        fit_height_diameter([(float("inf"), float("inf"))], min_points=1)


def test_curtis_rejects_nonpositive_params():
    curtis = MODELS["curtis"]
    with pytest.raises(ValueError):
        curtis.evaluate(10.0, (-5.0, 3.0))
    with pytest.raises(ValueError):
        curtis.evaluate(10.0, (30.0, -1.0))
    assert curtis.evaluate(10.0, (30.0, 3.0)) > 1.3


# --- #27: gap width stays inside its component -------------------------------------------

def _gap_layout(connected):
    clear = {(r, c) for r in range(2, 7) for c in range(2, 5)}
    if connected:
        clear |= {(r, c) for r in range(2, 7) for c in range(6, 11)}
        clear.add((4, 5))
    rows = []
    for r in range(11):
        rows.append([0.0 if (r, c) in clear else 10.0 for c in range(16)])
    return Grid.from_rows(rows, cellsize=1.0)


def test_gap_width_stays_within_its_component():
    isolated = find_gaps(_gap_layout(False), threshold=2.0, min_area=0.0,
                         min_width=2.0)
    assert len(isolated) == 1
    connected = find_gaps(_gap_layout(True), threshold=2.0, min_area=0.0,
                          min_width=2.0)
    narrow = [g for g in connected if g.area == isolated[0].area]
    assert len(narrow) == 1
    assert narrow[0].inscribed_radius == pytest.approx(
        isolated[0].inscribed_radius, abs=1e-6)

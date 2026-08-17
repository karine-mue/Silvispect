"""Tests for inventory parsing, serialisation and matching."""

from __future__ import annotations

import pytest

from silvispect.detect import Crown, TreeTop
from silvispect.inventory import (
    InventoryError,
    Plot,
    Tree,
    match_trees,
    parse_trees,
    read_trees,
    shift_trees,
    trees_from_crowns,
    trees_to_csv,
    write_trees,
)

CSV = """tree_id,x,y,species,dbh_cm,height_m,status
A1,10.0,20.0,PIAB,30.0,25.0,live
A2,12.5,22.5,FASY,25.0,22.0,live
A3,15.0,25.0,PIAB,,,dead
"""


def crown(tree_id: int, x: float, y: float, height: float, area: float = 12.0) -> Crown:
    top = TreeTop(tree_id, 0, 0, x, y, height)
    return Crown(tree_id, top, ((0, 0),), area, height * 0.7, 2.0)


def test_parse_basic():
    trees = parse_trees(CSV)
    assert len(trees) == 3
    assert trees[0].tree_id == "A1"
    assert trees[0].species == "PIAB"
    assert trees[0].dbh_cm == 30.0
    assert trees[2].dbh_cm is None
    assert trees[2].is_live is False


def test_parse_accepts_aliases_and_odd_case():
    trees = parse_trees("ID,Easting,Northing,DBH,Height,SP\nz1, 1.0 , 2.0 ,44,30,PISY\n")
    tree = trees[0]
    assert (tree.tree_id, tree.x, tree.y) == ("z1", 1.0, 2.0)
    assert tree.dbh_cm == 44.0
    assert tree.height_m == 30.0
    assert tree.species == "PISY"


def test_parse_handles_missing_markers_and_blank_rows():
    trees = parse_trees("x,y,dbh_cm,height_m\n1,1,NA,\n\n2,2,-,NULL\n")
    assert len(trees) == 2
    assert all(tree.dbh_cm is None and tree.height_m is None for tree in trees)


def test_parse_generates_ids_when_absent():
    trees = parse_trees("x,y\n1,1\n2,2\n")
    assert [tree.tree_id for tree in trees] == ["1", "2"]


def test_parse_empty_document():
    assert parse_trees("") == []


def test_parse_requires_coordinates():
    with pytest.raises(InventoryError, match="'x' and 'y'"):
        parse_trees("tree_id,dbh_cm\nA1,30\n")


def test_parse_rejects_non_numeric():
    with pytest.raises(InventoryError, match="not numeric"):
        parse_trees("x,y,dbh_cm\n1,1,thick\n")


def test_parse_rejects_missing_coordinate_value():
    with pytest.raises(InventoryError, match="missing coordinates"):
        parse_trees("x,y\n1,\n")


def test_csv_round_trip(tmp_path):
    trees = parse_trees(CSV)
    path = write_trees(trees, tmp_path / "out" / "trees.csv")
    reloaded = read_trees(path)
    assert [t.as_dict() for t in reloaded] == [t.as_dict() for t in trees]
    assert trees_to_csv(trees).splitlines()[0].startswith("tree_id,x,y,species")


def test_basal_area():
    assert Tree("a", 0, 0, dbh_cm=100.0).basal_area_m2 == pytest.approx(0.7853981, rel=1e-6)
    assert Tree("a", 0, 0, dbh_cm=None).basal_area_m2 is None
    assert Tree("a", 0, 0, dbh_cm=0.0).basal_area_m2 is None


def test_live_status_variants():
    assert Tree("a", 0, 0, status="LIVE").is_live
    assert Tree("a", 0, 0, status="alive").is_live
    assert not Tree("a", 0, 0, status="snag").is_live


def test_plot_requires_positive_area():
    with pytest.raises(InventoryError):
        Plot("p", [], 0.0)
    plot = Plot("p", parse_trees(CSV), 0.1)
    assert len(plot) == 3
    assert len(plot.live_trees()) == 2


def test_trees_from_crowns_without_model():
    trees = trees_from_crowns([crown(1, 5.0, 5.0, 20.0)], species="PIAB")
    assert trees[0].tree_id == "D1"
    assert trees[0].dbh_cm is None
    assert trees[0].height_m == 20.0
    assert trees[0].species == "PIAB"


def test_trees_from_crowns_with_model():
    from silvispect.allometry import default_model

    trees = trees_from_crowns([crown(1, 5.0, 5.0, 25.0)], dbh_model=default_model("PIAB"))
    assert trees[0].dbh_cm is not None
    assert 10.0 < trees[0].dbh_cm < 200.0


def test_matching_pairs_nearest():
    crowns = [crown(1, 10.0, 10.0, 25.0), crown(2, 20.0, 20.0, 20.0)]
    reference = [
        Tree("f1", 10.4, 10.2, height_m=25.5, dbh_cm=30.0),
        Tree("f2", 20.1, 19.9, height_m=20.6, dbh_cm=25.0),
    ]
    match = match_trees(crowns, reference, tolerance=2.5)
    assert match.matched == 2
    assert match.recall == 1.0
    assert match.precision == 1.0
    assert match.f1 == 1.0
    assert match.height_bias == pytest.approx(-0.55, abs=1e-6)
    assert match.height_rmse == pytest.approx(0.5522680508, abs=1e-6)
    assert match.mean_offset < 0.5
    assert match.as_dict()["matched"] == 2


def test_matching_reports_omissions_and_commissions():
    crowns = [crown(1, 10.0, 10.0, 25.0), crown(2, 50.0, 50.0, 18.0)]
    reference = [Tree("f1", 10.2, 10.1, height_m=25.0), Tree("f2", 80.0, 80.0)]
    match = match_trees(crowns, reference, tolerance=2.0)
    assert match.matched == 1
    assert [tree.tree_id for tree in match.omissions] == ["f2"]
    assert [c.tree_id for c in match.commissions] == [2]
    assert match.recall == 0.5
    assert match.precision == 0.5
    assert match.f1 == pytest.approx(0.5)


def test_matching_is_one_to_one():
    crowns = [crown(1, 10.0, 10.0, 25.0), crown(2, 10.5, 10.0, 24.0)]
    reference = [Tree("f1", 10.1, 10.0, height_m=25.0)]
    match = match_trees(crowns, reference, tolerance=3.0)
    assert match.matched == 1
    assert len(match.commissions) == 1
    assert match.matches[0][0].tree_id == 1


def test_matching_skips_dead_trees_by_default():
    crowns = [crown(1, 10.0, 10.0, 25.0)]
    reference = [Tree("f1", 10.0, 10.0, height_m=25.0, status="dead")]
    assert match_trees(crowns, reference).matched == 0
    assert match_trees(crowns, reference, live_only=False).matched == 1


def test_matching_with_no_data_is_empty():
    match = match_trees([], [])
    assert match.recall == 0.0
    assert match.precision == 0.0
    assert match.f1 == 0.0
    assert match.height_bias is None
    assert match.height_rmse is None
    assert match.mean_offset is None


def test_matching_rejects_bad_tolerance():
    with pytest.raises(InventoryError):
        match_trees([], [], tolerance=0.0)


def test_shift_trees():
    moved = shift_trees([Tree("a", 1.0, 2.0)], 10.0, -1.0)
    assert (moved[0].x, moved[0].y) == (11.0, 1.0)

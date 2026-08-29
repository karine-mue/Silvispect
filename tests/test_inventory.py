"""Tests for inventory parsing, serialisation and matching."""

from __future__ import annotations

import math

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


def test_matching_is_independent_of_record_order():
    """Equal distances must not let CSV row order decide the match count."""
    crowns = [crown(1, 0.0, 0.0, 10.0), crown(2, 2.0, 0.0, 20.0)]
    a = Tree("A", 1.0, 0.0, height_m=11.0)
    b = Tree("B", 0.0, 1.0, height_m=31.0)
    forward = match_trees(crowns, [a, b], tolerance=1.01)
    backward = match_trees(crowns, [b, a], tolerance=1.01)
    assert forward.as_dict() == backward.as_dict()


def test_matching_is_permutation_invariant_over_many_shuffles():
    import random

    rng = random.Random(4)
    crowns = [crown(i, float(i % 5), float(i // 5), 20.0) for i in range(1, 16)]
    trees = [Tree(f"f{i}", float(i % 5) + 0.2, float(i // 5), height_m=20.0) for i in range(15)]
    baseline = match_trees(crowns, trees).as_dict()
    for _ in range(25):
        shuffled = trees[:]
        rng.shuffle(shuffled)
        assert match_trees(crowns, shuffled).as_dict() == baseline


def _rigid_transforms():
    """Isometries of the plane that a stand's geometry must be indifferent to."""
    return {
        "identity": lambda x, y: (x, y),
        "rotate_180": lambda x, y: (-x, -y),
        "rotate_90": lambda x, y: (-y, x),
        "reflect_x": lambda x, y: (-x, y),
        "transpose": lambda x, y: (y, x),
        "translate": lambda x, y: (x + 137.5, y - 4021.25),
    }


def test_matching_is_invariant_under_rigid_motion():
    """Matching describes a forest, not a coordinate frame.

    Breaking equal-distance ties on position merely traded a dependence on CSV
    row order for a dependence on which way the plot happened to be oriented:
    the same two crowns and two stems, all exactly one metre apart, matched one
    pair in one frame and two after rotating the whole scene.  How many pairs a
    tie yields is a property of the geometry, so every isometry must agree.
    """
    import random

    rng = random.Random(2026)
    for trial in range(60):
        count = rng.randint(2, 7)
        # Place stems on a lattice and crowns exactly one metre away in a random
        # direction, which manufactures ties in bulk.
        stems = [(float(rng.randint(0, 4)), float(rng.randint(0, 4))) for _ in range(count)]
        tops = [
            (x + rng.choice([-1.0, 1.0, 0.0]), y + rng.choice([0.0, -1.0, 1.0])) for x, y in stems
        ]
        baseline = None
        for name, move in _rigid_transforms().items():
            crowns = [crown(i + 1, *move(x, y), 20.0 + i) for i, (x, y) in enumerate(tops)]
            trees = [
                Tree(f"S{i:02d}", *move(x, y), height_m=20.0 + i) for i, (x, y) in enumerate(stems)
            ]
            result = match_trees(crowns, trees, tolerance=1.01)
            summary = (len(result.matches), len(result.omissions), len(result.commissions))
            if baseline is None:
                baseline = summary
            assert summary == baseline, f"trial {trial}: {name} disagreed with identity"


def test_matching_never_gives_up_a_nearer_pair():
    """Resolving ties must not reroute a pair that no tie was involved in."""
    crowns = [crown(1, 0.0, 0.0, 20.0), crown(2, 1.0, 0.0, 20.0)]
    trees = [Tree("A", 0.1, 0.0, height_m=20.0), Tree("B", 2.0, 0.0, height_m=20.0)]
    result = match_trees(crowns, trees, tolerance=1.5)
    paired = {c.tree_id: t.tree_id for c, t, _ in result.matches}
    assert paired == {1: "A", 2: "B"}


def test_matching_pairs_as_many_ties_as_the_geometry_allows():
    """Two crowns and two stems mutually one metre apart admit two pairs."""
    crowns = [crown(1, 0.0, 0.0, 10.0), crown(2, 2.0, 0.0, 20.0)]
    trees = [Tree("A", 1.0, 0.0, height_m=11.0), Tree("B", 0.0, 1.0, height_m=31.0)]
    assert len(match_trees(crowns, trees, tolerance=1.01).matches) == 2


def test_matching_count_does_not_depend_on_who_the_stems_are():
    """Relabelling two stems cannot change how many crowns get matched.

    Two stems sit a metre either side of one crown; a second crown is two
    metres from the left-hand stem and out of reach of the right-hand one.
    The nearest distance can be honoured either way, but only one of those ways
    leaves a partner for the second crown.  Settling each distance on its own
    and never looking further picked by identifier, so swapping the two labels
    — the same forest, the same geometry — moved the match count from one to
    two.
    """
    crowns = [crown(1, 0.0, 0.0, 20.0), crown(2, -3.0, 0.0, 20.0)]
    counts = set()
    for left, right in (("A", "B"), ("B", "A")):
        trees = [Tree(left, -1.0, 0.0, height_m=20.0), Tree(right, 1.0, 0.0, height_m=20.0)]
        result = match_trees(crowns, trees, tolerance=2.5)
        counts.add(len(result.matches))
        assert sorted(round(d, 6) for _, _, d in result.matches) == [1.0, 2.0]
    assert counts == {2}


def _best_distance_profile(crowns, trees, tolerance):
    """Exhaustive optimum: the largest per-distance count vector, lexicographic."""
    import itertools

    edges = [
        (i, j, math.hypot(c.x - t.x, c.y - t.y))
        for i, c in enumerate(crowns)
        for j, t in enumerate(trees)
        if math.hypot(c.x - t.x, c.y - t.y) <= tolerance
    ]
    tiers = sorted({distance for _, _, distance in edges})
    rank = {distance: index for index, distance in enumerate(tiers)}
    best = tuple([0] * len(tiers))
    for size in range(len(crowns) + 1):
        for combination in itertools.combinations(range(len(edges)), size):
            used_crowns: set[int] = set()
            used_trees: set[int] = set()
            profile = [0] * len(tiers)
            for edge in combination:
                i, j, distance = edges[edge]
                if i in used_crowns or j in used_trees:
                    break
                used_crowns.add(i)
                used_trees.add(j)
                profile[rank[distance]] += 1
            else:
                best = max(best, tuple(profile))
    return best


def _profile(result, tolerance, crowns, trees):
    tiers = sorted(
        {
            math.hypot(c.x - t.x, c.y - t.y)
            for c in crowns
            for t in trees
            if math.hypot(c.x - t.x, c.y - t.y) <= tolerance
        }
    )
    rank = {distance: index for index, distance in enumerate(tiers)}
    profile = [0] * len(tiers)
    for _, _, distance in result.matches:
        profile[rank[distance]] += 1
    return tuple(profile)


def test_matching_reaches_the_exhaustive_optimum():
    """The stated objective, checked against every possible pairing.

    Honour as many pairs as possible at the shortest distance, then as many as
    possible at the next without giving any of those up, and so on.  Checked
    here against brute-force enumeration of every partial pairing, because a
    rule about ties is only worth as much as the cases nobody thought of.
    """
    import random

    rng = random.Random(20260819)
    for _ in range(300):
        count_c, count_t = rng.randint(1, 4), rng.randint(1, 4)
        crowns = [
            crown(i + 1, float(rng.randint(-3, 3)), float(rng.randint(-3, 3)), 20.0)
            for i in range(count_c)
        ]
        trees = [
            Tree(f"S{j}", float(rng.randint(-3, 3)), float(rng.randint(-3, 3)), height_m=20.0)
            for j in range(count_t)
        ]
        tolerance = rng.choice([1.0, 2.0, 2.5, 3.0])
        result = match_trees(crowns, trees, tolerance=tolerance)
        assert _profile(result, tolerance, crowns, trees) == _best_distance_profile(
            crowns, trees, tolerance
        )


def test_matching_count_survives_relabelling_and_reordering():
    """Neither the names on the stems nor the order of the rows may matter."""
    import random

    rng = random.Random(11)
    for _ in range(120):
        crowns = [
            crown(i + 1, float(rng.randint(-4, 4)), float(rng.randint(-4, 4)), 20.0 + i * 0.1)
            for i in range(rng.randint(2, 5))
        ]
        places = [(float(rng.randint(-4, 4)), float(rng.randint(-4, 4))) for _ in range(4)]
        heights = [20.0 + rng.random() for _ in places]
        names = [f"S{j}" for j in range(len(places))]

        def build(labels, order):
            return [
                Tree(labels[j], *places[j], height_m=heights[j])  # noqa: B023
                for j in order
            ]

        def summary(trees, crowns=crowns):
            outcome = match_trees(crowns, trees, tolerance=2.5)
            return (len(outcome.matches), len(outcome.omissions), len(outcome.commissions))

        baseline = summary(build(names, list(range(len(places)))))

        shuffled_rows = list(range(len(places)))
        rng.shuffle(shuffled_rows)
        assert summary(build(names, shuffled_rows)) == baseline

        relabelled = names[:]
        rng.shuffle(relabelled)
        assert summary(build(relabelled, list(range(len(places))))) == baseline


def test_matching_never_trades_a_nearer_pair_for_more_farther_ones():
    """One pair at half a metre outranks two at two metres."""
    crowns = [crown(1, 0.0, 0.0, 20.0)]
    trees = [
        Tree("near", 0.5, 0.0, height_m=20.0),
        Tree("far", 2.0, 0.0, height_m=20.0),
    ]
    result = match_trees(crowns, trees, tolerance=2.5)
    assert [t.tree_id for _, t, _ in result.matches] == ["near"]

    # Two crowns competing for one near stem: the nearer claim wins and the
    # other crown falls back to its own second-choice stem rather than
    # displacing it.
    crowns = [crown(1, 0.0, 0.0, 20.0), crown(2, 1.0, 0.0, 20.0)]
    trees = [Tree("a", 0.1, 0.0, height_m=20.0), Tree("b", 2.0, 0.0, height_m=20.0)]
    paired = {c.tree_id: t.tree_id for c, t, _ in match_trees(crowns, trees, tolerance=2.5).matches}
    assert paired == {1: "a", 2: "b"}


def test_matching_resolves_independent_tie_clusters_together():
    """Several tied clusters in one plot are each settled on their own merits."""
    crowns = []
    trees = []
    for cluster in range(3):
        offset = cluster * 100.0
        crowns.append(crown(2 * cluster + 1, offset + 0.0, 0.0, 20.0))
        crowns.append(crown(2 * cluster + 2, offset - 3.0, 0.0, 20.0))
        trees.append(Tree(f"L{cluster}", offset - 1.0, 0.0, height_m=20.0))
        trees.append(Tree(f"R{cluster}", offset + 1.0, 0.0, height_m=20.0))
    result = match_trees(crowns, trees, tolerance=2.5)
    assert len(result.matches) == 6
    assert sorted(round(d, 6) for _, _, d in result.matches) == [1.0, 1.0, 1.0, 2.0, 2.0, 2.0]


def test_inventory_csv_is_a_quantised_record():
    """CSV output rounds to ``precision`` decimals, and that is the contract.

    An inventory is a record of measurements, so it is written at the
    resolution a field instrument reports rather than at the resolution a float
    happens to hold.  The consequence is real and documented in
    ``docs/data-formats.md``: a stem four ten-thousandths outside a tolerance
    is inside it after a round trip.  Raising ``precision`` is the way out, and
    it has to keep working.
    """
    tree = Tree("A", 2.5004, 0.0, height_m=20.0)
    assert parse_trees(trees_to_csv([tree]))[0].x == 2.5
    assert parse_trees(trees_to_csv([tree], precision=6))[0].x == 2.5004

    crowns = [crown(1, 0.0, 0.0, 20.0)]
    assert len(match_trees(crowns, [tree], tolerance=2.5).matches) == 0
    rounded = parse_trees(trees_to_csv([tree]))
    assert len(match_trees(crowns, rounded, tolerance=2.5).matches) == 1
    kept = parse_trees(trees_to_csv([tree], precision=6))
    assert len(match_trees(crowns, kept, tolerance=2.5).matches) == 0

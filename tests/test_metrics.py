"""Tests for stand mensuration metrics."""

from __future__ import annotations

import decimal
import math
import sys
from decimal import Decimal

import pytest

from silvispect.inventory import MAX_DBH_CM, Tree
from silvispect.metrics import (
    DEFAULT_WOOD_DENSITY,
    MIN_CLASS_WIDTH_CM,
    MetricsError,
    above_ground_biomass,
    basal_area,
    diameter_distribution,
    dominant_height,
    gini_coefficient,
    lorey_height,
    quadratic_mean_diameter,
    reineke_sdi,
    shannon_index,
    simpson_index,
    species_composition,
    stand_metrics,
    stem_volume,
)


def tree(tree_id, dbh, height, species="PIAB", status="live"):
    return Tree(tree_id, 0.0, 0.0, species=species, dbh_cm=dbh, height_m=height, status=status)


def test_basal_area_of_a_known_diameter():
    # A 20 cm stem has a cross-section of pi * 0.1^2 m2.
    assert basal_area([tree("a", 20.0, 15.0)]) == pytest.approx(math.pi * 0.01)


def test_basal_area_ignores_unmeasured():
    assert basal_area([tree("a", None, 15.0)]) == 0.0


def test_quadratic_mean_diameter():
    # QMD is the root mean square of the diameters, so it exceeds the mean.
    trees = [tree("a", 10.0, 10.0), tree("b", 30.0, 20.0)]
    assert quadratic_mean_diameter(trees) == pytest.approx(math.sqrt(500.0))
    assert quadratic_mean_diameter([]) is None


def test_lorey_height_is_basal_area_weighted():
    trees = [tree("a", 10.0, 10.0), tree("b", 30.0, 30.0)]
    weights = [math.pi * (10 / 200) ** 2, math.pi * (30 / 200) ** 2]
    expected = (weights[0] * 10 + weights[1] * 30) / sum(weights)
    assert lorey_height(trees) == pytest.approx(expected)
    assert lorey_height(trees) > 25.0  # dominated by the thick stem
    assert lorey_height([tree("a", None, 10.0)]) is None


def test_dominant_height_takes_the_thickest():
    trees = [tree(f"t{i}", float(i + 10), float(i + 5)) for i in range(50)]
    # 100 stems/ha on 0.1 ha means the 10 thickest stems.
    assert dominant_height(trees, 0.1) == pytest.approx(sum(range(45, 55)) / 10.0)
    assert dominant_height([], 1.0) is None
    with pytest.raises(MetricsError):
        dominant_height(trees, 0.0)


def test_reineke_sdi_at_the_reference_diameter():
    assert reineke_sdi(500.0, 25.0) == pytest.approx(500.0)
    assert reineke_sdi(500.0, 50.0) > 500.0
    with pytest.raises(MetricsError):
        reineke_sdi(500.0, 0.0)


def test_gini_bounds():
    assert gini_coefficient([5.0, 5.0, 5.0]) == pytest.approx(0.0)
    assert gini_coefficient([0.0, 0.0, 0.0]) == 0.0
    assert gini_coefficient([1.0]) is None
    assert 0.0 < gini_coefficient([1.0, 2.0, 3.0, 40.0]) < 1.0


def test_species_diversity():
    trees = [tree("a", 20.0, 15.0, "PIAB"), tree("b", 20.0, 15.0, "FASY")]
    assert species_composition(trees) == {"PIAB": 0.5, "FASY": 0.5}
    assert shannon_index(trees) == pytest.approx(math.log(2))
    assert simpson_index(trees) == pytest.approx(0.5)

    pure = [tree("a", 20.0, 15.0, "PIAB")] * 3
    assert shannon_index(pure) == pytest.approx(0.0)
    assert simpson_index(pure) == pytest.approx(0.0)
    assert shannon_index([]) is None
    assert species_composition([]) == {}


def test_species_composition_labels_blanks():
    assert species_composition([tree("a", 20.0, 15.0, "")]) == {"UNKNOWN": 1.0}


def test_diameter_distribution():
    trees = [tree("a", 7.0, 5.0), tree("b", 12.0, 8.0), tree("c", 13.0, 9.0)]
    assert diameter_distribution(trees, class_width=5.0) == {"5-10": 1, "10-15": 2}
    with pytest.raises(MetricsError):
        diameter_distribution(trees, class_width=0.0)


def test_stem_volume_and_biomass():
    subject = tree("a", 40.0, 25.0)
    assert stem_volume(subject) == pytest.approx(0.5 * math.pi * 0.2**2 * 25.0)
    assert stem_volume(subject, form_factor=1.0) == pytest.approx(math.pi * 0.2**2 * 25.0)
    assert stem_volume(tree("a", None, 25.0)) is None
    assert stem_volume(tree("a", 40.0, None)) is None

    biomass = above_ground_biomass(subject)
    assert biomass is not None and 500.0 < biomass < 3000.0
    denser = above_ground_biomass(subject, wood_density=0.8)
    assert denser > biomass
    assert above_ground_biomass(tree("a", None, 25.0)) is None


def test_biomass_uses_the_species_density():
    spruce = above_ground_biomass(tree("a", 30.0, 20.0, "PIAB"))
    beech = above_ground_biomass(tree("a", 30.0, 20.0, "FASY"))
    assert beech > spruce  # beech wood is denser


def test_stand_metrics_end_to_end():
    trees = [
        tree("a", 20.0, 18.0, "PIAB"),
        tree("b", 30.0, 24.0, "PIAB"),
        tree("c", 40.0, 28.0, "FASY"),
        tree("d", None, None, "FASY"),
        tree("e", 25.0, 20.0, "PIAB", status="dead"),
    ]
    metrics = stand_metrics(trees, area_ha=0.25)
    assert metrics.tree_count == 4  # the dead stem is excluded
    assert metrics.measured_count == 3
    assert metrics.stems_per_ha == pytest.approx(16.0)
    assert metrics.basal_area_per_ha == pytest.approx(basal_area(trees[:3]) / 0.25)
    assert metrics.quadratic_mean_diameter_cm == pytest.approx(math.sqrt((400 + 900 + 1600) / 3))
    assert metrics.max_height_m == 28.0
    assert metrics.sdi > 0
    assert metrics.species_shares["PIAB"] == pytest.approx(0.5)
    assert sum(metrics.diameter_classes.values()) == 3
    payload = metrics.as_dict()
    assert payload["tree_count"] == 4
    assert payload["area_ha"] == 0.25


def test_stand_metrics_can_include_dead():
    trees = [tree("a", 20.0, 18.0), tree("b", 25.0, 20.0, status="dead")]
    assert stand_metrics(trees, 1.0, live_only=False).tree_count == 2


def test_stand_metrics_without_measurements():
    metrics = stand_metrics([Tree("a", 0.0, 0.0)], area_ha=1.0)
    assert metrics.measured_count == 0
    assert metrics.basal_area_per_ha is None
    assert metrics.quadratic_mean_diameter_cm is None
    assert metrics.sdi is None
    assert metrics.volume_per_ha_m3 is None


def test_stand_metrics_requires_positive_area():
    with pytest.raises(MetricsError):
        stand_metrics([], area_ha=0.0)


def test_metrics_of_a_synthetic_stand_are_realistic(stand):
    metrics = stand_metrics(stand.trees, stand.area_ha)
    assert 100 < metrics.stems_per_ha < 600
    assert 5 < metrics.basal_area_per_ha < 80
    assert 10 < metrics.lorey_height_m < 45
    assert metrics.dominant_height_m >= metrics.mean_height_m
    assert 0 < metrics.gini_basal_area < 1
    assert metrics.volume_per_ha_m3 > 0


def test_yield_totals_are_unknown_not_zero_without_heights():
    """A stand with diameters but no heights has an unknown volume, not zero."""
    metrics = stand_metrics([tree("a", 30.0, None)], area_ha=1.0)
    assert metrics.measured_count == 1
    assert metrics.yield_basis_count == 0
    assert metrics.volume_per_ha_m3 is None
    assert metrics.biomass_per_ha_t is None


def test_yield_basis_count_exposes_the_partial_subset():
    """Totals cover only the stems with both dimensions, and say so."""
    trees = [tree("a", 30.0, 22.0), tree("b", 30.0, None)]
    metrics = stand_metrics(trees, area_ha=1.0)
    assert metrics.measured_count == 2
    assert metrics.yield_basis_count == 1
    assert metrics.volume_per_ha_m3 == pytest.approx(stem_volume(trees[0]))
    assert metrics.as_dict()["yield_basis_count"] == 1


def test_dominant_height_tie_at_the_cutoff_ignores_record_order():
    """Two equal diameters at the cutoff must not resolve by row order."""
    short = Tree("a", 0.0, 0.0, dbh_cm=30.0, height_m=10.0)
    tall = Tree("b", 1.0, 1.0, dbh_cm=30.0, height_m=30.0)
    assert dominant_height([short, tall], 0.01) == dominant_height([tall, short], 0.01)


def test_lorey_height_does_not_depend_on_row_order():
    """The same stand summarised in any order is the same stand.

    A running total drops the small stems once a large one has been added, and
    how much it drops depends on when the large one arrived.  Summed exactly,
    reversing the list cannot move the answer.
    """
    import random

    trees = [Tree("t0", 0.0, 0.0, dbh_cm=30.0, height_m=1e16)] + [
        Tree(f"t{i}", 0.0, 0.0, dbh_cm=30.0, height_m=1.0) for i in range(1, 21)
    ]
    baseline = lorey_height(trees)
    assert lorey_height(list(reversed(trees))) == baseline
    rng = random.Random(9)
    for _ in range(30):
        shuffled = trees[:]
        rng.shuffle(shuffled)
        assert lorey_height(shuffled) == baseline

    ordinary = [
        Tree(f"n{i}", 0.0, 0.0, dbh_cm=20.0 + i, height_m=15.0 + i * 0.5) for i in range(25)
    ]
    assert lorey_height(list(reversed(ordinary))) == lorey_height(ordinary)


def test_a_diameter_class_is_a_band_not_a_rounded_centimetre():
    """A class narrower than a centimetre still has to be its own class.

    The lower bound was truncated to a whole number, so at a width of 0.5 cm
    two stems a whole band apart were counted together under a label that read
    ``0-0``.  The band number carries the width; the boundaries are computed
    back from it.
    """
    trees = [
        Tree(tree_id=str(i), x=0.0, y=0.0, dbh_cm=dbh) for i, dbh in enumerate([0.2, 0.7, 1.1])
    ]
    assert diameter_distribution(trees, class_width=0.5) == {"0-0.5": 1, "0.5-1": 1, "1-1.5": 1}

    # A whole-centimetre width still reads exactly as it always did.
    ordinary = [
        Tree(tree_id=str(i), x=0.0, y=0.0, dbh_cm=dbh)
        for i, dbh in enumerate([3.0, 12.0, 17.5, 44.0])
    ]
    assert diameter_distribution(ordinary) == {"0-5": 1, "10-15": 1, "15-20": 1, "40-45": 1}

    # Every stem lands in exactly one class, whatever the width.
    for width in (0.25, 0.5, 1.0, 2.5, 5.0, 7.5):
        classes = diameter_distribution(ordinary, class_width=width)
        assert sum(classes.values()) == len(ordinary)
        assert len(set(classes)) == len(classes)


def test_metrics_stay_finite_for_stems_the_inventory_accepts():
    """Nothing an inventory will hold may overflow the summaries that read it.

    Squaring a diameter before scaling it reaches the finite limit above about
    1.3e154 cm, and a plain total reaches it on a running sum whose answer is
    representable.  Both are far inside the range a CSV is allowed to state.
    """
    big = 1e150
    trees = [Tree(tree_id=str(i), x=float(i), y=0.0, dbh_cm=big, height_m=30.0) for i in range(4)]
    assert math.isfinite(quadratic_mean_diameter(trees))
    assert quadratic_mean_diameter(trees) == pytest.approx(big)
    assert math.isfinite(basal_area(trees))
    assert math.isfinite(lorey_height(trees))
    assert lorey_height(trees) == pytest.approx(30.0)
    assert gini_coefficient([tree.basal_area_m2 for tree in trees]) == pytest.approx(0.0)

    # Reineke's index really is out of range for these, and says so rather
    # than handing an infinity to a JSON writer.
    assert math.isfinite(reineke_sdi(1000.0, big))
    with pytest.raises(MetricsError, match="not representable"):
        reineke_sdi(1000.0, 1e200)


def test_lorey_height_still_ignores_row_order_after_scaling():
    """The exact-summation contract from the previous repair must survive."""
    trees = [
        Tree(tree_id="tall", x=0.0, y=0.0, dbh_cm=100.0, height_m=1e16),
        *(Tree(tree_id=str(i), x=float(i), y=0.0, dbh_cm=10.0, height_m=1.0) for i in range(20)),
    ]
    assert lorey_height(trees) == lorey_height(list(reversed(trees)))


def test_lorey_height_survives_heights_as_large_as_diameters():
    """Normalising one side of a weighted mean is not enough.

    Two equal weights become 1 and 1, and the heights then reach the finite
    limit in the sum instead: two stems of the largest finite height have
    exactly that height as their weighted mean, and asking for it raised.
    """
    peak = sys.float_info.max
    equal = [
        Tree(tree_id="a", x=0.0, y=0.0, dbh_cm=1.0, height_m=peak),
        Tree(tree_id="b", x=1.0, y=0.0, dbh_cm=1.0, height_m=peak),
    ]
    # To the last ulp a weighted mean is a division, not an identity; what it
    # must not do is leave the range.
    assert lorey_height(equal) == pytest.approx(peak, rel=1e-15)

    mixed = [
        Tree(tree_id="a", x=0.0, y=0.0, dbh_cm=1e150, height_m=peak / 2),
        Tree(tree_id="b", x=1.0, y=0.0, dbh_cm=1e-3, height_m=1.0),
    ]
    assert math.isfinite(lorey_height(mixed))

    # Ordinary stands are unaffected.
    plain = [
        Tree(tree_id="a", x=0.0, y=0.0, dbh_cm=20.0, height_m=15.0),
        Tree(tree_id="b", x=1.0, y=0.0, dbh_cm=40.0, height_m=25.0),
    ]
    assert lorey_height(plain) == pytest.approx(23.0, abs=0.01)


def test_biomass_does_not_overflow_on_a_diameter_it_accepts():
    """``(rho D^2 H)^0.976`` must not be computed by squaring first.

    The square reaches the finite limit above about 1.3e154 cm even when the
    answer, once the exponent is applied, is an ordinary finite number.
    """
    value = above_ground_biomass(Tree(tree_id="x", x=0.0, y=0.0, dbh_cm=1e155, height_m=1.0))
    assert math.isfinite(value)
    assert value == pytest.approx(1.2422522398e301, rel=1e-9)

    # The distributed form is the same equation, so ordinary stems are unmoved.
    ordinary = Tree(tree_id="t", x=0.0, y=0.0, dbh_cm=30.0, height_m=20.0, species="FASY")
    rho = 0.58
    assert above_ground_biomass(ordinary) == pytest.approx(
        0.0673 * (rho * 30.0**2 * 20.0) ** 0.976, rel=1e-12
    )


@pytest.mark.parametrize("class_width", [0.1, 0.25, 0.5, 1.0, 2.5, 5.0, 10.0])
def test_every_stem_lands_in_the_class_that_contains_it(class_width):
    """Checked against an independent decimal interval oracle.

    Binary floats put ``0.3 / 0.1`` a hair below three, so a stem sitting
    exactly on a class boundary fell into the band below and dragged a label
    like ``0.2-0.30000000000000004`` with it.  The oracle here reads the
    printed bounds back and asserts ``lower <= dbh < upper`` in decimal, which
    is the claim the label makes.
    """
    diameters = [0.1, 0.3, 0.5, 1.0, 2.5, 5.0, 7.5, 12.3, 30.0, 44.0]
    trees = [
        Tree(tree_id=str(index), x=0.0, y=0.0, dbh_cm=dbh) for index, dbh in enumerate(diameters)
    ]
    classes = diameter_distribution(trees, class_width=class_width)
    assert sum(classes.values()) == len(diameters)

    placed = 0
    for label, count in classes.items():
        low_text, high_text = label.split("-")
        low, high = Decimal(low_text), Decimal(high_text)
        assert high - low == Decimal(repr(class_width))
        inside = [d for d in diameters if low <= Decimal(repr(d)) < high]
        assert len(inside) == count, (label, inside)
        placed += count
    assert placed == len(diameters)


def test_a_stem_on_a_class_boundary_starts_the_upper_class():
    """The reported case, and the label it used to carry."""
    stem = [Tree(tree_id="t", x=0.0, y=0.0, dbh_cm=0.3)]
    assert diameter_distribution(stem, class_width=0.1) == {"0.3-0.4": 1}
    assert diameter_distribution(stem, class_width=0.3) == {"0.3-0.6": 1}


def _lorey_oracle(trees):
    """The weighted mean at 400 digits, far outside anything a float can lose."""
    with decimal.localcontext() as context:
        context.prec = 400
        pairs = [(t.basal_area_m2, t.height_m) for t in trees]
        top = sum(decimal.Decimal(w) * decimal.Decimal(h) for w, h in pairs)
        return top / sum(decimal.Decimal(w) for w, _ in pairs)


def _biomass_oracle(dbh_cm, height_m, wood_density=DEFAULT_WOOD_DENSITY):
    with decimal.localcontext() as context:
        context.prec = 120
        inner = (
            decimal.Decimal(wood_density) * decimal.Decimal(dbh_cm) ** 2 * decimal.Decimal(height_m)
        )
        return decimal.Decimal("0.0673") * (inner.ln() * decimal.Decimal("0.976")).exp()


@pytest.mark.parametrize(
    "first,second",
    [
        ((1e-159, sys.float_info.max), (1e156, 1e-308)),
        ((1.0, sys.float_info.max), (1.0, sys.float_info.max)),
        ((1e150, 1e-300), (1e-150, 1e300)),
        ((1e-160, 1e-320), (1e155, 1e-100)),
        ((20.0, 15.0), (40.0, 25.0)),
        ((1e-3, 1e-3), (1e150, 1e150)),
    ],
)
def test_lorey_height_matches_a_high_precision_oracle_across_the_range(first, second):
    """Scaling by the largest value protects one end and destroys the other.

    Normalising a subnormal basal area against one of 1e307, or a 1e-308 height
    against the largest float, sends the term to zero — so a pair whose
    weighted mean is 1e-308 came back as 0, having been repaired out of an
    overflow into an underflow.  The exponents are carried as integers now, so
    neither end can be lost.
    """
    trees = [
        Tree(tree_id="a", x=0.0, y=0.0, dbh_cm=first[0], height_m=first[1]),
        Tree(tree_id="b", x=1.0, y=0.0, dbh_cm=second[0], height_m=second[1]),
    ]
    value = lorey_height(trees)
    expected = float(_lorey_oracle(trees))
    assert value is not None and math.isfinite(value) and value > 0.0
    assert value == pytest.approx(expected, rel=1e-12)
    # The exact-summation contract from the earlier repair still holds.
    assert lorey_height(trees) == lorey_height(list(reversed(trees)))


@pytest.mark.parametrize(
    "dbh_cm,height_m",
    [
        (1e154, 5e-324),
        (1e155, 1.0),
        (1e156, 1e-320),
        (1e-100, 1e-100),
        (1.0, 1e300),
        (30.0, 20.0),
        (0.5, 0.5),
    ],
)
def test_biomass_matches_a_high_precision_oracle_across_the_range(dbh_cm, height_m):
    """Distributing the exponent traded an overflow for an underflow.

    ``(rho * H) ** 0.976`` sends a height of 5e-324 m to zero and the whole
    biomass with it, for a result an ordinary float holds without difficulty.
    Neither end is lost when the bracket is carried as a mantissa and a power
    of two.
    """
    value = above_ground_biomass(Tree(tree_id="x", x=0.0, y=0.0, dbh_cm=dbh_cm, height_m=height_m))
    expected = float(_biomass_oracle(dbh_cm, height_m))
    assert value is not None and math.isfinite(value) and value > 0.0
    assert value == pytest.approx(expected, rel=1e-12)


def test_a_diameter_class_narrower_than_the_record_is_refused():
    """The inventory keeps three decimals, so a finer band separates nothing.

    Below that the arithmetic could not answer either: at 3e-28 cm the decimal
    division ran out of digits and produced a class whose own interval excluded
    the stem inside it, and a shade narrower ``decimal.InvalidOperation``
    escaped from the library.
    """
    stem = [Tree(tree_id="t", x=0.0, y=0.0, dbh_cm=1.0)]
    for bad in (3e-29, 3e-28, MIN_CLASS_WIDTH_CM * 0.999, 0.0, -1.0):
        with pytest.raises(MetricsError, match="class_width"):
            diameter_distribution(stem, class_width=bad)

    # Exactly at the bound the answer is ordinary, and it contains the stem.
    at_bound = diameter_distribution(stem, class_width=MIN_CLASS_WIDTH_CM)
    (label,) = at_bound
    low, high = (Decimal(part) for part in label.split("-"))
    assert low <= Decimal("1.0") < high
    assert high - low == Decimal(repr(MIN_CLASS_WIDTH_CM))

    # The widest stem the inventory accepts, at the narrowest band it accepts.
    widest = [Tree(tree_id="w", x=0.0, y=0.0, dbh_cm=MAX_DBH_CM)]
    (extreme,) = diameter_distribution(widest, class_width=MIN_CLASS_WIDTH_CM)
    low, high = (Decimal(part) for part in extreme.split("-"))
    assert low <= Decimal(repr(MAX_DBH_CM)) < high


def _lorey_from_diameters(pairs):
    """The oracle from the definition, with the square taken in decimal."""
    with decimal.localcontext() as context:
        context.prec = 400
        pi = decimal.Decimal(math.pi)
        top = sum(pi * (decimal.Decimal(d) / 200) ** 2 * decimal.Decimal(h) for d, h in pairs)
        bottom = sum(pi * (decimal.Decimal(d) / 200) ** 2 for d, _ in pairs)
        return top / bottom


@pytest.mark.parametrize(
    "pairs",
    [
        [(2e-160, sys.float_info.max), (1.0, 1e-20)],
        [(1e-159, sys.float_info.max), (1e156, 1e-308)],
        [(2e-160, 1.0), (2e-160, 3.0)],
        [(3e-160, 1e-300), (1.0, 1e-320)],
        [(1e-170, 25.0), (30.0, 25.0), (1e150, 25.0)],
        [(1e155, sys.float_info.max), (1e155, sys.float_info.max)],
        [(1e-165, 1e300), (1e-165, 1e-300)],
        [(20.0, 15.0), (40.0, 25.0), (12.5, 9.0)],
    ],
)
def test_lorey_height_never_forms_the_basal_area_as_a_float(pairs):
    """The weight of a stem is carried as digits and an exponent, never a float.

    Squaring a diameter of 2e-160 cm underflows to zero before any careful
    summation can see it — and because that stem's enormous height still set
    the scale of the numerator, the other stem's real contribution was shifted
    out of range too, and a mean of 7e-12 came back as 0.  Every case here has
    a positive, representable answer, and the oracle takes the square in
    decimal at 400 digits.
    """
    trees = [
        Tree(tree_id=str(index), x=float(index), y=0.0, dbh_cm=d, height_m=h)
        for index, (d, h) in enumerate(pairs)
    ]
    value = lorey_height(trees)
    expected = float(_lorey_from_diameters(pairs))
    assert value is not None and math.isfinite(value) and value > 0.0
    assert value == pytest.approx(expected, rel=1e-12)
    assert lorey_height(list(reversed(trees))) == value


def test_lorey_height_reported_case():
    trees = [
        Tree(tree_id="a", x=0.0, y=0.0, dbh_cm=2e-160, height_m=sys.float_info.max),
        Tree(tree_id="b", x=1.0, y=0.0, dbh_cm=1.0, height_m=1e-20),
    ]
    assert lorey_height(trees) == pytest.approx(7.1907725494492624e-12, rel=1e-12)

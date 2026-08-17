"""Tests for stand mensuration metrics."""

from __future__ import annotations

import math

import pytest

from silvispect.inventory import Tree
from silvispect.metrics import (
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

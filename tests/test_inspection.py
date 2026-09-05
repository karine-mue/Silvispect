"""Tests for the inspection rules and report assembly."""

from __future__ import annotations

import json

import pytest

from silvispect.canopy import canopy_height_model
from silvispect.detect import DetectionConfig
from silvispect.grid import Grid, GridError
from silvispect.inspection import (
    RULES,
    InspectionConfig,
    Severity,
    inspect_stand,
)
from silvispect.inventory import Tree


def codes(report) -> set[str]:
    return {finding.code for finding in report.findings}


def flat_canopy(height: float = 20.0, *, nrows: int = 40, ncols: int = 40) -> Grid:
    return Grid.filled(nrows, ncols, height, cellsize=0.5)


def uniform_trees(count: int, *, dbh: float = 30.0, height: float = 25.0, species="PIAB"):
    return [
        Tree(f"t{i:03d}", 1.0 + i % 10, 1.0 + i // 10, species=species, dbh_cm=dbh, height_m=height)
        for i in range(count)
    ]


# ----------------------------------------------------------------------
# plumbing
# ----------------------------------------------------------------------
def test_rules_are_registered():
    assert len(RULES) >= 12


def test_requires_some_input():
    with pytest.raises(ValueError, match="needs a canopy height model"):
        inspect_stand()


def test_requires_area_without_a_raster():
    with pytest.raises(ValueError, match="area_ha"):
        inspect_stand(trees=uniform_trees(10))


def test_area_defaults_to_the_raster_extent():
    report = inspect_stand(chm=flat_canopy(), run_detection=False)
    assert report.area_ha == pytest.approx(0.04)


def test_severity_parsing():
    assert Severity.parse("warning") is Severity.WARNING
    assert Severity.parse(" CRITICAL ") is Severity.CRITICAL
    with pytest.raises(ValueError, match="unknown severity"):
        Severity.parse("catastrophic")
    assert Severity.WARNING.label == "warning"


def test_config_round_trip(tmp_path):
    config = InspectionConfig(min_stems_per_ha=50.0)
    payload = config.as_dict()
    assert payload["min_stems_per_ha"] == 50.0
    restored = InspectionConfig.from_dict(payload)
    assert restored == config

    path = tmp_path / "thresholds.json"
    path.write_text(json.dumps({"min_canopy_cover": 0.9}), encoding="utf-8")
    assert InspectionConfig.load(path).min_canopy_cover == 0.9


def test_config_rejects_unknown_keys():
    with pytest.raises(ValueError, match="unknown configuration keys"):
        InspectionConfig.from_dict({"min_trees": 5})


# ----------------------------------------------------------------------
# individual rules
# ----------------------------------------------------------------------
def test_understocking_is_flagged():
    report = inspect_stand(trees=uniform_trees(10), area_ha=1.0)
    assert "SV001" in codes(report)
    assert report.max_severity >= Severity.WARNING


def test_overstocking_is_flagged():
    report = inspect_stand(trees=uniform_trees(200), area_ha=0.1)
    assert "SV002" in codes(report)


def test_density_index_rules():
    thin = inspect_stand(
        trees=uniform_trees(150, dbh=10.0),
        area_ha=1.0,
        config=InspectionConfig(min_stems_per_ha=1.0, max_stems_per_ha=100_000.0),
    )
    assert "SV003" in codes(thin)

    dense = inspect_stand(
        trees=uniform_trees(900, dbh=40.0),
        area_ha=1.0,
        config=InspectionConfig(max_stems_per_ha=100_000.0),
    )
    assert "SV004" in codes(dense)


def test_low_canopy_cover_and_gap_fraction():
    chm = flat_canopy(20.0)
    for row in range(30):  # knock out three quarters of the canopy
        for col in range(40):
            chm.set(row, col, 0.0)
    report = inspect_stand(chm=chm, run_detection=False)
    assert "SV010" in codes(report)
    assert "SV012" in codes(report)


def test_large_gap_is_flagged():
    chm = flat_canopy(20.0, nrows=80, ncols=80)  # 40 m x 40 m
    for row in range(20, 70):  # a 25 m x 25 m opening, 625 m2
        for col in range(20, 70):
            chm.set(row, col, 0.0)
    report = inspect_stand(chm=chm, run_detection=False)
    gap_findings = [f for f in report.findings if f.code == "SV011"]
    assert gap_findings
    assert gap_findings[0].value > 400.0
    assert report.gaps[0].width > 5.0


def test_small_gaps_are_not_flagged():
    chm = flat_canopy(20.0, nrows=80, ncols=80)
    for row in range(40, 44):  # 2 m x 2 m, below the width criterion
        for col in range(40, 44):
            chm.set(row, col, 0.0)
    assert "SV011" not in codes(inspect_stand(chm=chm, run_detection=False))


def test_uniform_stand_has_low_structural_diversity():
    report = inspect_stand(trees=uniform_trees(300), area_ha=1.0)
    assert "SV020" in codes(report)


def test_species_diversity_rules():
    report = inspect_stand(trees=uniform_trees(300, species="PIAB"), area_ha=1.0)
    assert "SV021" in codes(report)
    assert "SV022" in codes(report)


def test_mixed_stand_passes_diversity(stand):
    report = inspect_stand(trees=stand.trees, area_ha=stand.area_ha)
    assert "SV021" not in codes(report)
    assert "SV022" not in codes(report)


def test_height_outlier_is_flagged(stand):
    trees = list(stand.trees)
    trees.append(Tree("SNAPPED", 5.0, 5.0, species="PIAB", dbh_cm=45.0, height_m=6.0))
    report = inspect_stand(trees=trees, area_ha=stand.area_ha)
    outliers = [f for f in report.findings if f.code == "SV030"]
    assert any(f.subject == "tree:SNAPPED" for f in outliers)


def test_diameter_outlier_is_flagged():
    trees = uniform_trees(40, dbh=30.0)
    trees[0] = Tree("FAT", 1.0, 1.0, species="PIAB", dbh_cm=120.0, height_m=30.0)
    report = inspect_stand(trees=trees, area_ha=1.0)
    assert any(f.code == "SV031" and f.subject == "tree:FAT" for f in report.findings)


def test_missing_measurements_are_reported():
    trees = uniform_trees(20)
    trees.append(Tree("BLANK", 1.0, 1.0, species="PIAB"))
    report = inspect_stand(trees=trees, area_ha=0.05)
    finding = next(f for f in report.findings if f.code == "SV040")
    assert finding.subject == "tree:BLANK"
    assert "dbh_cm" in finding.detail and "height_m" in finding.detail


def test_dead_trees_may_lack_measurements():
    trees = uniform_trees(20)
    trees.append(Tree("SNAG", 1.0, 1.0, species="PIAB", status="dead"))
    report = inspect_stand(trees=trees, area_ha=0.05)
    assert "SV040" not in codes(report)


def test_duplicate_ids_are_critical():
    trees = uniform_trees(20)
    trees.append(Tree("t000", 1.0, 1.0, species="PIAB", dbh_cm=30.0, height_m=25.0))
    report = inspect_stand(trees=trees, area_ha=0.05)
    duplicate = next(f for f in report.findings if f.code == "SV041")
    assert duplicate.severity is Severity.CRITICAL
    assert duplicate.value == 2.0
    assert report.max_severity is Severity.CRITICAL


def test_implausible_dimensions_are_critical():
    trees = uniform_trees(20)
    trees.append(Tree("BAD1", 1.0, 1.0, species="PIAB", dbh_cm=-5.0, height_m=25.0))
    trees.append(Tree("BAD2", 1.0, 1.0, species="PIAB", dbh_cm=30.0, height_m=0.0))
    trees.append(Tree("BAD3", 1.0, 1.0, species="PIAB", dbh_cm=900.0, height_m=500.0))
    report = inspect_stand(trees=trees, area_ha=0.05)
    flagged = {f.subject for f in report.findings if f.code == "SV042"}
    assert flagged == {"tree:BAD1", "tree:BAD2", "tree:BAD3"}


def test_coordinates_outside_the_raster_are_flagged():
    trees = uniform_trees(20)
    trees.append(Tree("FARAWAY", 5000.0, 5000.0, species="PIAB", dbh_cm=30.0, height_m=25.0))
    report = inspect_stand(chm=flat_canopy(), trees=trees, run_detection=False)
    assert any(f.code == "SV043" and f.subject == "tree:FARAWAY" for f in report.findings)


def test_detection_agreement_rules_fire_when_the_reference_is_shifted(stand):
    from silvispect.inventory import shift_trees

    shifted = shift_trees(stand.trees, 15.0, 15.0)
    report = inspect_stand(chm=stand.chm, trees=shifted, area_ha=stand.area_ha)
    assert "SV050" in codes(report)
    assert "SV051" in codes(report)


def test_good_detection_raises_no_agreement_findings(stand):
    report = inspect_stand(chm=stand.chm, trees=stand.trees, area_ha=stand.area_ha)
    assert {"SV050", "SV051", "SV052"} & codes(report) == set()
    assert report.match["recall"] > 0.8


def test_height_bias_rule_fires_when_field_heights_are_wrong(stand):
    from dataclasses import replace

    inflated = [replace(tree, height_m=(tree.height_m or 0) + 6.0) for tree in stand.trees]
    report = inspect_stand(chm=stand.chm, trees=inflated, area_ha=stand.area_ha)
    assert "SV052" in codes(report)


# ----------------------------------------------------------------------
# report assembly
# ----------------------------------------------------------------------
def test_report_without_field_data_uses_detections(stand):
    report = inspect_stand(chm=stand.chm)
    assert report.metrics_source == "detected"
    assert report.detection["tree_count"] > 0
    assert report.match is None
    assert report.metrics.tree_count == report.detection["tree_count"]
    # Species rules need field data and must stay quiet.
    assert {"SV021", "SV022", "SV040", "SV041"} & codes(report) == set()


def test_report_without_a_raster_has_no_canopy_section():
    report = inspect_stand(trees=uniform_trees(300), area_ha=1.0)
    assert report.canopy == {}
    assert report.detection is None
    assert report.gaps == []


def test_report_counts_and_filtering():
    report = inspect_stand(trees=uniform_trees(10), area_ha=1.0)
    counts = report.counts()
    assert set(counts) == {"info", "notice", "warning", "critical"}
    assert sum(counts.values()) == len(report.findings)
    assert report.findings_at_or_above(Severity.CRITICAL) == []
    assert report.findings_at_or_above(Severity.INFO) == report.findings


def test_report_serialises_to_json(stand):
    report = inspect_stand(chm=stand.chm, trees=stand.trees, area_ha=stand.area_ha)
    payload = json.loads(json.dumps(report.as_dict()))
    assert payload["summary"]["finding_count"] == len(report.findings)
    assert payload["metrics"]["stems_per_ha"] > 0
    assert "config" in payload and "canopy" in payload
    assert payload["allometry"]["fitted"] is True
    assert payload["allometry"]["by_species"]


def test_findings_are_sorted_by_severity():
    trees = uniform_trees(10)
    trees.append(Tree("t000", 1.0, 1.0, species="PIAB", dbh_cm=30.0, height_m=25.0))
    report = inspect_stand(trees=trees, area_ha=1.0)
    severities = [int(f.severity) for f in report.findings]
    assert severities == sorted(severities, reverse=True)


def test_finding_dict_shape():
    report = inspect_stand(trees=uniform_trees(10), area_ha=1.0)
    finding = report.findings[0].as_dict()
    assert set(finding) == {"code", "severity", "title", "detail", "subject", "value", "threshold"}


def test_clean_stand_produces_no_findings():
    """A configuration matched to the stand should report nothing."""
    config = InspectionConfig(
        min_stems_per_ha=1.0,
        max_stems_per_ha=100_000.0,
        min_sdi=0.0,
        max_sdi=100_000.0,
        min_gini=0.0,
        min_shannon=0.0,
        max_species_share=1.0,
        height_residual_z=99.0,
        dbh_outlier_z=99.0,
    )
    report = inspect_stand(trees=uniform_trees(300), area_ha=1.0, config=config)
    assert report.findings == []
    assert report.max_severity is Severity.INFO


def test_inspection_accepts_a_derived_canopy_model(stand):
    chm = canopy_height_model(stand.dsm, stand.dtm)
    report = inspect_stand(chm=chm, run_detection=False)
    assert report.canopy["cover"] > 0.0
    assert report.detection is None


# ----------------------------------------------------------------------
# absent inputs must stay absent, not become zero
# ----------------------------------------------------------------------
def test_skipping_detection_does_not_invent_an_empty_stand():
    """No detection and no inventory means unknown stocking, not zero stocking."""
    report = inspect_stand(chm=flat_canopy(), run_detection=False)
    assert report.metrics_source == "none"
    assert report.detection is None
    # Stem-based rules have no stem source and must stay quiet.
    assert {"SV001", "SV002", "SV003", "SV004", "SV020"} & codes(report) == set()


def test_detection_that_finds_nothing_still_reports_understocking():
    """A detector that ran and found nothing IS evidence of an empty stand."""
    report = inspect_stand(chm=flat_canopy(0.0), run_detection=True)
    assert report.metrics_source == "detected"
    assert report.detection["tree_count"] == 0
    assert "SV001" in codes(report)


def test_config_rejects_a_non_object_payload():
    with pytest.raises(ValueError, match="must be a JSON object"):
        InspectionConfig.from_dict(None)  # type: ignore[arg-type]


def test_config_rejects_non_numeric_values():
    with pytest.raises(ValueError, match="must be a number"):
        InspectionConfig.from_dict({"min_sdi": [1, 2]})
    with pytest.raises(ValueError, match="must be a number"):
        InspectionConfig.from_dict({"min_sdi": True})


def test_config_rejects_out_of_range_values():
    with pytest.raises(ValueError, match="min_allometry_points"):
        InspectionConfig(min_allometry_points=0)
    with pytest.raises(ValueError, match="proportion between 0 and 1"):
        InspectionConfig(min_canopy_cover=1.5)
    with pytest.raises(ValueError, match="must not be negative"):
        InspectionConfig(max_gap_area=-1.0)
    with pytest.raises(ValueError, match="must not exceed"):
        InspectionConfig(min_stems_per_ha=900.0, max_stems_per_ha=100.0)


def test_inspection_survives_an_inventory_without_heights():
    """Diameters with no heights must not divide by zero during fitting."""
    trees = [Tree(f"t{i}", i, i, species="PIAB", dbh_cm=30.0) for i in range(20)]
    report = inspect_stand(trees=trees, area_ha=0.05)
    assert report.allometry["fitted"] is False
    assert report.metrics.volume_per_ha_m3 is None
    assert "SV040" in codes(report)


def test_uncounted_stand_reports_unknown_rather_than_zero():
    """metrics_source "none" must not publish a confident zero stocking."""
    report = inspect_stand(chm=flat_canopy(), run_detection=False)
    assert report.metrics_source == "none"
    for value in (
        report.metrics.tree_count,
        report.metrics.measured_count,
        report.metrics.yield_basis_count,
        report.metrics.stems_per_ha,
    ):
        assert value is None
    payload = report.as_dict()["metrics"]
    assert payload["tree_count"] is None
    assert payload["stems_per_ha"] is None


def test_counted_empty_stand_still_reports_zero():
    """Having counted and found nothing is a measurement, and stays numeric."""
    report = inspect_stand(chm=flat_canopy(0.0), run_detection=True)
    assert report.metrics_source == "detected"
    assert report.metrics.tree_count == 0
    assert report.metrics.stems_per_ha == 0.0


def test_config_rejects_non_finite_and_fractional_counts():
    with pytest.raises(ValueError, match="finite"):
        InspectionConfig.from_dict({"max_gap_area": float("inf")})
    with pytest.raises(ValueError, match="finite"):
        InspectionConfig.from_dict({"min_allometry_points": float("inf")})
    with pytest.raises(ValueError, match="whole number"):
        InspectionConfig.from_dict({"min_allometry_points": 1.9})
    with pytest.raises(ValueError, match="finite"):
        InspectionConfig(max_gap_area=float("inf"))
    # A whole number expressed as a float is still a whole number.
    assert InspectionConfig.from_dict({"min_allometry_points": 4.0}).min_allometry_points == 4


def test_an_empty_inventory_is_data_not_absence_of_data():
    """Somebody counted and found nothing — that is a measurement."""
    report = inspect_stand(trees=[], area_ha=1.0)
    assert report.metrics_source == "field"
    assert report.metrics.tree_count == 0
    assert report.metrics.stems_per_ha == 0.0
    assert "SV001" in codes(report)


def test_an_empty_inventory_wins_over_detection_as_the_source():
    report = inspect_stand(chm=flat_canopy(), trees=[])
    assert report.metrics_source == "field"
    assert report.metrics.tree_count == 0


def test_missing_inventory_is_still_missing():
    with pytest.raises(ValueError, match="needs a canopy height model"):
        inspect_stand()


def test_a_limit_is_not_breached_by_sitting_exactly_on_it():
    """SV012 reads "above 30%" strictly, at exactly 30% too.

    Taking the open share as one minus the cover made an exactly-three-tenths
    plot come out as 0.30000000000000004, so a stand sitting precisely on the
    documented limit was reported as over it — and both numbers printed as
    "30%", leaving nothing on screen to explain the finding.
    """
    exactly = Grid.from_rows([[0.0] * 3 + [20.0] * 7])
    codes = {finding.code for finding in inspect_stand(chm=exactly, area_ha=0.001).findings}
    assert "SV012" not in codes

    over = Grid.from_rows([[0.0] * 4 + [20.0] * 6])
    assert "SV012" in {finding.code for finding in inspect_stand(chm=over, area_ha=0.001).findings}

    # The same reading at every limit the option can be set to.
    for open_cells in range(11):
        rows = [[0.0] * open_cells + [20.0] * (10 - open_cells)]
        config = InspectionConfig(max_gap_fraction=open_cells / 10)
        report = inspect_stand(chm=Grid.from_rows(rows), area_ha=0.001, config=config)
        assert "SV012" not in {finding.code for finding in report.findings}


# ----------------------------------------------------------------------
# no observations is not an observation of nothing
# ----------------------------------------------------------------------
def test_a_raster_of_nothing_is_not_an_empty_forest():
    """A raster with no valid cell has measured nothing, so it reports nothing.

    Detection over an all-nodata raster returns no crowns, which was recorded
    as a *count* of zero from "detected" data.  The rules then read that zero
    as a measurement: an understocked stand, and a canopy cover below target —
    two confident findings about a plot nobody had looked at.
    """
    blank = Grid.filled(4, 4, None, cellsize=0.5)
    report = inspect_stand(chm=blank, area_ha=0.01)

    assert report.metrics_source == "none"
    assert report.metrics.tree_count is None
    assert report.metrics.stems_per_ha is None
    assert codes(report) == set()

    payload = report.as_dict()
    assert payload["metrics_source"] == "none"
    assert payload["metrics"]["tree_count"] is None
    assert json.loads(json.dumps(payload, allow_nan=False))["findings"] == []


def test_one_valid_cell_is_enough_to_be_an_observation():
    """The distinction is "nothing was seen", not "little was seen"."""
    values = [None] * 16
    values[5] = 1.0
    sparse = Grid.filled(4, 4, None, cellsize=0.5)
    sparse.values = values
    report = inspect_stand(chm=sparse, area_ha=0.01)
    assert report.metrics_source == "detected"
    assert report.metrics.tree_count == 0
    assert "SV010" in codes(report)


def test_field_data_still_speaks_for_an_unobserved_canopy():
    """Stems counted in the field are evidence even when the raster is blank."""
    trees = [Tree(tree_id=str(i), x=float(i), y=0.0, dbh_cm=20.0, height_m=15.0) for i in range(3)]
    report = inspect_stand(chm=Grid.filled(4, 4, None, cellsize=0.5), trees=trees, area_ha=0.1)
    assert report.metrics_source == "field"
    assert report.metrics.tree_count == 3
    assert "SV001" in codes(report)


# ----------------------------------------------------------------------
# comparator strictness
# ----------------------------------------------------------------------
#: Every rule whose finding records the value it judged and the limit it was
#: judged against, with a setting that forces it to fire and any companion
#: field that has to move out of the way for that setting to be valid.
STRICT_RULES = [
    ("min_stems_per_ha", "SV001", 1e9, {"max_stems_per_ha": 1e12}),
    ("max_stems_per_ha", "SV002", 0.0, {"min_stems_per_ha": 0.0}),
    ("min_sdi", "SV003", 1e9, {"max_sdi": 1e12}),
    ("max_sdi", "SV004", 0.0, {"min_sdi": 0.0}),
    ("min_canopy_cover", "SV010", 1.0, {}),
    ("max_gap_area", "SV011", 0.0, {"min_gap_area": 0.0}),
    ("max_gap_fraction", "SV012", 0.0, {}),
    ("min_gini", "SV020", 1e9, {}),
    ("min_shannon", "SV021", 1e9, {}),
    ("max_species_share", "SV022", 0.0, {}),
    ("max_height_bias_m", "SV052", 0.0, {}),
]


@pytest.mark.parametrize(
    "field,code,forcing,companions", STRICT_RULES, ids=[r[1] for r in STRICT_RULES]
)
def test_no_rule_fires_on_a_value_sitting_exactly_on_its_limit(
    stand, field, code, forcing, companions
):
    """The rule table promises "below", "above" or "more than" everywhere.

    Each rule is first forced to fire so the stand's own value can be read off
    the finding, and then re-run with the limit set to exactly that value.  A
    strict comparison goes quiet; SV030 and SV031 did not, which made a
    threshold of zero standard deviations report a stem sitting exactly on the
    mean as an outlier.

    Reading the value from the finding is what keeps this test honest as the
    fixture changes: there is no table of expected numbers to drift.
    """

    def report_for(value):
        config = InspectionConfig(**{field: value}, **companions)
        return inspect_stand(chm=stand.chm, trees=stand.trees, config=config)

    forced = [f for f in report_for(forcing).findings if f.code == code and f.value is not None]
    assert forced, f"{code} could not be provoked, so its comparator is untested"

    observed = forced[0].value
    on_the_limit = report_for(observed)
    assert not [f for f in on_the_limit.findings if f.code == code and f.value == f.threshold], (
        f"{code} fires on a value equal to its own limit"
    )


def test_perfect_agreement_is_not_below_a_target_of_perfect_agreement(stand):
    """SV050 and SV051 cannot be provoked on a stand they match perfectly.

    That is the on-the-limit case itself: recall of 1.0 against a target of
    1.0 is not "fewer than", and a target it really is under must still speak.
    """
    perfect = inspect_stand(
        chm=stand.chm,
        trees=stand.trees,
        config=InspectionConfig(min_recall=1.0, min_precision=1.0),
    )
    assert perfect.match["recall"] == perfect.match["precision"] == 1.0
    assert {"SV050", "SV051"}.isdisjoint(codes(perfect))

    partial = inspect_stand(
        chm=stand.chm,
        trees=[*stand.trees[:-1], Tree(tree_id="ghost", x=1.0, y=1.0, height_m=25.0)],
        config=InspectionConfig(min_recall=1.0, min_precision=1.0),
    )
    assert "SV050" in codes(partial)


def test_no_rule_fires_on_a_z_score_sitting_exactly_on_its_limit():
    """SV030 and SV031 are documented as "more than z", so exactly z is not.

    They carry no companion field to trade against, so they get their own
    case: five diameters symmetric about their mean, one of which therefore
    has a z-score of exactly zero.
    """
    diameters = [10.0, 20.0, 30.0, 40.0, 50.0]
    trees = [
        Tree(tree_id=f"T{i}", x=float(i), y=0.0, dbh_cm=dbh, height_m=20.0)
        for i, dbh in enumerate(diameters)
    ]
    config = InspectionConfig(dbh_outlier_z=0.0, height_residual_z=0.0)
    report = inspect_stand(trees=trees, area_ha=0.1, config=config)
    assert [f for f in report.findings if f.code in {"SV030", "SV031"} and f.value == 0.0] == []


# ----------------------------------------------------------------------
# configuration domain
# ----------------------------------------------------------------------
@pytest.mark.parametrize("field", ["min_canopy_cover", "min_stems_per_ha", "min_allometry_points"])
@pytest.mark.parametrize("value", [True, False, None])
def test_a_profile_has_the_same_domain_however_it_is_built(field, value):
    """Parsing a profile and writing one in Python must accept the same things.

    ``from_dict`` refused booleans and non-numbers; the constructor did not, so
    ``InspectionConfig(min_canopy_cover=True)`` built a profile whose cover
    threshold was ``True`` — it compared as 1.0, printed as ``True`` in the
    report, and could not have come from a profile file.
    """
    with pytest.raises(ValueError):
        InspectionConfig(**{field: value})
    with pytest.raises(ValueError):
        InspectionConfig.from_dict({field: value})


def test_a_parsed_number_is_stored_as_a_number():
    """Parsing may read "0.5" from a file; it may not leave a string behind.

    That is the direction the two doors legitimately differ in — one of them
    is a parser — so what has to match is the domain of the field once it is
    set, not the shape of the argument.
    """
    assert InspectionConfig.from_dict({"min_canopy_cover": "0.5"}).min_canopy_cover == 0.5
    with pytest.raises(ValueError):
        InspectionConfig(min_canopy_cover="0.5")


def test_detection_configuration_has_a_domain_too():
    """The same parity for the detection knobs, which count whole cells."""
    for bad in ({"smooth_radius": True}, {"min_crown_cells": 1.5}, {"min_height": float("inf")}):
        with pytest.raises(GridError):
            DetectionConfig(**bad)
    assert DetectionConfig(smooth_radius=0, min_crown_cells=2).min_crown_cells == 2


# ----------------------------------------------------------------------
# a finding must not contradict itself
# ----------------------------------------------------------------------
def test_a_finding_never_prints_the_two_numbers_it_compared_as_one():
    """ "200 stems/ha is below the minimum of 200 stems/ha" explains nothing.

    Density is reported to the nearest stem because tenths of a stem per
    hectare are noise, but a stand just under the limit then printed the same
    number twice and read as a contradiction.  Precision is added only where
    it is needed to keep the sentence true.
    """
    trees = [
        Tree(tree_id=str(i), x=float(i), y=0.0, dbh_cm=20.0, height_m=15.0) for i in range(100)
    ]
    report = inspect_stand(
        trees=trees, area_ha=100 / 199.6, config=InspectionConfig(min_stems_per_ha=200.0)
    )
    detail = next(f.detail for f in report.findings if f.code == "SV001")
    assert "199.6 stems/ha is below the minimum of 200.0 stems/ha" in detail

    # An ordinary finding keeps the resolution the rule chose.
    plain = inspect_stand(trees=trees, area_ha=1.0, config=InspectionConfig(min_stems_per_ha=200.0))
    assert "100 stems/ha is below the minimum of 200 stems/ha" in next(
        f.detail for f in plain.findings if f.code == "SV001"
    )


# ----------------------------------------------------------------------
# nothing observed means nothing measured, everywhere
# ----------------------------------------------------------------------
def test_an_unobserved_raster_makes_no_claim_about_a_field_inventory():
    """Matching a field list against an empty detection is not a comparison.

    Detection over an all-nodata raster returns no crowns, and matching those
    against a real inventory reported every stem as an omission — a warning
    that the sensor missed them, from a sensor that had never seen anything.
    On a ``--fail-on warning`` run that failed the build.
    """
    blank = Grid.filled(2, 2, None, cellsize=0.5)
    trees = [Tree(tree_id="T", x=0.5, y=0.5, dbh_cm=10.0, height_m=10.0)]
    report = inspect_stand(chm=blank, trees=trees)

    assert report.match is None
    assert {"SV050", "SV051", "SV052"}.isdisjoint(codes(report))
    assert report.metrics_source == "field"


def test_an_unobserved_raster_reports_no_canopy_numbers():
    """Cover, open fraction, rugosity and gap count are observations too.

    Each of them came back as a numeric zero for a plot nobody had looked at,
    which reads as a measured absence rather than a missing measurement.
    """
    blank = Grid.filled(4, 4, None, cellsize=0.5)
    report = inspect_stand(chm=blank, area_ha=0.01)

    payload = report.as_dict()
    assert payload["canopy"] == {}
    assert "detection" not in payload
    assert "match" not in payload
    assert "gaps" not in payload
    assert report.gaps == []
    assert json.loads(json.dumps(payload, allow_nan=False))["metrics_source"] == "none"

    # One valid cell is an observation, and then every block is populated.
    seen = Grid.filled(4, 4, None, cellsize=0.5)
    seen.set(1, 1, 1.0)
    populated = inspect_stand(chm=seen, area_ha=0.01).as_dict()
    assert populated["canopy"]["cover"] == 0.0
    assert "detection" in populated


def test_a_reporting_floor_may_not_sit_above_the_limit_it_feeds():
    """SV011 cannot report a gap that ``min_gap_area`` has already deleted.

    A floor of 101 m2 against a limit of 50 m2 removed a 100 m2 opening before
    the rule saw it, and the inspection came back clean — the configuration
    silenced exactly the finding it was meant to raise.
    """
    with pytest.raises(ValueError, match="min_gap_area must not exceed max_gap_area"):
        InspectionConfig(min_gap_area=101.0, max_gap_area=50.0)
    with pytest.raises(ValueError, match="min_gap_area must not exceed max_gap_area"):
        InspectionConfig.from_dict({"min_gap_area": 101.0, "max_gap_area": 50.0})

    # The same plot, with a floor the limit can live with, is critical.
    plot = Grid.from_rows([[9.0, 9.0, 9.0], [9.0, 0.0, 9.0], [9.0, 9.0, 9.0]], cellsize=10.0)
    config = InspectionConfig(min_gap_area=50.0, max_gap_area=50.0, min_gap_width=0.0)
    report = inspect_stand(chm=plot, config=config, run_detection=False)
    assert "SV011" in codes(report)
    assert report.max_severity is Severity.CRITICAL


def test_coverage_is_judged_cell_by_cell_not_plot_by_plot():
    """A hole in the coverage is not a stem the detector missed.

    The whole-raster check only recognised a raster with no valid cell at all,
    so every partly-covered plot went on being judged as though it were fully
    covered: a plot observed in one corner reported the stems standing in its
    nodata as detector omissions, and SV050 warned about a sensor that had
    never looked there.
    """
    chm = Grid.filled(5, 5, None, cellsize=1.0)
    for row, col in ((0, 0), (0, 1), (1, 0)):
        chm.set(row, col, 10.0)
    seen = Tree(tree_id="seen", x=0.5, y=4.5, dbh_cm=20.0, height_m=10.0)
    blind = Tree(tree_id="blind", x=4.5, y=0.5, dbh_cm=20.0, height_m=10.0)

    report = inspect_stand(chm=chm, trees=[seen, blind])
    assert report.match["matched"] == 1
    assert report.match["omissions"] == 0
    assert report.match["recall"] == 1.0
    assert "SV050" not in codes(report)

    # A stem standing in the observed corner and genuinely undetected still
    # counts against recall — the gate is coverage, not convenience.
    missed = Tree(tree_id="missed", x=1.5, y=4.5, dbh_cm=20.0, height_m=10.0)
    chm.set(0, 1, 10.0)
    with_hole = inspect_stand(chm=chm, trees=[seen, blind, missed])
    assert with_hole.match["omissions"] == 1
    assert with_hole.match["recall"] == pytest.approx(0.5)


def test_a_stem_on_the_boundary_of_the_coverage_is_judged_by_its_own_cell():
    """Eligibility follows the cell a stem stands in, either side of the edge."""
    chm = Grid.filled(1, 4, None, cellsize=1.0)
    chm.set(0, 0, 12.0)
    chm.set(0, 1, 12.0)
    trees = [
        Tree(tree_id="c0", x=0.5, y=0.5, dbh_cm=20.0, height_m=12.0),
        Tree(tree_id="c1", x=1.5, y=0.5, dbh_cm=20.0, height_m=12.0),
        Tree(tree_id="c2", x=2.5, y=0.5, dbh_cm=20.0, height_m=12.0),
        Tree(tree_id="c3", x=3.5, y=0.5, dbh_cm=20.0, height_m=12.0),
    ]
    report = inspect_stand(chm=chm, trees=trees)
    # Only the two stems standing on observed cells can be judged at all.
    assert report.match["matched"] + report.match["omissions"] == 2

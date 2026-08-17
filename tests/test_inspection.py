"""Tests for the inspection rules and report assembly."""

from __future__ import annotations

import json

import pytest

from silvispect.canopy import canopy_height_model
from silvispect.grid import Grid
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

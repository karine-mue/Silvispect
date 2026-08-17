"""End-to-end tests for the command line interface."""

from __future__ import annotations

import io
import json

import pytest

from silvispect.cli import EXIT_ERROR, EXIT_FINDINGS, EXIT_OK, main
from silvispect.grid import Grid
from silvispect.inventory import read_trees, write_trees


def run(*argv) -> tuple[int, str]:
    buffer = io.StringIO()
    code = main([str(a) for a in argv], out=buffer)
    return code, buffer.getvalue()


@pytest.fixture(scope="module")
def workspace(tmp_path_factory, stand):
    directory = tmp_path_factory.mktemp("stand")
    stand.chm.write(directory / "chm.asc", precision=2)
    stand.dsm.write(directory / "dsm.asc", precision=2)
    stand.dtm.write(directory / "dtm.asc", precision=2)
    write_trees(stand.trees, directory / "trees.csv", precision=2)
    return directory


def test_version_and_help():
    with pytest.raises(SystemExit) as exit_info:
        main(["--version"])
    assert exit_info.value.code == 0
    with pytest.raises(SystemExit):
        main([])


def test_chm_command(workspace, tmp_path):
    out_path = tmp_path / "derived.asc"
    code, text = run(
        "chm", "--dsm", workspace / "dsm.asc", "--dtm", workspace / "dtm.asc", "-o", out_path
    )
    assert code == EXIT_OK
    assert "canopy height model:" in text
    assert out_path.exists()
    assert Grid.read(out_path).stats().minimum >= 0.0


def test_detect_command_writes_trees(workspace, tmp_path):
    trees_path = tmp_path / "detected.csv"
    labels_path = tmp_path / "labels.asc"
    code, text = run("detect", workspace / "chm.asc", "-o", trees_path, "--labels", labels_path)
    assert code == EXIT_OK
    assert "detected" in text
    detected = read_trees(trees_path)
    assert len(detected) > 10
    assert all(tree.height_m and tree.crown_diameter_m for tree in detected)
    assert Grid.read(labels_path).stats().maximum == float(len(detected))


def test_detect_command_json(workspace):
    code, text = run("detect", workspace / "chm.asc", "--json")
    payload = json.loads(text)
    assert code == EXIT_OK
    assert payload["tree_count"] == len(payload["trees"])
    assert payload["config"]["min_height"] == 2.0


def test_detect_respects_min_height(workspace):
    _, tall = run("detect", workspace / "chm.asc", "--min-height", 25.0, "--json")
    _, low = run("detect", workspace / "chm.asc", "--min-height", 2.0, "--json")
    assert json.loads(tall)["tree_count"] < json.loads(low)["tree_count"]


def test_metrics_command(workspace, stand):
    code, text = run("metrics", workspace / "trees.csv", "--area-ha", stand.area_ha)
    assert code == EXIT_OK
    assert "stems_per_ha" in text
    code, payload = run("metrics", workspace / "trees.csv", "--area-ha", 1.0, "--json")
    assert json.loads(payload)["tree_count"] == len(stand.trees)


def test_gaps_command(workspace):
    code, text = run("gaps", workspace / "chm.asc", "--min-area", 20.0)
    assert code == EXIT_OK
    assert "gaps at or above" in text
    code, payload = run("gaps", workspace / "chm.asc", "--json")
    assert isinstance(json.loads(payload), list)


def test_gaps_width_filter_reduces_merging(workspace):
    _, raw = run("gaps", workspace / "chm.asc", "--min-width", 0.0, "--json")
    _, opened = run("gaps", workspace / "chm.asc", "--min-width", 5.0, "--json")
    raw_gaps, opened_gaps = json.loads(raw), json.loads(opened)
    if raw_gaps:
        assert raw_gaps[0]["area_m2"] >= opened_gaps[0]["area_m2"]
        assert all(gap["width_m"] >= 5.0 for gap in opened_gaps)


def test_render_command(workspace):
    code, text = run("render", workspace / "chm.asc", "--width", 40)
    assert code == EXIT_OK
    assert "low" in text.splitlines()[-2] or "low" in text.splitlines()[-1]
    code, with_tops = run("render", workspace / "chm.asc", "--width", 40, "--tops")
    assert "^" in with_tops


def test_match_command(workspace):
    code, text = run("match", workspace / "trees.csv", "--chm", workspace / "chm.asc")
    assert code == EXIT_OK
    assert "recall" in text
    code, payload = run("match", workspace / "trees.csv", "--chm", workspace / "chm.asc", "--json")
    assert json.loads(payload)["recall"] > 0.5


def test_allometry_command(workspace):
    code, text = run("allometry", workspace / "trees.csv")
    assert code == EXIT_OK
    assert "r-squared" in text
    code, payload = run("allometry", workspace / "trees.csv", "--model", "curtis", "--json")
    assert json.loads(payload)["model"] == "curtis"


def test_inspect_text_and_markdown(workspace, stand):
    code, text = run(
        "inspect", "--chm", workspace / "chm.asc", "--inventory", workspace / "trees.csv"
    )
    assert code == EXIT_OK
    assert "Silvispect inspection" in text

    code, markdown = run(
        "inspect",
        "--chm",
        workspace / "chm.asc",
        "--inventory",
        workspace / "trees.csv",
        "--format",
        "markdown",
        "--with-map",
    )
    assert "# Silvispect stand inspection" in markdown
    assert "## Canopy map" in markdown


def test_inspect_json_and_output_file(workspace, tmp_path):
    target = tmp_path / "report.json"
    code, text = run("inspect", "--chm", workspace / "chm.asc", "--format", "json", "-o", target)
    assert code == EXIT_OK
    assert "wrote" in text
    payload = json.loads(target.read_text())
    assert payload["metrics_source"] == "detected"


def test_inspect_inventory_only_requires_area(workspace):
    code, _ = run("inspect", "--inventory", workspace / "trees.csv")
    assert code == EXIT_ERROR
    code, text = run("inspect", "--inventory", workspace / "trees.csv", "--area-ha", 0.36)
    assert code == EXIT_OK
    assert "field data" in text


def test_inspect_needs_an_input():
    with pytest.raises(SystemExit):
        main(["inspect"])


def test_inspect_fail_on_severity(workspace):
    code, _ = run(
        "inspect",
        "--inventory",
        workspace / "trees.csv",
        "--area-ha",
        0.001,  # absurd area makes the stand look wildly overstocked
        "--fail-on",
        "notice",
    )
    assert code == EXIT_FINDINGS

    code, _ = run(
        "inspect", "--chm", workspace / "chm.asc", "--fail-on", "critical", "--no-detection"
    )
    assert code == EXIT_OK


def test_inspect_with_a_config_file(workspace, tmp_path):
    config = tmp_path / "thresholds.json"
    # Above the sample plot's 340 stems/ha, but still under the default
    # maximum, so the profile stays internally consistent.
    config.write_text(json.dumps({"min_stems_per_ha": 500.0}), encoding="utf-8")
    code, text = run(
        "inspect",
        "--inventory",
        workspace / "trees.csv",
        "--area-ha",
        0.36,
        "--config",
        config,
    )
    assert code == EXIT_OK
    assert "SV001" in text


def test_inspect_rejects_a_bad_config(workspace, tmp_path):
    config = tmp_path / "bad.json"
    config.write_text(json.dumps({"nonsense": 1.0}), encoding="utf-8")
    code, _ = run(
        "inspect", "--inventory", workspace / "trees.csv", "--area-ha", 1.0, "--config", config
    )
    assert code == EXIT_ERROR


def test_synth_command(tmp_path):
    target = tmp_path / "generated"
    code, text = run("synth", "--out", target, "--width", 40, "--height", 40, "--seed", 5)
    assert code == EXIT_OK
    assert "generated" in text
    for name in ("chm.asc", "dsm.asc", "dtm.asc", "trees.csv"):
        assert (target / name).exists()
    assert len(read_trees(target / "trees.csv")) > 0


def test_synth_is_reproducible(tmp_path):
    first, second = tmp_path / "a", tmp_path / "b"
    run("synth", "--out", first, "--width", 30, "--height", 30, "--seed", 42)
    run("synth", "--out", second, "--width", 30, "--height", 30, "--seed", 42)
    assert (first / "chm.asc").read_text() == (second / "chm.asc").read_text()
    assert (first / "trees.csv").read_text() == (second / "trees.csv").read_text()


def test_demo_command():
    code, text = run("demo", "--seed", 5)
    assert code == EXIT_OK
    assert "Silvispect inspection" in text
    code, payload = run("demo", "--seed", 5, "--format", "json")
    assert json.loads(payload)["metrics_source"] == "field"


def test_missing_file_is_reported(tmp_path):
    code, _ = run("detect", tmp_path / "absent.asc")
    assert code == EXIT_ERROR


def test_malformed_raster_is_reported(tmp_path):
    broken = tmp_path / "broken.asc"
    broken.write_text("ncols 2\nnrows 2\ncellsize 1\n1 2 3\n", encoding="utf-8")
    code, _ = run("detect", broken)
    assert code == EXIT_ERROR


def test_malformed_inventory_is_reported(tmp_path):
    broken = tmp_path / "broken.csv"
    broken.write_text("tree_id,dbh_cm\nA,30\n", encoding="utf-8")
    code, _ = run("metrics", broken, "--area-ha", 1.0)
    assert code == EXIT_ERROR

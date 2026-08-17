"""Command line interface for Silvispect.

Run ``silvispect --help`` for the full list of subcommands.  Every command
reads and writes plain text, so the tool composes with ordinary shell
pipelines and with CI systems: ``silvispect inspect --fail-on warning``
returns a non-zero exit status when a stand breaks a rule.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from collections.abc import Sequence
from pathlib import Path

from . import __version__
from .allometry import DEFAULT_MODEL, MODELS, AllometryError, fit_height_diameter
from .canopy import canopy_cover, canopy_height_model, find_gaps
from .detect import DetectionConfig, detect_trees
from .grid import Grid, GridError
from .inspection import InspectionConfig, Severity, inspect_stand
from .inventory import InventoryError, match_trees, read_trees, trees_from_crowns, write_trees
from .metrics import stand_metrics
from .report import ascii_map, crown_markers, render_json, render_markdown, render_text
from .synth import SynthConfig, synthesize

EXIT_OK = 0
EXIT_FINDINGS = 1
EXIT_ERROR = 2


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="silvispect",
        description="Inspect forest canopy models and inventories.",
    )
    parser.add_argument("--version", action="version", version=f"silvispect {__version__}")
    subparsers = parser.add_subparsers(dest="command", required=True)

    chm = subparsers.add_parser("chm", help="derive a canopy height model from DSM and DTM")
    chm.add_argument("--dsm", required=True, type=Path, help="digital surface model (.asc)")
    chm.add_argument("--dtm", required=True, type=Path, help="digital terrain model (.asc)")
    chm.add_argument("-o", "--output", type=Path, help="write the canopy height model here")
    chm.set_defaults(handler=_cmd_chm)

    detect = subparsers.add_parser("detect", help="detect trees and delineate crowns")
    detect.add_argument("chm", type=Path, help="canopy height model (.asc)")
    detect.add_argument("--min-height", type=float, default=2.0, help="minimum tree height")
    detect.add_argument("--smooth-radius", type=int, default=1, help="pre-smoothing radius")
    detect.add_argument("--window-intercept", type=float, default=0.8)
    detect.add_argument("--window-slope", type=float, default=0.055)
    detect.add_argument("-o", "--output", type=Path, help="write detected trees as CSV")
    detect.add_argument("--labels", type=Path, help="write the crown label raster here")
    detect.add_argument("--json", action="store_true", help="emit JSON instead of text")
    detect.set_defaults(handler=_cmd_detect)

    metrics = subparsers.add_parser("metrics", help="summarise an inventory")
    metrics.add_argument("inventory", type=Path, help="tree list (.csv)")
    metrics.add_argument("--area-ha", type=float, required=True, help="plot area in hectares")
    metrics.add_argument("--include-dead", action="store_true", help="keep non-live stems")
    metrics.add_argument("--json", action="store_true")
    metrics.set_defaults(handler=_cmd_metrics)

    gaps = subparsers.add_parser("gaps", help="map canopy openings")
    gaps.add_argument("chm", type=Path)
    gaps.add_argument("--threshold", type=float, default=2.0, help="canopy height threshold")
    gaps.add_argument("--min-area", type=float, default=25.0, help="smallest reported gap")
    gaps.add_argument(
        "--min-width", type=float, default=5.0, help="narrowest reported gap (Brokaw criterion)"
    )
    gaps.add_argument("--json", action="store_true")
    gaps.set_defaults(handler=_cmd_gaps)

    render = subparsers.add_parser("render", help="draw a raster as an ASCII map")
    render.add_argument("chm", type=Path)
    render.add_argument("--width", type=int, default=72)
    render.add_argument("--tops", action="store_true", help="overlay detected treetops")
    render.set_defaults(handler=_cmd_render)

    match = subparsers.add_parser("match", help="compare detections with field trees")
    match.add_argument("inventory", type=Path, help="field tree list (.csv)")
    match.add_argument("--chm", required=True, type=Path, help="canopy height model (.asc)")
    match.add_argument("--tolerance", type=float, default=2.5, help="matching radius in metres")
    match.add_argument("--json", action="store_true")
    match.set_defaults(handler=_cmd_match)

    allometry = subparsers.add_parser("allometry", help="fit a height-diameter curve")
    allometry.add_argument("inventory", type=Path)
    allometry.add_argument("--model", default=DEFAULT_MODEL, choices=sorted(MODELS))
    allometry.add_argument("--json", action="store_true")
    allometry.set_defaults(handler=_cmd_allometry)

    inspect_cmd = subparsers.add_parser("inspect", help="run the inspection rules")
    inspect_cmd.add_argument("--chm", type=Path, help="canopy height model (.asc)")
    inspect_cmd.add_argument("--inventory", type=Path, help="field tree list (.csv)")
    inspect_cmd.add_argument("--area-ha", type=float, help="plot area, defaults to the raster")
    inspect_cmd.add_argument("--config", type=Path, help="threshold overrides (.json)")
    inspect_cmd.add_argument(
        "--format", default="text", choices=("text", "markdown", "json"), help="output format"
    )
    inspect_cmd.add_argument("-o", "--output", type=Path, help="write the report here")
    inspect_cmd.add_argument("--no-detection", action="store_true", help="skip tree detection")
    inspect_cmd.add_argument("--with-map", action="store_true", help="append an ASCII map")
    inspect_cmd.add_argument(
        "--fail-on",
        default="none",
        choices=("none", "notice", "warning", "critical"),
        help="exit non-zero when a finding reaches this severity",
    )
    inspect_cmd.set_defaults(handler=_cmd_inspect)

    synth = subparsers.add_parser("synth", help="generate a synthetic stand")
    synth.add_argument("--out", type=Path, default=Path("stand"), help="output directory")
    synth.add_argument("--seed", type=int, default=SynthConfig.seed)
    synth.add_argument("--width", type=float, default=SynthConfig.width)
    synth.add_argument("--height", type=float, default=SynthConfig.height)
    synth.add_argument("--cellsize", type=float, default=SynthConfig.cellsize)
    synth.add_argument("--stems-per-ha", type=float, default=SynthConfig.stems_per_ha)
    synth.add_argument("--gaps", type=int, default=SynthConfig.gap_count)
    synth.set_defaults(handler=_cmd_synth)

    demo = subparsers.add_parser("demo", help="generate a stand and inspect it end to end")
    demo.add_argument("--seed", type=int, default=SynthConfig.seed)
    demo.add_argument("--format", default="text", choices=("text", "markdown", "json"))
    demo.set_defaults(handler=_cmd_demo)

    return parser


# ----------------------------------------------------------------------
# command handlers
# ----------------------------------------------------------------------
def _cmd_chm(args: argparse.Namespace, out) -> int:
    dsm = Grid.read(args.dsm)
    dtm = Grid.read(args.dtm)
    chm = canopy_height_model(dsm, dtm)
    if args.output:
        chm.write(args.output)
        print(f"wrote {args.output}", file=out)
    stats = chm.stats()
    print(
        f"canopy height model: {chm.nrows} x {chm.ncols} cells at {chm.cellsize} m, "
        f"mean {_fmt(stats.mean)} m, max {_fmt(stats.maximum)} m, "
        f"cover {canopy_cover(chm):.1%}",
        file=out,
    )
    return EXIT_OK


def _cmd_detect(args: argparse.Namespace, out) -> int:
    chm = Grid.read(args.chm)
    config = DetectionConfig(
        min_height=args.min_height,
        smooth_radius=args.smooth_radius,
        window_intercept=args.window_intercept,
        window_slope=args.window_slope,
    )
    result = detect_trees(chm, config)
    if args.output:
        write_trees(trees_from_crowns(result.crowns), args.output)
    if args.labels:
        result.label_grid.write(args.labels)
    if args.json:
        print(json.dumps(result.as_dict(), indent=2), file=out)
        return EXIT_OK
    print(
        f"detected {len(result.crowns)} trees "
        f"({result.density_per_ha():.0f} stems/ha over {chm.area_ha:.3f} ha)",
        file=out,
    )
    if result.crowns:
        heights = [crown.height for crown in result.crowns]
        crowns = [crown.diameter for crown in result.crowns]
        print(
            f"height  min {min(heights):.1f} m  mean "
            f"{sum(heights) / len(heights):.1f} m  max {max(heights):.1f} m",
            file=out,
        )
        print(
            f"crown   min {min(crowns):.1f} m  mean "
            f"{sum(crowns) / len(crowns):.1f} m  max {max(crowns):.1f} m",
            file=out,
        )
    if args.output:
        print(f"wrote {args.output}", file=out)
    if args.labels:
        print(f"wrote {args.labels}", file=out)
    return EXIT_OK


def _cmd_metrics(args: argparse.Namespace, out) -> int:
    trees = read_trees(args.inventory)
    metrics = stand_metrics(trees, args.area_ha, live_only=not args.include_dead)
    if args.json:
        print(json.dumps(metrics.as_dict(), indent=2), file=out)
        return EXIT_OK
    for key, value in metrics.as_dict().items():
        if isinstance(value, dict):
            if not value:
                continue
            print(f"{key}:", file=out)
            for sub_key, sub_value in value.items():
                print(f"  {sub_key:<12} {_fmt(sub_value)}", file=out)
        else:
            print(f"{key:<22} {_fmt(value)}", file=out)
    return EXIT_OK


def _cmd_gaps(args: argparse.Namespace, out) -> int:
    chm = Grid.read(args.chm)
    gaps = find_gaps(
        chm, threshold=args.threshold, min_area=args.min_area, min_width=args.min_width
    )
    if args.json:
        print(json.dumps([gap.as_dict() for gap in gaps], indent=2), file=out)
        return EXIT_OK
    open_fraction = 1.0 - canopy_cover(chm, args.threshold)
    print(
        f"{len(gaps)} gaps at or above {args.min_area:.0f} m2 and {args.min_width:.1f} m wide; "
        f"{open_fraction:.1%} of the plot is below {args.threshold:.1f} m",
        file=out,
    )
    for gap in gaps:
        print(
            f"  gap {gap.gap_id:<3} {gap.area:8.0f} m2  width {gap.width:5.1f} m  "
            f"centre ({gap.centroid_x:7.1f}, {gap.centroid_y:7.1f})"
            f"{'  [edge]' if gap.touches_edge else ''}",
            file=out,
        )
    return EXIT_OK


def _cmd_render(args: argparse.Namespace, out) -> int:
    chm = Grid.read(args.chm)
    markers: list[tuple[float, float, str]] = []
    if args.tops:
        markers = crown_markers(detect_trees(chm).crowns)
    print(ascii_map(chm, width=args.width, markers=markers), file=out)
    return EXIT_OK


def _cmd_match(args: argparse.Namespace, out) -> int:
    chm = Grid.read(args.chm)
    reference = read_trees(args.inventory)
    result = detect_trees(chm)
    match = match_trees(result.crowns, reference, tolerance=args.tolerance)
    if args.json:
        print(json.dumps(match.as_dict(), indent=2), file=out)
        return EXIT_OK
    for key, value in match.as_dict().items():
        print(f"{key:<16} {_fmt(value)}", file=out)
    return EXIT_OK


def _cmd_allometry(args: argparse.Namespace, out) -> int:
    trees = read_trees(args.inventory)
    model = fit_height_diameter(trees, model=args.model)
    if args.json:
        print(json.dumps(model.as_dict(), indent=2), file=out)
        return EXIT_OK
    print(f"model      {model.name}", file=out)
    print(f"equation   {model.spec.description}", file=out)
    for name, value in zip(model.spec.param_names, model.params, strict=True):
        print(f"  {name:<8} {value:.4f}", file=out)
    print(f"n          {model.n}", file=out)
    print(f"rmse       {model.rmse:.3f} m", file=out)
    print(f"r-squared  {model.r_squared:.4f}", file=out)
    return EXIT_OK


def _cmd_inspect(args: argparse.Namespace, out) -> int:
    if not args.chm and not args.inventory:
        raise SystemExit("inspect needs --chm, --inventory, or both")
    chm = Grid.read(args.chm) if args.chm else None
    trees = read_trees(args.inventory) if args.inventory else None
    config = InspectionConfig.load(args.config) if args.config else InspectionConfig()
    report = inspect_stand(
        chm=chm,
        trees=trees,
        area_ha=args.area_ha,
        config=config,
        run_detection=not args.no_detection,
    )
    if args.format == "json":
        text = render_json(report)
    elif args.format == "markdown":
        text = render_markdown(report, chm=chm if args.with_map else None)
    else:
        text = render_text(report)
        if args.with_map and chm is not None:
            text += "\n" + ascii_map(chm) + "\n"
    if args.output:
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        Path(args.output).write_text(text, encoding="utf-8")
        print(f"wrote {args.output}", file=out)
    else:
        print(text, end="", file=out)

    if args.fail_on != "none":
        threshold = Severity.parse(args.fail_on)
        if report.findings_at_or_above(threshold):
            return EXIT_FINDINGS
    return EXIT_OK


def _cmd_synth(args: argparse.Namespace, out) -> int:
    config = SynthConfig(
        width=args.width,
        height=args.height,
        cellsize=args.cellsize,
        stems_per_ha=args.stems_per_ha,
        seed=args.seed,
        gap_count=args.gaps,
    )
    stand = synthesize(config)
    directory = Path(args.out)
    directory.mkdir(parents=True, exist_ok=True)
    stand.chm.write(directory / "chm.asc", precision=2)
    stand.dsm.write(directory / "dsm.asc", precision=2)
    stand.dtm.write(directory / "dtm.asc", precision=2)
    write_trees(stand.trees, directory / "trees.csv", precision=2)
    print(
        f"generated {len(stand.trees)} trees over {stand.area_ha:.3f} ha "
        f"({len(stand.trees) / stand.area_ha:.0f} stems/ha) in {directory}/",
        file=out,
    )
    return EXIT_OK


def _cmd_demo(args: argparse.Namespace, out) -> int:
    stand = synthesize(SynthConfig(seed=args.seed))
    report = inspect_stand(chm=stand.chm, trees=stand.trees, area_ha=stand.area_ha)
    if args.format == "json":
        print(render_json(report), end="", file=out)
    elif args.format == "markdown":
        print(render_markdown(report, chm=stand.chm), end="", file=out)
    else:
        print(render_text(report), end="", file=out)
        print(ascii_map(stand.chm, width=72, markers=[]), file=out)
    return EXIT_OK


def main(argv: Sequence[str] | None = None, *, out=None) -> int:
    """Entry point.  Returns a process exit status."""
    out = out or sys.stdout
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return args.handler(args, out)
    except (GridError, InventoryError, AllometryError, ValueError) as exc:
        print(f"silvispect: {exc}", file=sys.stderr)
        return EXIT_ERROR
    except FileNotFoundError as exc:
        print(f"silvispect: {exc.filename}: file not found", file=sys.stderr)
        return EXIT_ERROR
    except BrokenPipeError:
        # A downstream reader such as `head` closed the pipe.  Redirect stdout
        # to devnull so the interpreter does not report the failure again while
        # flushing at exit, which is what well-behaved text tools do.
        devnull = os.open(os.devnull, os.O_WRONLY)
        os.dup2(devnull, sys.stdout.fileno())
        return EXIT_OK


def _fmt(value: object) -> str:
    if value is None:
        return "-"
    if isinstance(value, float):
        return f"{value:,.3f}".rstrip("0").rstrip(".")
    return str(value)


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())

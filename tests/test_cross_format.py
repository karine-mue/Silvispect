"""Cross-format equivalence of ``.asc`` and ``.tif`` (SPEC-GIS.md section 7.2).

The same raster read from an ``.asc`` and from a GeoTIFF must be the same
grid, field for field, and must give the same analysis results.  The GeoTIFF
fixtures are assembled in the tests with ``tiff_builder`` (nothing is
committed); fixture (i) is a small hand-made grid of exactly representable
values, fixture (ii) is the repository's ``data/plot-a1-chm.asc``.  Every
``.asc`` written here uses ``precision=None``.
"""

from __future__ import annotations

import dataclasses
import io
import json
import struct
from pathlib import Path

import pytest

import tiff_builder
from silvispect.cli import EXIT_OK, main
from silvispect.grid import Grid

PLOT_ASC = Path(__file__).resolve().parent.parent / "data" / "plot-a1-chm.asc"
EPSG = 6677


# --- helpers ------------------------------------------------------------------


def _f32_exact(value: float) -> bool:
    try:
        return struct.unpack("<f", struct.pack("<f", value))[0] == value
    except OverflowError:
        return False


def _write_tif(path: Path, grid: Grid, *, epsg: int | None = None) -> Path:
    """Write ``grid`` as a GeoTIFF using only ``tiff_builder`` (32-bit if exact, else 64-bit)."""
    nodata = grid.nodata_value
    exact = _f32_exact(nodata) and all(v is None or _f32_exact(v) for v in grid.values)
    rows = [
        [nodata if v is None else v for v in grid.values[r * grid.ncols : (r + 1) * grid.ncols]]
        for r in range(grid.nrows)
    ]
    y_ul = grid.yllcorner + grid.nrows * grid.cellsize
    tags = [
        tiff_builder.scale_tag(grid.cellsize, grid.cellsize),
        tiff_builder.tiepoint_tag([0, 0, 0, grid.xllcorner, y_ul, 0]),
        tiff_builder.nodata_tag(repr(float(nodata))),
    ]
    if epsg is not None:
        tags.append(
            tiff_builder.geokey_tag((1024, 0, 1, 1), (1025, 0, 1, 1), (3072, 0, 1, epsg))
        )
    data = tiff_builder.build_tiff(rows, bits=32 if exact else 64, tags=tags)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(data)
    return path


def _save(path: Path, grid: Grid, *, epsg: int | None = None) -> Path:
    if path.suffix == ".tif":
        return _write_tif(path, grid, epsg=epsg)
    path.parent.mkdir(parents=True, exist_ok=True)
    grid.write(path, precision=None)
    return path


def _small_grid() -> Grid:
    """Fixture (i): quarter-multiples, two nodata cells, exact ``y_ul`` arithmetic."""
    values = [0.25, 1.5, 12.75, None, 3.0, 0.0, None, 7.25, 20.5, 2.0, 0.5, 9.75]
    return Grid(4, 3, 1000.0, 2000.0, 0.5, -9999.0, values)


def _plot_grid() -> Grid:
    """Fixture (ii): the repository CHM."""
    return Grid.read(PLOT_ASC)


def _dsm_dtm() -> tuple[Grid, Grid]:
    ground = [100.0 + 0.25 * i for i in range(12)]
    heights = [0.0, 12.5, 3.0, 0.25, 8.75, 20.0, 1.5, 0.0, 15.25, 6.0, 2.25, 9.5]
    dtm_values = list(ground)
    dsm_values = [g + h for g, h in zip(ground, heights, strict=True)]
    dtm_values[5] = None
    dsm_values[5] = None
    dtm = Grid(4, 3, 500.0, 1000.0, 1.0, -9999.0, dtm_values)
    dsm = Grid(4, 3, 500.0, 1000.0, 1.0, -9999.0, dsm_values)
    return dsm, dtm


def _run(argv: list[str]) -> tuple[int, str]:
    buffer = io.StringIO()
    code = main(argv, out=buffer)
    return code, buffer.getvalue()


def _scrub(node, replacements: list[str]):
    if isinstance(node, dict):
        return {key: _scrub(value, replacements) for key, value in node.items()}
    if isinstance(node, list):
        return [_scrub(value, replacements) for value in node]
    if isinstance(node, str):
        for text in replacements:
            node = node.replace(text, "<path>")
    return node


def _assert_same_fields(a: Grid, b: Grid, *, exact_yll: bool) -> None:
    assert a.ncols == b.ncols
    assert a.nrows == b.nrows
    assert a.cellsize == b.cellsize
    assert a.xllcorner == b.xllcorner
    if exact_yll:
        assert a.yllcorner == b.yllcorner
    else:
        assert abs(a.yllcorner - b.yllcorner) <= 1e-6 * a.cellsize
    assert a.nodata_value == b.nodata_value
    assert a.values == b.values


# --- field-wise equality: fixture (i) -----------------------------------------


def test_small_grid_asc_and_tif_agree_field_by_field(tmp_path):
    grid = _small_grid()
    a = Grid.read(_save(tmp_path / "g.asc", grid))
    b = Grid.read(_save(tmp_path / "g.tif", grid))
    _assert_same_fields(a, b, exact_yll=True)


def test_small_grid_without_geokeys_is_fully_equal(tmp_path):
    grid = _small_grid()
    a = Grid.read(_save(tmp_path / "g.asc", grid))
    b = Grid.read(_save(tmp_path / "g.tif", grid))
    assert b.crs is None
    assert a == b


def test_small_grid_tif_values_keep_nodata_positions(tmp_path):
    b = Grid.read(_write_tif(tmp_path / "g.tif", _small_grid()))
    assert [i for i, v in enumerate(b.values) if v is None] == [3, 6]


# --- field-wise equality: fixture (ii) ----------------------------------------


def test_plot_asc_and_tif_agree_field_by_field(tmp_path):
    a = _plot_grid()
    b = Grid.read(_write_tif(tmp_path / "plot.tif", a))
    _assert_same_fields(a, b, exact_yll=False)


def test_plot_asc_and_tif_are_aligned(tmp_path):
    a = _plot_grid()
    b = Grid.read(_write_tif(tmp_path / "plot.tif", a))
    a.assert_aligned(b)
    b.assert_aligned(a)


def test_plot_with_epsg_keeps_crs_only_on_the_tif(tmp_path):
    a = _plot_grid()
    b = Grid.read(_write_tif(tmp_path / "plot.tif", a, epsg=EPSG))
    assert a.crs is None
    assert b.crs == "EPSG:6677"
    _assert_same_fields(a, b, exact_yll=False)


# --- CRS in the equivalence ----------------------------------------------------


def test_epsg_tif_has_crs_and_asc_has_none(tmp_path):
    grid = _small_grid()
    a = Grid.read(_save(tmp_path / "g.asc", grid))
    b = Grid.read(_write_tif(tmp_path / "g.tif", grid, epsg=EPSG))
    assert a.crs is None
    assert b.crs == "EPSG:6677"


def test_epsg_tif_differs_from_asc_only_by_crs(tmp_path):
    grid = _small_grid()
    a = Grid.read(_save(tmp_path / "g.asc", grid))
    b = Grid.read(_write_tif(tmp_path / "g.tif", grid, epsg=EPSG))
    assert a != b
    assert dataclasses.replace(a, crs="EPSG:6677") == b


# --- round trips ---------------------------------------------------------------


def test_round_trip_through_tif_keeps_every_field_including_crs(tmp_path):
    values = [0.1, 0.2, None, 1 / 3, 12.3456789, 0.25]
    grid = Grid(3, 2, 10.0, 20.0, 0.5, -9999.0, values, crs="EPSG:6677")
    path = tmp_path / "rt.tif"
    grid.write(path)
    assert Grid.read(path) == grid


def test_round_trip_through_tif_is_exact_for_float32_exact_values(tmp_path):
    grid = _small_grid()
    path = tmp_path / "rt.tif"
    grid.write(path)
    assert Grid.read(path) == grid


def test_round_trip_through_tif_with_awkward_geometry_and_nodata(tmp_path):
    values = [0.1, None, 7.7, 1e-7, 3.3, 123456.789]
    grid = Grid(3, 2, 1000.7, 2000.123, 0.3, -9999.1234, values, crs="EPSG:6678")
    path = tmp_path / "rt.tif"
    grid.write(path)
    back = Grid.read(path)
    assert back.values == values
    assert back.nodata_value == -9999.1234
    assert back.xllcorner == grid.xllcorner
    assert back.cellsize == grid.cellsize
    assert abs(back.yllcorner - grid.yllcorner) <= 1e-6 * grid.cellsize
    assert back.crs == "EPSG:6678"


def test_round_trip_through_tif_of_float32_inexact_value_has_zero_tolerance(tmp_path):
    grid = Grid(2, 1, 0.0, 0.0, 1.0, values=[0.1, 16777217.0])
    path = tmp_path / "rt.tif"
    grid.write(path)
    assert Grid.read(path).values == [0.1, 16777217.0]


def test_round_trip_through_asc_keeps_data_and_drops_crs(tmp_path):
    grid = Grid(3, 2, 10.0, 20.0, 0.5, -9999.0, [0.1, 0.2, None, 1 / 3, 4.5, 6.0], crs="EPSG:6677")
    path = tmp_path / "rt.asc"
    grid.write(path, precision=None)
    back = Grid.read(path)
    assert back.crs is None
    assert back.values == grid.values
    assert back.nodata_value == grid.nodata_value
    assert (back.ncols, back.nrows) == (3, 2)
    assert (back.xllcorner, back.yllcorner, back.cellsize) == (10.0, 20.0, 0.5)


def test_tif_and_tiff_outputs_are_byte_identical(tmp_path):
    grid = Grid(3, 2, 10.0, 20.0, 0.5, -9999.0, [0.1, 0.2, None, 1 / 3, 4.5, 6.0], crs="EPSG:6677")
    grid.write(tmp_path / "a.tif")
    grid.write(tmp_path / "a.tiff")
    assert (tmp_path / "a.tif").read_bytes() == (tmp_path / "a.tiff").read_bytes()


# --- analysis results ------------------------------------------------------------


@pytest.mark.parametrize("make", [_small_grid, _plot_grid], ids=["small", "plot"])
def test_stats_agree_between_formats(tmp_path, make):
    grid = make()
    a = Grid.read(_save(tmp_path / "g.asc", grid))
    b = Grid.read(_write_tif(tmp_path / "g.tif", grid))
    assert a.stats() == b.stats()


@pytest.mark.parametrize("make", [_small_grid, _plot_grid], ids=["small", "plot"])
def test_histogram_agrees_between_formats(tmp_path, make):
    grid = make()
    a = Grid.read(_save(tmp_path / "g.asc", grid))
    b = Grid.read(_write_tif(tmp_path / "g.tif", grid))
    assert a.histogram(bins=7) == b.histogram(bins=7)


# --- CLI: the same report from either format ------------------------------------

_COMMANDS = {
    "inspect": lambda p: ["inspect", "--chm", p, "--format", "json"],
    "gaps": lambda p: ["gaps", p, "--json"],
    "detect": lambda p: ["detect", p, "--json"],
}


@pytest.mark.parametrize("command", sorted(_COMMANDS))
def test_cli_json_output_is_identical_for_asc_and_tif(tmp_path, command):
    grid = _plot_grid()
    asc = _save(tmp_path / "asc" / "chm.asc", grid)
    tif = _write_tif(tmp_path / "tif" / "chm.tif", grid)
    code_a, text_a = _run(_COMMANDS[command](str(asc)))
    code_b, text_b = _run(_COMMANDS[command](str(tif)))
    assert code_a == EXIT_OK
    assert code_b == EXIT_OK
    names = [str(asc), str(tif), asc.name, tif.name]
    assert _scrub(json.loads(text_a), names) == _scrub(json.loads(text_b), names)


def test_cli_inspect_accepts_a_tiff_extension(tmp_path):
    tiff = _write_tif(tmp_path / "chm.tiff", _plot_grid())
    code, text = _run(_COMMANDS["inspect"](str(tiff)))
    assert code == EXIT_OK
    assert json.loads(text)


# --- CLI: mixed inputs -------------------------------------------------------------


@pytest.mark.parametrize(
    ("dsm_name", "dtm_name"),
    [("dsm.asc", "dtm.tif"), ("dsm.tif", "dtm.asc"), ("dsm.tif", "dtm.tif")],
)
def test_chm_with_mixed_formats_matches_the_asc_only_result(tmp_path, dsm_name, dtm_name):
    dsm, dtm = _dsm_dtm()
    _save(tmp_path / "dsm.asc", dsm)
    _save(tmp_path / "dtm.asc", dtm)
    reference = tmp_path / "ref.asc"
    code, _ = _run(["chm", "--dsm", str(tmp_path / "dsm.asc"), "--dtm", str(tmp_path / "dtm.asc"), "-o", str(reference)])
    assert code == EXIT_OK

    dsm_path = _save(tmp_path / "mix" / dsm_name, dsm)
    dtm_path = _save(tmp_path / "mix" / dtm_name, dtm)
    out = tmp_path / "out.tif"
    code, _ = _run(["chm", "--dsm", str(dsm_path), "--dtm", str(dtm_path), "-o", str(out)])
    assert code == EXIT_OK

    expected = Grid.read(reference)
    got = Grid.read(out)
    assert got.values == expected.values
    assert got.nodata_value == expected.nodata_value
    got.assert_aligned(expected)


# --- CLI: CRS checks -----------------------------------------------------------------


def test_chm_with_differing_crs_exits_with_error_and_writes_nothing(tmp_path):
    dsm, dtm = _dsm_dtm()
    dsm_path = _write_tif(tmp_path / "dsm.tif", dsm, epsg=6677)
    dtm_path = _write_tif(tmp_path / "dtm.tif", dtm, epsg=6678)
    out = tmp_path / "out.tif"
    code, _ = _run(["chm", "--dsm", str(dsm_path), "--dtm", str(dtm_path), "-o", str(out)])
    assert code == 2
    assert not out.exists()


@pytest.mark.parametrize(
    ("dsm_epsg", "dtm_epsg"),
    [(6677, None), (None, 6677), (6677, 6677), (None, None)],
    ids=["dsm-only", "dtm-only", "same", "neither"],
)
def test_chm_with_unknown_or_equal_crs_succeeds(tmp_path, dsm_epsg, dtm_epsg):
    dsm, dtm = _dsm_dtm()
    dsm_path = _write_tif(tmp_path / "dsm.tif", dsm, epsg=dsm_epsg)
    dtm_path = _write_tif(tmp_path / "dtm.tif", dtm, epsg=dtm_epsg)
    out = tmp_path / "out.tif"
    code, _ = _run(["chm", "--dsm", str(dsm_path), "--dtm", str(dtm_path), "-o", str(out)])
    assert code == EXIT_OK
    assert out.is_file()


def test_chm_output_keeps_the_dsm_crs(tmp_path):
    # SPEC 5.4: the CHM derived from an EPSG-tagged DSM carries its CRS.
    dsm, dtm = _dsm_dtm()
    dsm_path = _write_tif(tmp_path / "dsm.tif", dsm, epsg=6677)
    dtm_path = _write_tif(tmp_path / "dtm.tif", dtm)
    out = tmp_path / "out.tif"
    code, _ = _run(["chm", "--dsm", str(dsm_path), "--dtm", str(dtm_path), "-o", str(out)])
    assert code == EXIT_OK
    assert Grid.read(out).crs == "EPSG:6677"

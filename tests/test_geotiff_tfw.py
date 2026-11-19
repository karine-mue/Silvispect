"""Tests for the ``.tfw`` world-file fallback (SPEC-GIS.md section 3).

A TIFF that carries none of ModelPixelScale / ModelTiepoint /
ModelTransformation falls back to the ``.tfw`` sidecar next to it.  GeoTIFF
tags always win and then the sidecar is not even opened.  A world file locates
the *centre* of the upper-left pixel, so the reader converts to the lower-left
corner of the lower-left cell (3.5).  Fixture TIFFs come from the independent
``tiff_builder`` so the reader shares no code with the files it is fed.
"""

from __future__ import annotations

import pytest
import tiff_builder as tb

from silvispect.grid import Grid, GridError

# Appendix B: A=0.5, D=0, B=0, E=-0.5, C=612500.25, F=6642049.75, nrows=100.
NROWS = 100
NCOLS = 2


def _rows(nrows: int = NROWS, ncols: int = NCOLS) -> list[list[float]]:
    return [[float(c + 1) for c in range(ncols)] for _ in range(nrows)]


def _tif(nrows: int = NROWS, **kwargs) -> bytes:
    return tb.build_tiff(_rows(nrows), **kwargs)


def _tfw(a=0.5, d=0.0, b=0.0, e=-0.5, c=612500.25, f=6642049.75) -> str:
    # File order is A, D, B, E, C, F (3.4), one value per line.
    return "".join(f"{value!r}\n" for value in (a, d, b, e, c, f))


def _put(tmp_path, tfw_text: str | bytes | None, *, stem: str = "x", tif: bytes | None = None):
    """Write ``<stem>.tif`` and, if given, ``<stem>.tfw``; return the tif path."""
    path = tmp_path / f"{stem}.tif"
    path.write_bytes(_tif() if tif is None else tif)
    if tfw_text is not None:
        side = tmp_path / f"{stem}.tfw"
        if isinstance(tfw_text, bytes):
            side.write_bytes(tfw_text)
        else:
            side.write_text(tfw_text, encoding="utf-8")
    return path


def _case_sensitive(directory) -> bool:
    probe = directory / "probe.lower"
    probe.write_bytes(b"")
    try:
        return not (directory / "PROBE.LOWER").exists()
    finally:
        probe.unlink()


# --- 3.4 / 3.5 / appendix B: conversion -----------------------------------


def test_tfw_appendix_b_example(tmp_path):
    grid = Grid.read(_put(tmp_path, _tfw()))
    assert grid.xllcorner == 612500.0
    assert grid.yllcorner == 6642000.0
    assert grid.cellsize == 0.5
    assert (grid.ncols, grid.nrows) == (NCOLS, NROWS)


def test_tfw_second_numeric_example_of_the_formula(tmp_path):
    # xll = C - A/2 = 11 - 1 = 10;  yll = F - (nrows - 0.5) * A = 19 - 2.5*2 = 14
    path = _put(tmp_path, _tfw(a=2.0, e=-2.0, c=11.0, f=19.0), tif=_tif(3))
    grid = Grid.read(path)
    assert (grid.xllcorner, grid.yllcorner, grid.cellsize) == (10.0, 14.0, 2.0)
    assert grid.nrows == 3


def test_tfw_leaves_cell_values_untouched(tmp_path):
    grid = Grid.read(_put(tmp_path, _tfw()))
    assert grid.values == [1.0, 2.0] * NROWS


def test_tfw_with_crlf_line_endings(tmp_path):
    # 3.4: tokens are split on any whitespace.
    grid = Grid.read(_put(tmp_path, _tfw().replace("\n", "\r\n")))
    assert grid.xllcorner == 612500.0
    assert grid.yllcorner == 6642000.0


def test_tfw_with_tiff_extension_uses_same_stem(tmp_path):
    path = tmp_path / "plot.a1.tiff"
    path.write_bytes(_tif())
    (tmp_path / "plot.a1.tfw").write_text(_tfw(), encoding="utf-8")
    grid = Grid.read(path)
    assert grid.xllcorner == 612500.0


def test_tfw_negative_zero_rotation_terms_are_zero(tmp_path):
    # 3.5: "-0.0 はゼロ".
    grid = Grid.read(_put(tmp_path, _tfw(d=-0.0, b=-0.0)))
    assert grid.cellsize == 0.5


def test_tfw_with_utf8_bom_is_accepted(tmp_path):
    # 3.4: decoded as utf-8-sig.
    grid = Grid.read(_put(tmp_path, b"\xef\xbb\xbf" + _tfw().encode("ascii")))
    assert grid.xllcorner == 612500.0
    assert grid.yllcorner == 6642000.0


# --- 3.2: GeoTIFF tags take priority; the sidecar is not opened ------------

_TAG_SETS = {
    "scale+tiepoint": [
        tb.scale_tag(0.5, 0.5),
        tb.tiepoint_tag((0, 0, 0, 612500.0, 6642050.0, 0)),
    ],
    "transformation": [
        tb.transform_tag(
            [0.5, 0, 0, 612500.0, 0, -0.5, 0, 6642050.0, 0, 0, 1, 0, 0, 0, 0, 1]
        )
    ],
}


@pytest.mark.parametrize("tag_set", sorted(_TAG_SETS))
@pytest.mark.parametrize(
    "bad_sidecar",
    [
        "not a world file at all",
        _tfw(b=0.1),  # rotated
        _tfw(e=0.5),  # south-up
        "1 2 3\n",  # wrong token count
        b"\xff\xfe\x00bad",  # undecodable
    ],
    ids=["garbage", "rotated", "south-up", "short", "undecodable"],
)
def test_tags_win_over_an_invalid_tfw(tmp_path, tag_set, bad_sidecar):
    path = _put(tmp_path, bad_sidecar, tif=_tif(tags=_TAG_SETS[tag_set]))
    grid = Grid.read(path)
    assert grid.xllcorner == 612500.0
    assert grid.yllcorner == 6642000.0
    assert grid.cellsize == 0.5


def test_tags_win_over_a_valid_but_different_tfw(tmp_path):
    other = _tfw(a=2.0, e=-2.0, c=1.0, f=1.0)
    path = _put(tmp_path, other, tif=_tif(tags=_TAG_SETS["scale+tiepoint"]))
    grid = Grid.read(path)
    assert (grid.xllcorner, grid.yllcorner, grid.cellsize) == (612500.0, 6642000.0, 0.5)


def test_partial_tags_do_not_fall_back_to_tfw(tmp_path):
    # 2.6 / 3.1: a lone PixelScale is corrupt, not "no georeferencing".
    path = _put(tmp_path, _tfw(), tif=_tif(tags=[tb.scale_tag(0.5, 0.5)]))
    with pytest.raises(GridError, match="corrupt TIFF"):
        Grid.read(path)


# --- 3.3: sidecar lookup ---------------------------------------------------


def test_uppercase_tfw_sidecar_is_found(tmp_path):
    path = tmp_path / "x.tif"
    path.write_bytes(_tif())
    (tmp_path / "x.TFW").write_text(_tfw(), encoding="utf-8")
    grid = Grid.read(path)
    assert grid.xllcorner == 612500.0
    assert grid.yllcorner == 6642000.0


def test_lowercase_tfw_preferred_over_uppercase(tmp_path):
    if not _case_sensitive(tmp_path):
        pytest.skip("case-insensitive file system cannot hold both x.tfw and x.TFW")
    path = tmp_path / "x.tif"
    path.write_bytes(_tif())
    (tmp_path / "x.tfw").write_text(_tfw(), encoding="utf-8")
    (tmp_path / "x.TFW").write_text(_tfw(c=1000.25, f=2000.25), encoding="utf-8")
    grid = Grid.read(path)
    assert grid.xllcorner == 612500.0


def test_invalid_lowercase_tfw_is_not_skipped_for_uppercase(tmp_path):
    # The lowercase file exists, so it is the one used (and found corrupt).
    if not _case_sensitive(tmp_path):
        pytest.skip("case-insensitive file system cannot hold both x.tfw and x.TFW")
    path = tmp_path / "x.tif"
    path.write_bytes(_tif())
    (tmp_path / "x.tfw").write_text("garbage", encoding="utf-8")
    (tmp_path / "x.TFW").write_text(_tfw(), encoding="utf-8")
    with pytest.raises(GridError, match="corrupt world file"):
        Grid.read(path)


# --- 3.5: rejections -------------------------------------------------------


@pytest.mark.parametrize("field", ["b", "d"])
def test_tfw_rotation_terms_are_rejected(tmp_path, field):
    path = _put(tmp_path, _tfw(**{field: 1e-12}))  # no tolerance (3.5)
    with pytest.raises(GridError, match="unsupported georeferencing: rotated"):
        Grid.read(path)


def test_tfw_both_rotation_terms_are_rejected(tmp_path):
    path = _put(tmp_path, _tfw(b=0.25, d=-0.25))
    with pytest.raises(GridError, match="unsupported georeferencing: rotated"):
        Grid.read(path)


def test_tfw_south_up_is_rejected(tmp_path):
    path = _put(tmp_path, _tfw(e=0.5))
    with pytest.raises(GridError, match="unsupported georeferencing: non-north-up"):
        Grid.read(path)


def test_tfw_zero_e_is_rejected(tmp_path):
    path = _put(tmp_path, _tfw(e=0.0))
    with pytest.raises(GridError, match="unsupported georeferencing: non-north-up"):
        Grid.read(path)


def test_tfw_non_positive_a_is_rejected(tmp_path):
    path = _put(tmp_path, _tfw(a=-0.5, e=-0.5))
    with pytest.raises(GridError, match="unsupported georeferencing: non-north-up"):
        Grid.read(path)


def test_tfw_non_square_cells_are_rejected(tmp_path):
    path = _put(tmp_path, _tfw(a=0.5, e=-0.6))
    with pytest.raises(GridError, match="unsupported georeferencing: non-square"):
        Grid.read(path)


def test_tfw_nearly_square_cells_are_accepted(tmp_path):
    # rel_tol=1e-9 (3.5): a difference far below it is not "non-square".
    grid = Grid.read(_put(tmp_path, _tfw(a=0.5, e=-0.5 * (1 + 1e-12))))
    assert grid.cellsize == 0.5


# --- 3.4 / 6.2 #28: corrupt world files -----------------------------------


@pytest.mark.parametrize(
    "text",
    [
        "",
        "0.5\n0\n0\n-0.5\n612500.25\n",  # five tokens
        "0.5\n0\n0\n-0.5\n612500.25\n6642049.75\n1.0\n",  # seven tokens
        "0.5 0 0 -0.5 612500.25 6642049.75 0",
        "0.5\n0\n0\n-0.5\nabc\n6642049.75\n",  # non-numeric
        "0.5\n0\n0\n-0.5\n612500.25\n6642049.75x\n",
        "0.5\n0\n0\n-0.5\nnan\n6642049.75\n",  # non-finite
        "0.5\n0\n0\n-0.5\n612500.25\ninf\n",
        "-inf\n0\n0\n-0.5\n612500.25\n6642049.75\n",
    ],
    ids=[
        "empty",
        "five-tokens",
        "seven-tokens",
        "seven-tokens-one-line",
        "non-numeric",
        "trailing-garbage",
        "nan",
        "inf",
        "minus-inf",
    ],
)
def test_corrupt_tfw_is_rejected(tmp_path, text):
    with pytest.raises(GridError, match="corrupt world file"):
        Grid.read(_put(tmp_path, text))


def test_undecodable_tfw_is_corrupt(tmp_path):
    with pytest.raises(GridError, match="corrupt world file"):
        Grid.read(_put(tmp_path, b"0.5\n0\n0\n-0.5\n\xff\xfe\n6642049.75\n"))


# --- 3.1 / 2.7: GeoKeys alone still fall back to .tfw, keeping the CRS -----

_GEOKEYS_6677 = tb.geokey_tag((1024, 0, 1, 1), (1025, 0, 1, 1), (3072, 0, 1, 6677))


def test_geokey_only_tiff_falls_back_to_tfw_and_keeps_crs(tmp_path):
    path = _put(tmp_path, _tfw(), tif=_tif(tags=[_GEOKEYS_6677]))
    grid = Grid.read(path)
    assert grid.crs == "EPSG:6677"
    assert grid.xllcorner == 612500.0
    assert grid.yllcorner == 6642000.0
    assert grid.cellsize == 0.5


def test_tfw_route_ignores_pixel_is_point_geokey(tmp_path):
    # 2.7: a world file always refers to pixel centres; GTRasterType is ignored.
    point_keys = tb.geokey_tag((1024, 0, 1, 1), (1025, 0, 1, 2), (3072, 0, 1, 6677))
    path = _put(tmp_path, _tfw(), tif=_tif(tags=[point_keys]))
    grid = Grid.read(path)
    assert grid.xllcorner == 612500.0
    assert grid.yllcorner == 6642000.0


def test_tfw_without_geokeys_gives_no_crs(tmp_path):
    assert Grid.read(_put(tmp_path, _tfw())).crs is None


def test_geographic_geokey_is_still_rejected_with_a_tfw(tmp_path):
    keys = tb.geokey_tag((1024, 0, 1, 2), (2048, 0, 1, 4326))
    path = _put(tmp_path, _tfw(), tif=_tif(tags=[keys]))
    with pytest.raises(GridError, match="unsupported CRS"):
        Grid.read(path)


# --- 3.3 / 6.2 #20: nothing to georeference with ---------------------------


def test_no_tags_and_no_tfw_is_missing_georeferencing(tmp_path):
    with pytest.raises(GridError, match="missing georeferencing"):
        Grid.read(_put(tmp_path, None))


def test_geokey_only_without_tfw_is_missing_georeferencing(tmp_path):
    path = _put(tmp_path, None, tif=_tif(tags=[_GEOKEYS_6677]))
    with pytest.raises(GridError, match="missing georeferencing"):
        Grid.read(path)


def test_tfw_of_another_stem_is_not_used(tmp_path):
    path = _put(tmp_path, None)
    (tmp_path / "y.tfw").write_text(_tfw(), encoding="utf-8")
    with pytest.raises(GridError, match="missing georeferencing"):
        Grid.read(path)


def test_directory_named_like_the_sidecar_is_not_a_sidecar(tmp_path):
    # 3.3: the lookup is is_file(), not exists().
    path = _put(tmp_path, None)
    (tmp_path / "x.tfw").mkdir()
    with pytest.raises(GridError, match="missing georeferencing"):
        Grid.read(path)

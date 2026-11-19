"""Tests for ``Grid.crs`` (SPEC-GIS.md section 5).

The field is an optional, keyword-only ``"EPSG:<code>"`` string.  These tests
pin its validation, its propagation through the derived-grid constructors, the
CRS rules of ``combine`` and ``assert_aligned``, the ``repr`` and equality
contract, positional-construction compatibility, and the fact that the
``.asc`` text format neither stores nor reads a CRS.

Tests that fix current behaviour (positional construction, the ``repr`` of a
grid without a CRS, an ``.asc`` round trip) are expected to pass today; the
rest fail until ``crs`` exists.
"""

from __future__ import annotations

import dataclasses

import pytest

from silvispect.grid import Grid, GridError

EPSG_6677 = "EPSG:6677"
EPSG_6678 = "EPSG:6678"


def _grid(**kwargs) -> Grid:
    return Grid(2, 1, 0.0, 0.0, 1.0, values=[1.0, 2.0], **kwargs)


# --- 5.1: default, type, validation ---------------------------------------


def test_crs_defaults_to_none():
    assert _grid().crs is None


@pytest.mark.parametrize("code", ["EPSG:1", "EPSG:4", "EPSG:6677", "EPSG:32766"])
def test_valid_epsg_strings_are_accepted(code):
    assert _grid(crs=code).crs == code


@pytest.mark.parametrize(
    "bad",
    [
        "epsg:6677",  # lower case is not normalised
        "Epsg:6677",
        "EPSG:0",
        "EPSG:06677",  # leading zero
        "EPSG:32767",  # out of range
        "EPSG:99999",
        "EPSG:-1",
        "EPSG:+6677",
        "EPSG:",
        "EPSG: 6677",
        " EPSG:6677",
        "EPSG:6677 ",
        "EPSG:6677\n",
        "EPSG:66 77",
        "EPSG:6677.0",
        "6677",
        "",
        "urn:ogc:def:crs:EPSG::6677",
    ],
)
def test_malformed_crs_string_is_rejected(bad):
    with pytest.raises(GridError, match=r"crs must be None or 'EPSG:<code>'"):
        _grid(crs=bad)


@pytest.mark.parametrize("bad", [6677, True, False, 0, 6677.0, b"EPSG:6677", ["EPSG:6677"]])
def test_non_string_crs_is_rejected(bad):
    with pytest.raises(GridError, match=r"crs must be None or 'EPSG:<code>'"):
        _grid(crs=bad)


def test_from_rows_validates_crs():
    with pytest.raises(GridError, match="crs must be None"):
        Grid.from_rows([[1.0]], crs="epsg:6677")


def test_filled_validates_crs():
    with pytest.raises(GridError, match="crs must be None"):
        Grid.filled(1, 1, 1.0, crs="EPSG:0")


# --- 5.1: positional-construction compatibility ----------------------------


def test_seven_positional_arguments_keep_their_order_and_meaning():
    grid = Grid(2, 1, 10.0, 20.0, 0.5, -1.0, [1.0, None])
    assert grid.ncols == 2
    assert grid.nrows == 1
    assert grid.xllcorner == 10.0
    assert grid.yllcorner == 20.0
    assert grid.cellsize == 0.5
    assert grid.nodata_value == -1.0
    assert grid.values == [1.0, None]


def test_eighth_positional_argument_is_a_type_error():
    with pytest.raises(TypeError):
        Grid(2, 1, 10.0, 20.0, 0.5, -1.0, [1.0, None], EPSG_6677)


def test_crs_is_the_last_field_and_keyword_only():
    fields = dataclasses.fields(Grid)
    assert [f.name for f in fields[:7]] == [
        "ncols",
        "nrows",
        "xllcorner",
        "yllcorner",
        "cellsize",
        "nodata_value",
        "values",
    ]
    assert fields[-1].name == "crs"
    assert fields[-1].kw_only is True
    assert fields[-1].default is None


# --- 5.4: propagation ------------------------------------------------------


def test_like_carries_crs():
    assert _grid(crs=EPSG_6677).like().crs == EPSG_6677


def test_like_of_grid_without_crs_has_none():
    assert _grid().like().crs is None


def test_copy_carries_crs_and_values():
    original = _grid(crs=EPSG_6677)
    clone = original.copy()
    assert clone.crs == EPSG_6677
    assert clone.values == original.values
    assert clone == original


def test_copy_is_deep_for_values():
    original = _grid(crs=EPSG_6677)
    clone = original.copy()
    clone.values[0] = 99.0
    assert original.values[0] == 1.0


def test_from_rows_passes_crs_through():
    grid = Grid.from_rows([[1.0, 2.0]], cellsize=0.5, crs=EPSG_6677)
    assert grid.crs == EPSG_6677


def test_from_rows_default_crs_is_none():
    assert Grid.from_rows([[1.0, 2.0]]).crs is None


def test_filled_passes_crs_through():
    grid = Grid.filled(2, 2, 1.0, crs=EPSG_6678)
    assert grid.crs == EPSG_6678


def test_filled_default_crs_is_none():
    assert Grid.filled(2, 2, 1.0).crs is None


def test_crs_is_keyword_only_for_constructors_helpers():
    with pytest.raises(TypeError):
        Grid.from_rows([[1.0]], 1.0, 0.0, 0.0, -9999.0, EPSG_6677)


@pytest.mark.parametrize(
    "derive",
    [
        lambda g: g.map_values(lambda v: v + 1.0),
        lambda g: g.clip(0.0, 1.5),
        lambda g: g.replace_nodata(0.0),
        lambda g: g.smooth_mean(1),
        lambda g: g.smooth_median(1),
        lambda g: g.focal_max(1),
        lambda g: g.focal(max, 1),
    ],
    ids=["map_values", "clip", "replace_nodata", "smooth_mean", "smooth_median", "focal_max", "focal"],
)
def test_derived_grids_keep_crs(derive):
    assert derive(_grid(crs=EPSG_6677)).crs == EPSG_6677


# --- 5.4: combine ----------------------------------------------------------


def _add(a: float, b: float) -> float:
    return a + b


def test_combine_of_equal_crs_keeps_it():
    assert _grid(crs=EPSG_6677).combine(_grid(crs=EPSG_6677), _add).crs == EPSG_6677


def test_combine_takes_other_crs_when_self_has_none():
    assert _grid().combine(_grid(crs=EPSG_6677), _add).crs == EPSG_6677


def test_combine_keeps_self_crs_when_other_has_none():
    assert _grid(crs=EPSG_6677).combine(_grid(), _add).crs == EPSG_6677


def test_combine_of_two_unknown_crs_is_none():
    assert _grid().combine(_grid(), _add).crs is None


def test_combine_refuses_differing_crs():
    with pytest.raises(GridError, match="grids differ in CRS"):
        _grid(crs=EPSG_6677).combine(_grid(crs=EPSG_6678), _add)


def test_combine_values_are_unchanged_by_crs():
    out = _grid(crs=EPSG_6677).combine(_grid(), _add)
    assert out.values == [2.0, 4.0]


# --- 5.5: assert_aligned ---------------------------------------------------


def test_assert_aligned_refuses_two_different_crs():
    with pytest.raises(GridError, match="grids differ in CRS"):
        _grid(crs=EPSG_6677).assert_aligned(_grid(crs=EPSG_6678))


def test_assert_aligned_passes_for_equal_crs():
    _grid(crs=EPSG_6677).assert_aligned(_grid(crs=EPSG_6677))


def test_assert_aligned_passes_when_only_self_has_crs():
    _grid(crs=EPSG_6677).assert_aligned(_grid())


def test_assert_aligned_passes_when_only_other_has_crs():
    _grid().assert_aligned(_grid(crs=EPSG_6677))


def test_assert_aligned_passes_when_neither_has_crs():
    _grid().assert_aligned(_grid())


def test_assert_aligned_crs_check_is_symmetric():
    a, b = _grid(crs=EPSG_6677), _grid()
    a.assert_aligned(b)
    b.assert_aligned(a)
    c = _grid(crs=EPSG_6678)
    with pytest.raises(GridError):
        a.assert_aligned(c)
    with pytest.raises(GridError):
        c.assert_aligned(a)


def test_shape_error_comes_before_crs_error():
    a = Grid(2, 1, 0.0, 0.0, 1.0, values=[1.0, 2.0], crs=EPSG_6677)
    b = Grid(1, 1, 0.0, 0.0, 1.0, values=[1.0], crs=EPSG_6678)
    with pytest.raises(GridError, match="grids differ in shape"):
        a.assert_aligned(b)


def test_cellsize_error_comes_before_crs_error():
    a = Grid(2, 1, 0.0, 0.0, 1.0, values=[1.0, 2.0], crs=EPSG_6677)
    b = Grid(2, 1, 0.0, 0.0, 2.0, values=[1.0, 2.0], crs=EPSG_6678)
    with pytest.raises(GridError, match="grids differ in cellsize"):
        a.assert_aligned(b)


def test_origin_error_comes_before_crs_error():
    a = Grid(2, 1, 0.0, 0.0, 1.0, values=[1.0, 2.0], crs=EPSG_6677)
    b = Grid(2, 1, 5.0, 0.0, 1.0, values=[1.0, 2.0], crs=EPSG_6678)
    with pytest.raises(GridError, match="grids differ in origin"):
        a.assert_aligned(b)


# --- 5.2: repr and equality -------------------------------------------------


def test_repr_without_crs_is_the_historic_string():
    grid = Grid(2, 1, 10.0, 20.0, 0.5, values=[1.0, 2.0])
    assert repr(grid) == "Grid(nrows=1, ncols=2, cellsize=0.5, origin=(10.0, 20.0))"


def test_repr_with_crs_appends_it_before_the_closing_parenthesis():
    grid = Grid(2, 1, 10.0, 20.0, 0.5, values=[1.0, 2.0], crs=EPSG_6677)
    assert repr(grid) == (
        "Grid(nrows=1, ncols=2, cellsize=0.5, origin=(10.0, 20.0), crs='EPSG:6677')"
    )


def test_equality_includes_crs():
    assert _grid(crs=EPSG_6677) == _grid(crs=EPSG_6677)
    assert _grid(crs=EPSG_6677) != _grid(crs=EPSG_6678)


def test_none_crs_does_not_equal_an_epsg_crs():
    assert _grid() != _grid(crs=EPSG_6677)


def test_grids_without_crs_still_compare_equal():
    assert _grid() == _grid()


def test_replace_can_set_crs_for_comparison():
    assert dataclasses.replace(_grid(), crs=EPSG_6677) == _grid(crs=EPSG_6677)


# --- 5.3: the .asc text format ignores crs ----------------------------------


def test_parse_returns_none_crs():
    text = "ncols 2\nnrows 1\nxllcorner 0\nyllcorner 0\ncellsize 1\n1 2\n"
    assert Grid.parse(text).crs is None


def test_to_text_ignores_crs():
    with_crs = _grid(crs=EPSG_6677)
    without = _grid()
    assert with_crs.to_text() == without.to_text()
    assert with_crs.to_text(precision=None) == without.to_text(precision=None)


def test_text_round_trip_drops_crs():
    back = Grid.parse(_grid(crs=EPSG_6677).to_text())
    assert back.crs is None
    assert back.values == [1.0, 2.0]


def test_asc_write_drops_crs_and_writes_no_sidecar(tmp_path):
    path = tmp_path / "chm.asc"
    _grid(crs=EPSG_6677).write(path)
    assert [p.name for p in tmp_path.iterdir()] == ["chm.asc"]
    assert Grid.read(path).crs is None

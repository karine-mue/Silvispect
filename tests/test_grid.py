"""Tests for the ASCII grid raster."""

from __future__ import annotations

import decimal
import math
import sys

import pytest

from silvispect.grid import Grid, GridError, mean_of, rms_of, stdev_of, sum_of

MAX = sys.float_info.max

SAMPLE = """\
ncols 3
nrows 2
xllcorner 100.0
yllcorner 200.0
cellsize 2.0
NODATA_value -9999
1.0 2.0 3.0
4.0 -9999 6.0
"""


def test_parse_round_trip():
    grid = Grid.parse(SAMPLE)
    assert (grid.ncols, grid.nrows) == (3, 2)
    assert grid.cellsize == 2.0
    assert grid.get(0, 0) == 1.0
    assert grid.get(1, 1) is None
    reparsed = Grid.parse(grid.to_text())
    assert reparsed.values == grid.values
    assert reparsed.extent == grid.extent


def test_write_and_read(tmp_path):
    grid = Grid.parse(SAMPLE)
    path = grid.write(tmp_path / "nested" / "grid.asc")
    assert path.exists()
    assert Grid.read(path).values == grid.values


def test_header_is_case_insensitive_and_accepts_centers():
    text = SAMPLE.replace("xllcorner 100.0", "XLLCENTER 101.0").replace(
        "NODATA_value", "nodata_value"
    )
    grid = Grid.parse(text)
    assert grid.xllcorner == pytest.approx(100.0)


def test_parse_rejects_wrong_cell_count():
    with pytest.raises(GridError, match="expected 6 cell values"):
        Grid.parse(SAMPLE.replace("4.0 -9999 6.0", "4.0 -9999"))


def test_parse_rejects_missing_header():
    with pytest.raises(GridError, match="missing required header"):
        Grid.parse("nrows 2\ncellsize 1\n1 2\n3 4\n")


def test_parse_rejects_non_numeric_cell():
    with pytest.raises(GridError, match="non-numeric"):
        Grid.parse(SAMPLE.replace("2.0 3.0", "2.0 tree"))


def test_geometry():
    grid = Grid.parse(SAMPLE)
    # Row 0 is the northern row.
    assert grid.cell_center(0, 0) == pytest.approx((101.0, 203.0))
    assert grid.cell_center(1, 2) == pytest.approx((105.0, 201.0))
    assert grid.cell_of(101.0, 203.0) == (0, 0)
    assert grid.cell_of(105.0, 201.0) == (1, 2)
    assert grid.cell_of(0.0, 0.0) is None
    assert grid.extent == (100.0, 200.0, 106.0, 204.0)
    assert grid.area == pytest.approx(24.0)
    assert grid.area_ha == pytest.approx(0.0024)


def test_cell_center_round_trips_through_cell_of(tiny_grid):
    for row, col in tiny_grid.cells():
        x, y = tiny_grid.cell_center(row, col)
        assert tiny_grid.cell_of(x, y) == (row, col)


def test_index_bounds(tiny_grid):
    with pytest.raises(IndexError):
        tiny_grid.get(4, 0)
    with pytest.raises(IndexError):
        tiny_grid.cell_center(-1, 0)


def test_stats_skips_nodata():
    grid = Grid.parse(SAMPLE)
    stats = grid.stats()
    assert stats.count == 5
    assert stats.nodata_count == 1
    assert stats.minimum == 1.0
    assert stats.maximum == 6.0
    assert stats.mean == pytest.approx(3.2)
    assert stats.stdev == pytest.approx(math.sqrt(3.7), rel=1e-9)
    assert stats.as_dict()["count"] == 5


def test_stats_of_empty_grid():
    grid = Grid.filled(2, 2, None)
    stats = grid.stats()
    assert stats.count == 0
    assert stats.mean is None


def test_combine_propagates_nodata():
    grid = Grid.parse(SAMPLE)
    doubled = grid.combine(grid, lambda a, b: a + b)
    assert doubled.get(0, 0) == 2.0
    assert doubled.get(1, 1) is None


def test_combine_rejects_misaligned():
    grid = Grid.parse(SAMPLE)
    other = Grid.filled(2, 3, 1.0, cellsize=5.0)
    with pytest.raises(GridError, match="cellsize"):
        grid.combine(other, lambda a, b: a + b)
    with pytest.raises(GridError, match="shape"):
        grid.combine(Grid.filled(3, 3, 1.0, cellsize=2.0), lambda a, b: a + b)


def test_map_and_clip():
    grid = Grid.parse(SAMPLE)
    assert grid.map_values(lambda v: v * 2).get(0, 1) == 4.0
    clipped = grid.clip(minimum=2.0, maximum=5.0)
    assert clipped.get(0, 0) == 2.0
    assert clipped.get(1, 2) == 5.0
    assert clipped.get(1, 1) is None


def test_replace_nodata():
    grid = Grid.parse(SAMPLE).replace_nodata(0.0)
    assert grid.get(1, 1) == 0.0
    assert grid.stats().nodata_count == 0


def test_focal_mean_smooths_the_peak(tiny_grid):
    smoothed = tiny_grid.smooth_mean(1)
    assert smoothed.get(1, 1) < tiny_grid.get(1, 1)
    # The kernel is circular, so radius 1 covers the four edge neighbours only.
    assert smoothed.get(1, 1) == pytest.approx((5.0 + 0.0 + 4.0 + 0.0 + 4.0) / 5.0)


def test_focal_median_removes_a_spike():
    grid = Grid.from_rows([[1.0] * 5 for _ in range(5)])
    grid.set(2, 2, 99.0)
    assert grid.smooth_median(1).get(2, 2) == 1.0


def test_focal_max_and_window_shapes(tiny_grid):
    assert tiny_grid.focal_max(1).get(1, 0) == 5.0
    assert tiny_grid.focal_max(1).get(0, 0) == 0.0  # (1, 1) is outside a circular r=1
    square = set(tiny_grid.window(1, 1, 1))
    circular = set(tiny_grid.window(1, 1, 1, circular=True))
    assert circular < square
    assert (0, 0) in square and (0, 0) not in circular


def test_neighbors_connectivity(tiny_grid):
    assert len(list(tiny_grid.neighbors(1, 1, 4))) == 4
    assert len(list(tiny_grid.neighbors(1, 1, 8))) == 8
    assert len(list(tiny_grid.neighbors(0, 0, 8))) == 3
    with pytest.raises(GridError):
        list(tiny_grid.neighbors(0, 0, 6))


def test_histogram():
    grid = Grid.from_rows([[0.0, 1.0], [2.0, 3.0]])
    edges, counts = grid.histogram(bins=3, lo=0.0, hi=3.0)
    assert len(edges) == 4
    assert sum(counts) == 4


def test_from_rows_validates():
    with pytest.raises(GridError, match="same length"):
        Grid.from_rows([[1.0, 2.0], [3.0]])
    with pytest.raises(GridError, match="zero rows"):
        Grid.from_rows([])


def test_construction_validates():
    with pytest.raises(GridError, match="positive number"):
        Grid(0, 1, 0.0, 0.0, 1.0)
    with pytest.raises(GridError, match="cellsize"):
        Grid(1, 1, 0.0, 0.0, 0.0)
    with pytest.raises(GridError, match="expected 4 values"):
        Grid(2, 2, 0.0, 0.0, 1.0, values=[1.0])


def test_copy_is_independent(tiny_grid):
    clone = tiny_grid.copy()
    clone.set(0, 0, 42.0)
    assert tiny_grid.get(0, 0) == 0.0


def test_item_access(tiny_grid):
    tiny_grid[0, 0] = 9.0
    assert tiny_grid[0, 0] == 9.0


def test_valid_cells_matches_finite_values():
    grid = Grid.parse(SAMPLE)
    assert [value for _, _, value in grid.valid_cells()] == list(grid.finite_values())


def test_serialisation_trims_trailing_zeros():
    text = Grid.from_rows([[1.5, 2.0]]).to_text()
    assert "1.5 2" in text
    assert "NODATA_value -9999" in text


@pytest.mark.parametrize("origin", [0.0, 1.0, 1e3, 1e6, 1e9, 1e12])
def test_origin_tolerance_does_not_grow_with_coordinates(origin):
    """A half-cell shift is a different raster at every easting.

    The allowance for origins is a fraction of a cell, so it means the same
    thing on a projected grid a billion metres from the false origin as it does
    at zero.  A relative tolerance on the coordinates instead let the allowance
    grow with them, and two rasters half a cell apart at an easting of 1e9 were
    subtracted cell for cell into a quietly wrong canopy model.
    """
    reference = Grid.from_rows([[10.0]], cellsize=1.0, xllcorner=origin)
    for shift in (0.5, 0.1, 0.01):
        shifted = Grid.from_rows([[1.0]], cellsize=1.0, xllcorner=origin + shift)
        with pytest.raises(GridError, match="origin"):
            reference.assert_aligned(shifted)


@pytest.mark.parametrize("cellsize", [1e-9, 1e-3, 1.0, 10.0, 1e6])
def test_origin_tolerance_is_a_fraction_of_a_cell(cellsize):
    """Rounding noise is forgiven; a visible fraction of a cell is not."""
    reference = Grid.from_rows([[10.0]], cellsize=cellsize)
    noisy = Grid.from_rows([[1.0]], cellsize=cellsize, xllcorner=cellsize * 1e-9)
    reference.assert_aligned(noisy)  # must not raise
    displaced = Grid.from_rows([[1.0]], cellsize=cellsize, xllcorner=cellsize * 0.01)
    with pytest.raises(GridError, match="origin"):
        reference.assert_aligned(displaced)


def test_window_yields_exactly_the_in_bounds_offsets():
    """Clipping the walk must not change which cells come back, or their order."""
    import random

    def unclipped(grid, row, col, radius, circular):
        out = []
        for drow in range(-radius, radius + 1):
            for dcol in range(-radius, radius + 1):
                if circular and drow * drow + dcol * dcol > radius * radius:
                    continue
                if grid.in_bounds(row + drow, col + dcol):
                    out.append((row + drow, col + dcol))
        return out

    rng = random.Random(4242)
    for _ in range(400):
        grid = Grid.filled(rng.randint(1, 7), rng.randint(1, 7), 0.0)
        row, col = rng.randrange(grid.nrows), rng.randrange(grid.ncols)
        radius, circular = rng.randint(0, 9), rng.random() < 0.5
        assert list(grid.window(row, col, radius, circular=circular)) == unclipped(
            grid, row, col, radius, circular
        )


def test_window_cost_follows_the_raster_not_the_requested_radius():
    """A one-cell raster is one cell of work however wide the window asks to be.

    The detector sizes its search window from the height of the cell it is
    testing, so an accepted — if implausible — height buys an enormous radius.
    Walking the whole requested square made a single-cell raster quadratic in a
    number that describes nothing about it: at 1e5 m the search took sixteen
    seconds to look at one cell.
    """
    import time

    grid = Grid.filled(1, 1, 5.0)
    for radius in (16, 1024, 1_000_000):
        assert list(grid.window(0, 0, radius, circular=True)) == [(0, 0)]

    def elapsed(radius: int) -> float:
        start = time.perf_counter()
        for _ in range(200):
            list(grid.window(0, 0, radius, circular=True))
        return time.perf_counter() - start

    small = max(elapsed(16), 1e-6)
    assert elapsed(4096) / small < 20.0


@pytest.mark.parametrize("sentinel", [0.0, -9999.0, 1e300, -1e-300])
def test_only_the_sentinel_itself_reads_back_as_absent(sentinel):
    """Absence is an exact value, not a neighbourhood.

    Treating nearby numbers as absent deleted measurements the document had
    written out in full — with a sentinel of ``0`` the smallest positive float
    a raster can hold vanished — and the band of swallowed values widened with
    the magnitude of the sentinel.
    """
    neighbour = math.nextafter(sentinel, math.inf)
    grid = Grid.from_rows([[neighbour, 1.0]], nodata_value=sentinel)
    assert Grid.parse(grid.to_text(precision=None)).values == [neighbour, 1.0]

    absent = Grid.from_rows([[None, 1.0]], nodata_value=sentinel)
    assert Grid.parse(absent.to_text(precision=None)).values == [None, 1.0]


def test_statistics_stay_finite_for_accepted_finite_values():
    """Every value a raster accepts must survive being described.

    A raster holding 0 and 1e200 is finite and is read, written and analysed
    without complaint, but squaring the deviations overflowed — after the
    command had already written its output — and the run reported failure over
    a number nobody had asked for.
    """
    stats = Grid.from_rows([[0.0, 1e200]]).stats()
    assert stats.mean == pytest.approx(5e199)
    assert math.isfinite(stats.stdev) and stats.stdev > 0.0

    extreme = Grid.from_rows([[-1e308, 1e308]]).stats()
    assert extreme.mean == pytest.approx(0.0, abs=1.0)
    assert math.isfinite(extreme.stdev)

    for value in (1e308, 1e-308):
        pair = Grid.from_rows([[value, value]]).stats()
        assert math.isfinite(pair.mean) and pair.stdev == 0.0


@pytest.mark.parametrize("header", ["ncols 1.9\nnrows 1\n", "ncols 1\nnrows 2.5\n"])
def test_fractional_dimensions_are_malformed(header):
    """A raster has a whole number of cells, so 1.9 columns is bad input."""
    with pytest.raises(GridError, match="whole number"):
        Grid.parse(header + "cellsize 1\n0\n")
    assert Grid.parse("ncols 2\nnrows 1\ncellsize 1\n0 1\n").ncols == 2


def _oracle(values):
    """Mean and sample standard deviation computed far outside float range."""
    with decimal.localcontext() as ctx:
        ctx.prec = 400
        exact = [decimal.Decimal(value) for value in values]
        mean = sum(exact) / len(exact)
        if len(exact) < 2:
            return mean, decimal.Decimal(0)
        variance = sum((value - mean) ** 2 for value in exact) / (len(exact) - 1)
        return mean, variance.sqrt()


@pytest.mark.parametrize(
    "values",
    [
        [MAX, MAX, MAX],
        [-MAX, MAX, MAX, MAX],
        [MAX / 4, -MAX / 4],
        [MAX / 2, MAX / 2, -MAX / 2],
        [1e308, 1e-308, 5.0],
        [0.0, 1e200],
        [5e-324, 1e-320, 3e-320],
        [-3.0, -3.0, -3.0],
        [7.5],
    ],
)
def test_statistics_match_an_independent_high_precision_oracle(values):
    """A statistic a raster's own values can express must come out exactly.

    Every list here is made of values the raster accepts and has a mean and a
    deviation inside the finite range, so there is a right answer to give.
    Adding before dividing reached the limit part-way through three copies of
    the largest float and raised, and subtracting a mean of the same size from
    ``-MAX`` produced ``NaN`` — both for statistics a Decimal at 120 digits
    prints without difficulty.

    The oracle runs at 400 significant digits: the largest float needs 309 of
    them, so a shorter context rounds the reference itself and invents a
    deviation where the true one is zero.
    """
    stats = Grid.from_rows([values], cellsize=1.0).stats()
    mean, stdev = _oracle(values)
    assert math.isfinite(stats.mean) and math.isfinite(stats.stdev)
    assert stats.mean == pytest.approx(float(mean), rel=1e-12, abs=1e-320)
    assert stats.stdev == pytest.approx(float(stdev), rel=1e-12, abs=1e-320)


def test_reductions_scale_before_they_sum():
    """The shared helpers are the reason the statistics survive, so test them.

    ``sum`` reaches the finite limit on a running total whose answer is
    representable, and squaring reaches it for any value past about 1.3e154.
    """
    assert mean_of([MAX, MAX, MAX]) == MAX
    assert stdev_of([MAX, MAX, MAX]) == 0.0
    assert sum_of([MAX, MAX, -MAX]) == MAX
    assert sum_of([]) == 0.0
    assert rms_of([3.0, 4.0]) == pytest.approx(math.sqrt(12.5))
    assert rms_of([1e200, 1e200]) == pytest.approx(1e200)
    assert math.isfinite(rms_of([MAX, MAX]))
    with pytest.raises(GridError):
        mean_of([])


def test_smoothing_a_saturated_raster_stays_finite():
    """Smoothing runs before anything else in detection, so it cannot overflow.

    The mean filter summed its window and divided afterwards.  Two cells at the
    largest finite float made that sum infinite, and detection then tried to
    turn an infinite search radius into a whole number of cells.
    """
    smoothed = Grid.from_rows([[MAX, MAX], [MAX, MAX]], cellsize=1.0).smooth_mean(1)
    assert all(value == MAX for value in smoothed.values)


def test_alignment_does_not_depend_on_which_grid_is_asked():
    """``a`` lines up with ``b`` exactly when ``b`` lines up with ``a``.

    The origin tolerance is a fraction of a cell, and cellsizes only have to
    agree to a relative 1e-9.  Reading that fraction off the receiver alone
    left a band in which ``a.combine(b)`` refused the pair that
    ``b.combine(a)`` accepted, so whether two rasters could be subtracted
    depended on the order they were written in.
    """
    a = Grid.filled(1, 1, 1.0, cellsize=1.0)
    b = Grid.filled(1, 1, 2.0, cellsize=1.0 + 2**-33)
    b.xllcorner = 1.00000000005e-06
    with pytest.raises(GridError):
        a.assert_aligned(b)
    with pytest.raises(GridError):
        b.assert_aligned(a)

    same = Grid.filled(1, 1, 3.0, cellsize=1.0)
    same.xllcorner = 1e-9
    a.assert_aligned(same)
    same.assert_aligned(a)


@pytest.mark.parametrize("origin", [0.0, 1e6, 1e9])
def test_symmetric_tolerance_still_does_not_grow_with_coordinates(origin):
    """The pair-wise tolerance must not reintroduce a magnitude-relative one."""
    a = Grid.filled(1, 1, 1.0, cellsize=1.0)
    b = Grid.filled(1, 1, 1.0, cellsize=1.0)
    a.xllcorner = origin
    b.xllcorner = origin + 0.5
    with pytest.raises(GridError):
        a.assert_aligned(b)
    with pytest.raises(GridError):
        b.assert_aligned(a)


@pytest.mark.parametrize(
    "header",
    [
        "ncols 5\nnrows 1\nncols 1\ncellsize 1\n",
        "ncols 5\nnrows 1\ncellsize 1\ncellsize 2\n",
        "ncols 5\nnrows 1\ncellsize 1\nNODATA_value -9999\nnodata_value -1\n",
    ],
)
def test_a_repeated_header_field_is_malformed(header):
    """Two answers to one question is bad input, not a correction.

    Keeping the last silently discarded the first: ``ncols 5`` followed by
    ``ncols 1`` read a one-column raster out of a five-value body and reported
    no problem at all.
    """
    with pytest.raises(GridError, match="more than once"):
        Grid.parse(header + "1 2 3 4 5\n")
    assert Grid.parse("ncols 5\nnrows 1\ncellsize 1\n1 2 3 4 5\n").ncols == 5


@pytest.mark.parametrize("axis", ["x", "y"])
def test_corner_and_centre_origins_cannot_both_be_given(axis):
    """The two spellings of the origin differ by half a cell.

    Preferring the corner threw away whichever of the two the writer had
    actually measured, and moved the raster half a cell without saying so.
    """
    header = f"ncols 1\nnrows 1\ncellsize 2\n{axis}llcorner 10\n{axis}llcenter 10\n"
    with pytest.raises(GridError, match="both"):
        Grid.parse(header + "1\n")
    centred = Grid.parse(f"ncols 1\nnrows 1\ncellsize 2\n{axis}llcenter 10\n1\n")
    assert getattr(centred, f"{axis}llcorner") == 9.0

"""Raster grids in the ESRI ASCII Grid (``.asc``) format.

Silvispect works on plain-text rasters so that the whole toolkit stays free of
binary geospatial dependencies.  The ESRI ASCII Grid format is a real, widely
supported interchange format for digital surface models (DSM), digital terrain
models (DTM) and canopy height models (CHM), and it is human readable, which
makes fixtures and regression tests easy to reason about.

Layout follows the format specification: row ``0`` is the *northern* row and
column ``0`` is the *western* column.  Cells equal to the ``NODATA_value`` are
stored as ``None`` so that arithmetic never silently propagates sentinel
numbers such as ``-9999`` into a mean height.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Iterable, Iterator, Sequence
from dataclasses import dataclass, field
from pathlib import Path

__all__ = ["Grid", "GridError", "GridStats", "mean_of", "rms_of", "stdev_of", "sum_of"]

NODATA_DEFAULT = -9999.0

Cell = tuple[int, int]
Value = float | None


class GridError(ValueError):
    """Raised when a grid cannot be parsed or two grids are incompatible."""


def sum_of(values: Sequence[float]) -> float:
    """Total that stays exact, and that overflows only if the answer does.

    ``sum`` rounds once per term and can pass the finite limit on a running
    total whose final value is representable.  Accumulating in units of the
    largest magnitude present removes both problems: the scaling is undone
    once, at the end, so infinity is returned only when the true total really
    is out of range.
    """
    if not values:
        return 0.0
    scale = max(abs(value) for value in values)
    if scale == 0.0 or not math.isfinite(scale):
        return math.fsum(values)
    return scale * math.fsum(value / scale for value in values)


def rms_of(values: Sequence[float]) -> float:
    """Root mean square that cannot overflow on finite input.

    Squaring first overflows for any magnitude above about ``1.3e154``, which
    is far inside the range of diameters a stand table accepts, and yields
    infinity for a quadratic mean that is perfectly ordinary.
    """
    count = len(values)
    if count == 0:
        raise GridError("root mean square of no values is undefined")
    scale = max(abs(value) for value in values)
    if scale == 0.0:
        return 0.0
    return scale * math.sqrt(math.fsum((value / scale) ** 2 for value in values) / count)


def _check_finite(value: Value) -> None:
    """Refuse a cell value the writer could not put in a file."""
    if value is not None and not math.isfinite(value):
        raise GridError(f"cell value {value!r} is not a finite number")


def mean_of(values: Sequence[float]) -> float:
    """Arithmetic mean that cannot overflow on finite input.

    Every value a raster accepts must survive being summarised.  Adding first
    and dividing afterwards overflows on three cells holding the largest finite
    float, and dividing each term first rounds the running total just past the
    limit on the same input — in both cases for a mean that is perfectly
    representable.  Working in units of the largest magnitude present keeps
    every term inside ``[-1, 1]``, and the scale is multiplied back in last.
    """
    count = len(values)
    if count == 0:
        raise GridError("mean of no values is undefined")
    scale = max(abs(value) for value in values)
    if scale == 0.0:
        return 0.0
    return scale * (math.fsum(value / scale for value in values) / count)


def stdev_of(values: Sequence[float], mean: float | None = None) -> float:
    """Sample standard deviation that cannot overflow on finite input.

    The deviations are taken in units of the largest magnitude present, so
    ``value - mean`` cannot overflow before the scaling that was meant to
    protect it: a raster of ``-MAX`` and three ``MAX`` has a deviation of
    exactly ``MAX``, which subtracting first turned into ``NaN``.
    """
    count = len(values)
    if count < 2:
        return 0.0
    scale = max(abs(value) for value in values)
    if scale == 0.0:
        return 0.0
    centre = (mean_of(values) if mean is None else mean) / scale
    total = math.fsum((value / scale - centre) ** 2 for value in values)
    return scale * math.sqrt(total / (count - 1))


@dataclass(frozen=True)
class GridStats:
    """Summary statistics over the valid (non-nodata) cells of a grid."""

    count: int
    nodata_count: int
    minimum: float | None
    maximum: float | None
    mean: float | None
    stdev: float | None

    def as_dict(self) -> dict[str, float | int | None]:
        return {
            "count": self.count,
            "nodata_count": self.nodata_count,
            "min": self.minimum,
            "max": self.maximum,
            "mean": self.mean,
            "stdev": self.stdev,
        }


@dataclass
class Grid:
    """A north-up, square-celled raster.

    Attributes:
        ncols: Number of columns.
        nrows: Number of rows.
        xllcorner: Easting of the lower-left corner of the lower-left cell.
        yllcorner: Northing of the lower-left corner of the lower-left cell.
        cellsize: Edge length of a cell in map units (metres by convention).
        nodata_value: Sentinel written for missing cells.
        values: Row-major cell values, ``None`` for nodata.
    """

    ncols: int
    nrows: int
    xllcorner: float
    yllcorner: float
    cellsize: float
    nodata_value: float = NODATA_DEFAULT
    values: list[Value] = field(default_factory=list)

    def __post_init__(self) -> None:
        # The geometry is held to the same standard as the cells: everything
        # here is written into the header of the file, and ``parse`` refuses a
        # header it cannot read as a finite number.  A raster built in memory
        # with a cell size of ``inf`` passed every analysis, was written out,
        # and was then refused by the reader that had just been told about it.
        for name in ("ncols", "nrows"):
            count = getattr(self, name)
            if isinstance(count, bool) or count != int(count):
                raise GridError(f"{name} must be a whole number, got {count!r}")
            setattr(self, name, int(count))
        for name in ("xllcorner", "yllcorner", "cellsize", "nodata_value"):
            value = getattr(self, name)
            if isinstance(value, bool) or not math.isfinite(value):
                raise GridError(f"{name} must be a finite number, got {value!r}")
        if self.ncols <= 0 or self.nrows <= 0:
            raise GridError("grid must have a positive number of rows and columns")
        if self.cellsize <= 0:
            raise GridError("cellsize must be positive")
        if not self.values:
            self.values = [None] * (self.ncols * self.nrows)
        if len(self.values) != self.ncols * self.nrows:
            raise GridError(f"expected {self.ncols * self.nrows} values, got {len(self.values)}")

    # ------------------------------------------------------------------
    # construction
    # ------------------------------------------------------------------
    @classmethod
    def filled(
        cls,
        nrows: int,
        ncols: int,
        value: Value = None,
        *,
        cellsize: float = 1.0,
        xllcorner: float = 0.0,
        yllcorner: float = 0.0,
        nodata_value: float = NODATA_DEFAULT,
    ) -> Grid:
        """Return a grid of ``nrows`` x ``ncols`` cells all set to ``value``.

        Raises:
            GridError: If ``value`` is not finite.  See :meth:`from_rows`.
        """
        _check_finite(value)
        return cls(
            ncols=ncols,
            nrows=nrows,
            xllcorner=xllcorner,
            yllcorner=yllcorner,
            cellsize=cellsize,
            nodata_value=nodata_value,
            values=[value] * (nrows * ncols),
        )

    @classmethod
    def from_rows(
        cls,
        rows: Sequence[Sequence[Value]],
        *,
        cellsize: float = 1.0,
        xllcorner: float = 0.0,
        yllcorner: float = 0.0,
        nodata_value: float = NODATA_DEFAULT,
    ) -> Grid:
        """Build a grid from a sequence of equal-length rows (north row first).

        Raises:
            GridError: If any cell is not finite.  :meth:`parse` has always
                refused ``inf`` and ``NaN`` in a file, so accepting them here
                left a raster that could be built and analysed but not written
                down and read back — which is how a canopy model came to be
                written with ``inf`` in it and then rejected by its own reader.
                A cell with no measurement is ``None``, not a non-number.
        """
        if not rows:
            raise GridError("cannot build a grid from zero rows")
        ncols = len(rows[0])
        if any(len(row) != ncols for row in rows):
            raise GridError("all rows must have the same length")
        values: list[Value] = []
        for row in rows:
            values.extend(row)
        for value in values:
            _check_finite(value)
        return cls(
            ncols=ncols,
            nrows=len(rows),
            xllcorner=xllcorner,
            yllcorner=yllcorner,
            cellsize=cellsize,
            nodata_value=nodata_value,
            values=values,
        )

    def like(self, value: Value = None) -> Grid:
        """Return an empty grid with identical georeferencing."""
        return Grid(
            ncols=self.ncols,
            nrows=self.nrows,
            xllcorner=self.xllcorner,
            yllcorner=self.yllcorner,
            cellsize=self.cellsize,
            nodata_value=self.nodata_value,
            values=[value] * (self.ncols * self.nrows),
        )

    def copy(self) -> Grid:
        """Return a deep copy of this grid."""
        clone = self.like()
        clone.values = list(self.values)
        return clone

    # ------------------------------------------------------------------
    # geometry
    # ------------------------------------------------------------------
    @property
    def has_observations(self) -> bool:
        """Whether any cell of this raster was actually measured.

        A raster of nothing but nodata covers ground that was never observed.
        It is not a measurement of an empty place, and nothing downstream may
        read it as one.
        """
        return any(value is not None for value in self.values)

    @property
    def cell_area(self) -> float:
        """Area of a single cell in square map units."""
        return self.cellsize * self.cellsize

    @property
    def extent(self) -> tuple[float, float, float, float]:
        """Return ``(xmin, ymin, xmax, ymax)`` of the grid envelope."""
        return (
            self.xllcorner,
            self.yllcorner,
            self.xllcorner + self.ncols * self.cellsize,
            self.yllcorner + self.nrows * self.cellsize,
        )

    @property
    def area(self) -> float:
        """Total footprint of the grid in square map units."""
        return self.ncols * self.nrows * self.cell_area

    @property
    def area_ha(self) -> float:
        """Total footprint of the grid in hectares (assumes metric units)."""
        return self.area / 10_000.0

    def in_bounds(self, row: int, col: int) -> bool:
        return 0 <= row < self.nrows and 0 <= col < self.ncols

    def index(self, row: int, col: int) -> int:
        if not self.in_bounds(row, col):
            raise IndexError(f"cell ({row}, {col}) is outside the grid")
        return row * self.ncols + col

    def cell_center(self, row: int, col: int) -> tuple[float, float]:
        """Return the map coordinates of the centre of a cell."""
        if not self.in_bounds(row, col):
            raise IndexError(f"cell ({row}, {col}) is outside the grid")
        x = self.xllcorner + (col + 0.5) * self.cellsize
        y = self.yllcorner + (self.nrows - row - 0.5) * self.cellsize
        return (x, y)

    def cell_of(self, x: float, y: float) -> Cell | None:
        """Return the cell containing map coordinate ``(x, y)`` or ``None``."""
        col = math.floor((x - self.xllcorner) / self.cellsize)
        row = math.floor((self.yllcorner + self.nrows * self.cellsize - y) / self.cellsize)
        if not self.in_bounds(row, col):
            return None
        return (row, col)

    # ------------------------------------------------------------------
    # element access
    # ------------------------------------------------------------------
    def get(self, row: int, col: int) -> Value:
        return self.values[self.index(row, col)]

    def set(self, row: int, col: int, value: Value) -> None:
        self.values[self.index(row, col)] = value

    def __getitem__(self, cell: Cell) -> Value:
        return self.get(cell[0], cell[1])

    def __setitem__(self, cell: Cell, value: Value) -> None:
        self.set(cell[0], cell[1], value)

    def cells(self) -> Iterator[Cell]:
        """Iterate over every ``(row, col)`` pair in row-major order."""
        for row in range(self.nrows):
            for col in range(self.ncols):
                yield (row, col)

    def valid_cells(self) -> Iterator[tuple[int, int, float]]:
        """Iterate over ``(row, col, value)`` for cells that hold data."""
        for row in range(self.nrows):
            base = row * self.ncols
            for col in range(self.ncols):
                value = self.values[base + col]
                if value is not None:
                    yield (row, col, value)

    def finite_values(self) -> Iterator[float]:
        for value in self.values:
            if value is not None:
                yield value

    def neighbors(self, row: int, col: int, connectivity: int = 8) -> Iterator[Cell]:
        """Yield in-bounds neighbours using 4- or 8-connectivity."""
        offsets: tuple[tuple[int, int], ...]
        if connectivity == 4:
            offsets = ((-1, 0), (1, 0), (0, -1), (0, 1))
        elif connectivity == 8:
            # fmt: off
            offsets = (
                (-1, -1), (-1, 0), (-1, 1),
                (0, -1), (0, 1),
                (1, -1), (1, 0), (1, 1),
            )
            # fmt: on
        else:
            raise GridError("connectivity must be 4 or 8")
        for drow, dcol in offsets:
            nrow, ncol = row + drow, col + dcol
            if self.in_bounds(nrow, ncol):
                yield (nrow, ncol)

    def window(self, row: int, col: int, radius: int, *, circular: bool = False) -> Iterator[Cell]:
        """Yield in-bounds cells within ``radius`` cells of ``(row, col)``.

        The offsets are clipped to the raster before they are walked, so the
        cost is the number of cells actually yielded rather than the area of
        the requested disc.  It has to be: the detector sizes its search window
        from the height of the cell it is testing, and an accepted — if
        implausible — height of 1e5 m asks for a radius of 5,501 cells.  On a
        one-cell raster that is a single comparison worth of real work, but
        walking the unclipped square spent thirty million.
        """
        if radius < 0:
            raise GridError("radius must be non-negative")
        squared = radius * radius
        for drow in range(max(-radius, -row), min(radius, self.nrows - 1 - row) + 1):
            span = radius
            if circular:
                # The widest in-disc offset on this row, so the column range is
                # already circular and needs no per-cell rejection test.
                span = math.isqrt(squared - drow * drow)
            base = row + drow
            for dcol in range(max(-span, -col), min(span, self.ncols - 1 - col) + 1):
                yield (base, col + dcol)

    # ------------------------------------------------------------------
    # whole-grid operations
    # ------------------------------------------------------------------
    def map_values(self, fn: Callable[[float], Value]) -> Grid:
        """Apply ``fn`` to every valid cell, leaving nodata untouched."""
        out = self.like()
        out.values = [None if v is None else fn(v) for v in self.values]
        return out

    def combine(self, other: Grid, fn: Callable[[float, float], Value]) -> Grid:
        """Combine two aligned grids cell by cell.

        Cells where either input is nodata become nodata in the result.
        """
        self.assert_aligned(other)
        out = self.like()
        out.values = [
            None if (a is None or b is None) else fn(a, b)
            for a, b in zip(self.values, other.values, strict=True)
        ]
        return out

    def assert_aligned(self, other: Grid) -> None:
        """Raise :class:`GridError` unless two grids share the same geometry.

        Origins must agree to within a millionth of a cell.  The tolerance is a
        fraction of the cell rather than a fixed distance because that is the
        only scale the answer means anything on: a centimetre is nothing on a
        10 m cell and a thousand cells on a millimetre one.  Expressing it in
        metres, or leaving a relative tolerance on the coordinates themselves,
        lets the allowance grow with the coordinate values — at an easting of
        1e9 a half-cell shift passed as aligned, and two rasters a half cell
        apart were subtracted cell for cell into a quietly wrong canopy model.

        Alignment is a property of the pair, so the cell the tolerance is a
        fraction of is taken from both grids rather than from whichever one the
        call was written on.  Cellsizes only have to agree to a relative 1e-9,
        and reading the fraction off the receiver alone left a band in which
        ``a.combine(b)`` refused the very same pair that ``b.combine(a)``
        accepted.  The smaller cell is the honest choice: it is the finer of
        the two resolutions being claimed to line up.
        """
        if (self.ncols, self.nrows) != (other.ncols, other.nrows):
            raise GridError("grids differ in shape")
        if not math.isclose(self.cellsize, other.cellsize, rel_tol=1e-9):
            raise GridError("grids differ in cellsize")
        tolerance = min(abs(self.cellsize), abs(other.cellsize)) * 1e-6
        if (
            abs(self.xllcorner - other.xllcorner) > tolerance
            or abs(self.yllcorner - other.yllcorner) > tolerance
        ):
            raise GridError("grids differ in origin")

    def clip(self, minimum: float | None = None, maximum: float | None = None) -> Grid:
        """Clamp valid values into ``[minimum, maximum]``.

        Raises:
            GridError: If a bound is not a finite number.  Clamping is a floor
                and a ceiling, and neither ``inf`` nor ``NaN`` is one: a floor
                of ``inf`` replaced every measured height with infinity, which
                the writer then put in the file and the reader refused.  A
                raster is only worth returning if it can be written down.
        """
        for name, bound in (("minimum", minimum), ("maximum", maximum)):
            if bound is not None and not math.isfinite(bound):
                raise GridError(f"{name} must be a finite number, got {bound!r}")

        def _clip(value: float) -> float:
            if minimum is not None:
                value = max(minimum, value)
            if maximum is not None:
                value = min(maximum, value)
            return value

        return self.map_values(_clip)

    def replace_nodata(self, value: float) -> Grid:
        out = self.like()
        out.values = [value if v is None else v for v in self.values]
        return out

    def focal(
        self,
        fn: Callable[[list[float]], Value],
        radius: int = 1,
        *,
        circular: bool = False,
        min_valid: int = 1,
    ) -> Grid:
        """Apply a moving-window reduction over the grid.

        ``fn`` receives the list of valid values inside the window.  Windows
        holding fewer than ``min_valid`` values yield nodata.
        """
        out = self.like()
        for row, col, _ in self.valid_cells():
            sample = [
                self.values[r * self.ncols + c]
                for r, c in self.window(row, col, radius, circular=circular)
            ]
            valid = [v for v in sample if v is not None]
            if len(valid) < min_valid:
                continue
            out.values[row * self.ncols + col] = fn(valid)
        return out

    def smooth_mean(self, radius: int = 1) -> Grid:
        """Mean filter — suppresses fine noise before tree detection."""
        return self.focal(mean_of, radius, circular=True)

    def smooth_median(self, radius: int = 1) -> Grid:
        """Median filter — removes spikes while preserving crown edges."""
        return self.focal(_median, radius, circular=True)

    def focal_max(self, radius: int = 1) -> Grid:
        return self.focal(max, radius, circular=True)

    def stats(self) -> GridStats:
        values = list(self.finite_values())
        nodata = len(self.values) - len(values)
        if not values:
            return GridStats(0, nodata, None, None, None, None)
        count = len(values)
        mean = mean_of(values)
        if count > 1:
            stdev = stdev_of(values, mean)
        else:
            stdev = 0.0
        return GridStats(count, nodata, min(values), max(values), mean, stdev)

    def histogram(self, bins: int = 10, lo: float | None = None, hi: float | None = None):
        """Return ``(edges, counts)`` for the valid values of the grid."""
        values = list(self.finite_values())
        if bins <= 0:
            raise GridError("bins must be positive")
        if not values:
            return ([], [])
        lo = min(values) if lo is None else lo
        hi = max(values) if hi is None else hi
        if math.isclose(hi, lo):
            hi = lo + 1.0
        width = (hi - lo) / bins
        counts = [0] * bins
        for value in values:
            if value < lo or value > hi:
                continue
            slot = min(bins - 1, int((value - lo) / width))
            counts[slot] += 1
        edges = [lo + i * width for i in range(bins + 1)]
        return (edges, counts)

    # ------------------------------------------------------------------
    # serialisation
    # ------------------------------------------------------------------
    @classmethod
    def parse(cls, text: str) -> Grid:
        """Parse an ESRI ASCII Grid document."""
        header: dict[str, float] = {}
        tokens: list[str] = []
        lines = text.splitlines()
        cursor = 0
        header_keys = {
            "ncols",
            "nrows",
            "xllcorner",
            "yllcorner",
            "xllcenter",
            "yllcenter",
            "cellsize",
            "nodata_value",
        }
        while cursor < len(lines):
            parts = lines[cursor].split()
            if not parts:
                cursor += 1
                continue
            key = parts[0].lower()
            if key not in header_keys:
                break
            if len(parts) < 2:
                raise GridError(f"header line {cursor + 1} has no value")
            # Repeating a field is not a correction: the reader cannot know
            # which of the two the writer meant, and keeping the last silently
            # threw the first away.  "ncols 5" followed by "ncols 50" read a
            # different raster out of the same body with no complaint.
            if key in header:
                raise GridError(f"header field {key!r} is given more than once")
            try:
                header[key] = float(parts[1])
            except ValueError as exc:
                raise GridError(f"bad header value on line {cursor + 1}") from exc
            if not math.isfinite(header[key]):
                raise GridError(f"header field {key!r} must be a finite number, got {parts[1]!r}")
            cursor += 1
        for missing in ("ncols", "nrows", "cellsize"):
            if missing not in header:
                raise GridError(f"missing required header field {missing!r}")
        for line in lines[cursor:]:
            tokens.extend(line.split())

        # A raster has a whole number of cells.  Truncating "ncols 1.9" to 1
        # accepts a header nobody meant to write and reads the wrong number of
        # values out of the body, so it is malformed input rather than a
        # rounding question.
        dimensions: dict[str, int] = {}
        for name in ("ncols", "nrows"):
            size = header[name]
            if size != int(size):
                raise GridError(f"header field {name!r} must be a whole number, got {size!r}")
            dimensions[name] = int(size)
        ncols = dimensions["ncols"]
        nrows = dimensions["nrows"]
        cellsize = header["cellsize"]
        nodata = header.get("nodata_value", NODATA_DEFAULT)
        # Corner and centre are two ways of stating the same origin, and they
        # disagree by half a cell.  Preferring the corner discarded whichever
        # of the two the writer had actually measured, so a header that states
        # both is asked to state one.
        for corner, centre in (("xllcorner", "xllcenter"), ("yllcorner", "yllcenter")):
            if corner in header and centre in header:
                raise GridError(f"header gives both {corner!r} and {centre!r}")
        if "xllcorner" in header:
            xll = header["xllcorner"]
        elif "xllcenter" in header:
            xll = header["xllcenter"] - cellsize / 2.0
        else:
            xll = 0.0
        if "yllcorner" in header:
            yll = header["yllcorner"]
        elif "yllcenter" in header:
            yll = header["yllcenter"] - cellsize / 2.0
        else:
            yll = 0.0

        if len(tokens) != ncols * nrows:
            raise GridError(f"expected {ncols * nrows} cell values, found {len(tokens)}")
        values: list[Value] = []
        for token in tokens:
            try:
                number = float(token)
            except ValueError as exc:
                raise GridError(f"non-numeric cell value {token!r}") from exc
            if not math.isfinite(number):
                # A height of infinity is not a measurement, and it would travel
                # all the way to a report that cannot serialise it.
                raise GridError(f"cell value {token!r} is not a finite number")
            # Absence is written as the sentinel the header declares and is read
            # back by matching it exactly.  Treating merely *nearby* values as
            # absent silently deleted measurements: with a sentinel of 0 the
            # smallest positive number a float can hold, 5e-324, is written out
            # in full and read back as nodata, and the band of swallowed values
            # widened with the magnitude of the sentinel.
            values.append(None if number == nodata else number)
        return cls(ncols, nrows, xll, yll, cellsize, nodata, values)

    @classmethod
    def read(cls, path: str | Path) -> Grid:
        return cls.parse(Path(path).read_text(encoding="utf-8"))

    def to_text(self, *, precision: int | None = 3) -> str:
        """Serialise to an ESRI ASCII Grid document.

        ``precision`` is the number of decimals to write; ``None`` writes the
        shortest decimal string that reads back as the same float, so the
        document round-trips exactly.  Any fixed number of decimals quantises,
        and quantisation is not a matter of degree here: two heights that differ
        below the last decimal written become one plateau, which is a different
        canopy surface.  Choosing a *larger* fixed precision only moves the
        magnitude at which that happens.
        """
        head = [
            f"ncols {self.ncols}",
            f"nrows {self.nrows}",
            f"xllcorner {_fmt(self.xllcorner, precision)}",
            f"yllcorner {_fmt(self.yllcorner, precision)}",
            f"cellsize {_fmt(self.cellsize, precision)}",
            f"NODATA_value {_fmt(self.nodata_value, precision)}",
        ]
        body = []
        for row in range(self.nrows):
            start = row * self.ncols
            chunk = self.values[start : start + self.ncols]
            body.append(
                " ".join(_fmt(self.nodata_value if v is None else v, precision) for v in chunk)
            )
        return "\n".join(head + body) + "\n"

    def write(self, path: str | Path, *, precision: int | None = 3) -> Path:
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(self.to_text(precision=precision), encoding="utf-8")
        return target

    def __repr__(self) -> str:  # pragma: no cover - debugging aid
        return (
            f"Grid(nrows={self.nrows}, ncols={self.ncols}, "
            f"cellsize={self.cellsize}, origin=({self.xllcorner}, {self.yllcorner}))"
        )


def _median(values: Iterable[float]) -> float:
    ordered = sorted(values)
    n = len(ordered)
    mid = n // 2
    if n % 2:
        return ordered[mid]
    return (ordered[mid - 1] + ordered[mid]) / 2.0


def _fmt(value: float, precision: int | None) -> str:
    if precision is None:
        # ``repr`` of a float is the shortest string that reads back identically,
        # so the value survives the round trip whatever its magnitude.
        return repr(float(value))
    text = f"{value:.{precision}f}"
    if "." in text:
        text = text.rstrip("0").rstrip(".")
    return text or "0"

"""Canopy height model derivation and canopy-level descriptors.

A canopy height model (CHM) is the height of vegetation above the ground,
normally obtained by subtracting a digital terrain model (DTM) from a digital
surface model (DSM).  Everything downstream in Silvispect — tree detection,
gap analysis, structural metrics — reads a CHM.
"""

from __future__ import annotations

import math
from collections import deque
from dataclasses import dataclass

from .geometry import Walls, distance_field, inscribed_circle, opening_field
from .grid import Grid, GridError, mean_of, stdev_of

__all__ = [
    "STRATA_BREAKS",
    "CanopyGap",
    "canopy_cover",
    "canopy_height_model",
    "distance_to_canopy",
    "find_gaps",
    "gap_fraction",
    "rugosity",
    "vertical_strata",
]

#: Default vertical strata breaks in metres: ground, regeneration, understorey,
#: mid-storey, canopy, emergent.
STRATA_BREAKS: tuple[float, ...] = (0.5, 2.0, 5.0, 15.0, 25.0)

STRATA_NAMES: tuple[str, ...] = (
    "ground",
    "regeneration",
    "understorey",
    "midstorey",
    "canopy",
    "emergent",
)


def canopy_height_model(dsm: Grid, dtm: Grid, *, floor: float = 0.0) -> Grid:
    """Return ``dsm - dtm`` clamped at ``floor``.

    Negative heights are a normalisation artefact rather than a physical
    quantity, so they are clamped instead of being discarded: a cell that
    measures slightly below ground is still a measured, treeless cell.

    Raises:
        GridError: If a cell's height is not representable.  A surface at the
            largest finite float over a terrain at its negative has a height
            that no float can hold; the subtraction returned infinity, the
            writer put ``inf`` in the file, and the reader — which forbids
            non-finite cells — then refused the model the command had just
            said it had written.  The refusal belongs here, before anything is
            written.
    """

    def height(surface: float, terrain: float) -> float:
        value = surface - terrain
        if not math.isfinite(value):
            raise GridError(
                f"canopy height {surface!r} - {terrain!r} is not a representable number"
            )
        return value

    chm = dsm.combine(dtm, height)
    return chm.clip(minimum=floor)


def canopy_cover(chm: Grid, threshold: float = 2.0) -> float:
    """Fraction of valid cells whose height is at or above ``threshold``.

    Returns ``0.0`` for an entirely empty grid.
    """
    total = 0
    covered = 0
    for _, _, height in chm.valid_cells():
        total += 1
        if height >= threshold:
            covered += 1
    if total == 0:
        return 0.0
    return covered / total


@dataclass(frozen=True)
class CanopyGap:
    """A connected patch of canopy below the cover threshold."""

    gap_id: int
    cells: tuple[tuple[int, int], ...]
    area: float
    centroid_x: float
    centroid_y: float
    mean_height: float
    max_height: float
    touches_edge: bool
    inscribed_radius: float = 0.0

    @property
    def area_ha(self) -> float:
        return self.area / 10_000.0

    @property
    def equivalent_radius(self) -> float:
        """Radius of a circle with the same area, in map units."""
        return math.sqrt(self.area / math.pi)

    @property
    def width(self) -> float:
        """Diameter of the largest circle that fits inside the gap."""
        return 2.0 * self.inscribed_radius

    def as_dict(self) -> dict[str, float | int | bool]:
        return {
            "gap_id": self.gap_id,
            "cell_count": len(self.cells),
            "area_m2": round(self.area, 2),
            "equivalent_radius_m": round(self.equivalent_radius, 2),
            "width_m": round(self.width, 2),
            "centroid_x": round(self.centroid_x, 2),
            "centroid_y": round(self.centroid_y, 2),
            "mean_height_m": round(self.mean_height, 2),
            "max_height_m": round(self.max_height, 2),
            "touches_edge": self.touches_edge,
        }


def _blocked(chm: Grid, threshold: float) -> list[bool]:
    """Cells a circle cannot enter: canopy at or above ``threshold``, or unknown."""
    return [value is None or value >= threshold for value in chm.values]


def _canopy_distance_field(chm: Grid, threshold: float) -> list[float]:
    """Exact distance from every cell centre to the nearest canopy cell or plot edge.

    **The plot edge counts as a boundary.**  Measuring only to canopy would let
    an opening that runs off the side of the raster claim an inscribed circle
    larger than the raster itself — a 10 m x 10 m plot with one corner tree
    reported a 25 m gap width, which no circle inside it could have.  Nothing is
    known about the forest beyond the extent, so the honest bound is the extent:
    the reported width describes the opening *as mapped*, and ``touches_edge``
    marks it as possibly continuing further.  It also keeps every value finite
    on a raster with no canopy at all.

    **Distances are Euclidean, to the canopy's boundary.**  Two earlier forms
    of this field were approximations.  A chamfer walk measured centre to
    centre, so a one-cell opening claimed a two-cell circle; seeding the walk at
    the boundary fixed its first step and left every later one octile — a
    knight's move counted as ``1 + sqrt(2)`` where the crow flies ``sqrt(5)``,
    which let a plus-shaped opening with a 3.16 m circle survive a 3.3 m
    minimum.  The field is now a separable exact transform against the cell
    squares (see :mod:`silvispect.geometry`), in map units.
    """
    return [
        value * chm.cellsize
        for value in distance_field(_blocked(chm, threshold), chm.nrows, chm.ncols)
    ]


def distance_to_canopy(chm: Grid, *, threshold: float = 2.0) -> Grid:
    """Distance from every cell to the nearest canopy cell or plot edge.

    Cells at or above ``threshold`` are walls and get distance zero.  Cells
    with no data are treated as unknown rather than as openings, so they are
    walls too.  The distance is exact and Euclidean, measured to the nearest
    wall's *boundary*, and the plot edge bounds it — see
    :func:`_canopy_distance_field` for why.
    """
    out = chm.like()
    # Always finite: the plot edge bounds every cell, so there is no unreachable
    # cell left to report as absent.
    out.values = list(_canopy_distance_field(chm, threshold))
    return out


def _opening_values(
    chm: Grid,
    threshold: float,
    *,
    at_least: float = 0.0,
    walls: Walls | None = None,
    cells: list[float] | None = None,
) -> list[float]:
    """The opening function in cell units — see :func:`silvispect.geometry.opening_field`."""
    if walls is None:
        walls = Walls(_blocked(chm, threshold), chm.nrows, chm.ncols)
    if cells is None:
        cells = distance_field(walls.blocked, chm.nrows, chm.ncols)
    return opening_field(walls, cells, at_least=at_least)


def opening_radius_field(chm: Grid, *, threshold: float = 2.0) -> list[float]:
    """For every sub-canopy cell, the radius of the largest opening covering it.

    This is the *opening function* of the sub-canopy area: cell ``c`` gets the
    largest ``r`` such that some disc of radius ``r`` fits entirely inside the
    opening and still covers ``c``.  Thresholding it at ``r`` yields exactly the
    morphological opening by a disc of that radius.

    Computing the field once and thresholding it is not merely tidier than an
    erosion followed by a dilation — it is the only way to get a *monotone*
    answer.  Composing the two steps mixes two quantisations (a cut on the
    distance field, then a fresh walk out of the eroded core), and the result
    oscillates: on a treeless 10 m x 10 m plot the reported opening ran 64, 96,
    64, 88, 60 cells as the requested width rose 1 m at a time, so *tightening*
    ``min_width`` could enlarge a gap and raise a finding that a looser setting
    did not.  A single scalar field per cell cannot do that: a larger threshold
    always selects a subset.

    The field is exact: the union of every circle that fits, with the circles
    enumerated in closed form from the walls rather than read off cell centres
    — see :func:`silvispect.geometry.opening_field`.  Values are in map units.
    """
    return [value * chm.cellsize for value in _opening_values(chm, threshold)]


def _opened_gap_mask(
    chm: Grid,
    threshold: float,
    radius: float,
    *,
    walls: Walls | None = None,
    cells: list[float] | None = None,
) -> tuple[list[bool], list[float]]:
    """Sub-canopy cells covered by an opening of at least ``radius``.

    Returns:
        The opened mask and the distance-to-boundary values behind it.
    """
    if walls is None:
        walls = Walls(_blocked(chm, threshold), chm.nrows, chm.ncols)
    if cells is None:
        cells = distance_field(walls.blocked, chm.nrows, chm.ncols)
    distances = [value * chm.cellsize for value in cells]
    # Circles that cannot reach the requested radius never cross the
    # threshold, so the search for them is skipped; the mask is identical to
    # thresholding the full field.
    radii = _opening_values(
        chm, threshold, at_least=radius / chm.cellsize, walls=walls, cells=cells
    )
    # Equality against the requested radius has to absorb the last digit of a
    # square root.  The tolerance is a fraction of a cell, not a fixed number
    # of metres: an absolute ``1e-9`` is invisible on a 1 m raster but swamps
    # a raster whose cells are nanometres wide, where it admitted every cell
    # regardless of radius.  The comparison is made in cell units, where the
    # field is computed, so the same cells pass at every cell size.
    wanted = radius / chm.cellsize
    mask = [
        value is not None and value < threshold and radii[position] >= wanted - 1e-9
        for position, value in enumerate(chm.values)
    ]
    return mask, distances


def _components_reaching_border(chm: Grid, threshold: float, connectivity: int = 8) -> list[bool]:
    """Flag every sub-canopy cell whose untrimmed component touches the border.

    Computed on the raw threshold mask, before any morphological trimming, so
    the answer describes the opening on the ground rather than what survived
    the erosion.  ``connectivity`` must match the connectivity the gaps
    themselves are labelled with: under four-connectivity two openings that meet
    only at a corner are separate gaps, and inheriting edge contact across that
    corner would mark an interior gap as running off the plot.
    """
    reaches = [False] * len(chm.values)
    seen = [False] * len(chm.values)
    for row, col, value in chm.valid_cells():
        start = row * chm.ncols + col
        if seen[start] or value >= threshold:
            continue
        queue: deque[tuple[int, int]] = deque([(row, col)])
        seen[start] = True
        members = [start]
        touches = False
        while queue:
            r, c = queue.popleft()
            if r in (0, chm.nrows - 1) or c in (0, chm.ncols - 1):
                touches = True
            for nrow, ncol in chm.neighbors(r, c, connectivity=connectivity):
                index = nrow * chm.ncols + ncol
                if seen[index]:
                    continue
                neighbour = chm.values[index]
                if neighbour is None or neighbour >= threshold:
                    continue
                seen[index] = True
                members.append(index)
                queue.append((nrow, ncol))
        if touches:
            for index in members:
                reaches[index] = True
    return reaches


def find_gaps(
    chm: Grid,
    *,
    threshold: float = 2.0,
    min_area: float = 25.0,
    min_width: float = 0.0,
    connectivity: int = 8,
    include_edge_gaps: bool = True,
) -> list[CanopyGap]:
    """Label canopy openings as connected components below ``threshold``.

    A pure threshold rule is not enough on a fine raster: the ground visible
    between neighbouring crowns forms a connected web that would be reported as
    one enormous opening.  ``min_width`` applies the Brokaw criterion by
    morphologically opening the sub-canopy mask with a disc of that diameter,
    which trims the inter-crown web — including the parts of it attached to a
    genuine opening — before the components are labelled.

    Args:
        chm: Canopy height model.
        threshold: Height at or above which a cell counts as canopy.
        min_area: Discard gaps smaller than this area in square map units.
        min_width: Narrowest opening to report, in map units.  ``0`` disables
            the morphological filter and reports the raw threshold components.
        connectivity: ``4`` or ``8`` cell connectivity.
        include_edge_gaps: Keep gaps that touch the grid border.  Such gaps are
            truncated by the extent, so their area is a lower bound.

    Returns:
        Gaps sorted by descending area.
    """
    if min_area < 0:
        raise GridError("min_area must be non-negative")
    if min_width < 0:
        raise GridError("min_width must be non-negative")

    walls = Walls(_blocked(chm, threshold), chm.nrows, chm.ncols)
    cells = distance_field(walls.blocked, chm.nrows, chm.ncols)
    if min_width > 0:
        mask, _ = _opened_gap_mask(chm, threshold, min_width / 2.0, walls=walls, cells=cells)
    else:
        mask = [value is not None and value < threshold for value in chm.values]

    # Morphological opening insets the mask from the border, so a gap that
    # genuinely runs off the plot keeps no cell *on* the border.  Proximity to
    # the border is not a usable stand-in — an interior opening separated from
    # the edge by a single row of trees is just as close.  Ask the untrimmed
    # sub-canopy area instead: this gap continues past the extent exactly when
    # the raw opening it was carved from reaches the border.
    reaches_border = _components_reaching_border(chm, threshold, connectivity)

    seen = [False] * len(chm.values)
    gaps: list[CanopyGap] = []
    next_id = 1
    for row, col, _ in chm.valid_cells():
        start = row * chm.ncols + col
        if seen[start] or not mask[start]:
            continue
        queue: deque[tuple[int, int]] = deque([(row, col)])
        seen[start] = True
        members: list[tuple[int, int]] = []
        heights: list[float] = []
        while queue:
            r, c = queue.popleft()
            members.append((r, c))
            value = chm.values[r * chm.ncols + c]
            if value is not None:  # the mask only covers valid cells
                heights.append(value)
            for nr, nc in chm.neighbors(r, c, connectivity):
                idx = nr * chm.ncols + nc
                if seen[idx] or not mask[idx]:
                    continue
                seen[idx] = True
                queue.append((nr, nc))
        area = len(members) * chm.cell_area
        if area < min_area:
            continue
        touches_edge = any(reaches_border[r * chm.ncols + c] for r, c in members)
        if touches_edge and not include_edge_gaps:
            continue
        # The largest circle that fits inside the open ground with its centre
        # in this gap's cells, found exactly rather than read off cell centres.
        circle = inscribed_circle(walls, cells, [r * chm.ncols + c for r, c in members])
        inscribed = circle.radius * chm.cellsize
        xs, ys = zip(*(chm.cell_center(r, c) for r, c in members), strict=True)
        gaps.append(
            CanopyGap(
                gap_id=next_id,
                cells=tuple(members),
                area=area,
                centroid_x=mean_of(xs),
                centroid_y=mean_of(ys),
                mean_height=mean_of(heights),
                max_height=max(heights),
                touches_edge=touches_edge,
                inscribed_radius=inscribed,
            )
        )
        next_id += 1
    gaps.sort(key=lambda gap: gap.area, reverse=True)
    return [
        CanopyGap(
            gap_id=position,
            cells=gap.cells,
            area=gap.area,
            centroid_x=gap.centroid_x,
            centroid_y=gap.centroid_y,
            mean_height=gap.mean_height,
            max_height=gap.max_height,
            touches_edge=gap.touches_edge,
            inscribed_radius=gap.inscribed_radius,
        )
        for position, gap in enumerate(gaps, start=1)
    ]


def gap_fraction(chm: Grid, *, threshold: float = 2.0) -> float:
    """Fraction of valid cells below ``threshold``.

    Counted from the open cells directly rather than as ``1 - cover``.  The two
    are the same number in arithmetic but not in floating point: three open
    cells in ten are exactly ``0.3`` counted, and ``0.30000000000000004`` taken
    as one minus seven tenths.  Rules phrased as *above* a limit are read
    strictly, so that last bit decided whether a plot sitting exactly on a
    30 % limit was reported as breaching it.
    """
    total = 0
    open_cells = 0
    for _, _, height in chm.valid_cells():
        total += 1
        if height < threshold:
            open_cells += 1
    if total == 0:
        return 0.0
    return open_cells / total


def rugosity(chm: Grid, *, threshold: float = 0.0) -> float:
    """Standard deviation of canopy heights above ``threshold``.

    Canopy rugosity is a compact proxy for vertical structural complexity:
    even-aged plantations score low, multi-layered stands score high.
    """
    heights = [h for _, _, h in chm.valid_cells() if h > threshold]
    if len(heights) < 2:
        return 0.0
    return stdev_of(heights)


def vertical_strata(chm: Grid, breaks: tuple[float, ...] = STRATA_BREAKS) -> dict[str, float]:
    """Return the proportion of canopy cells falling in each vertical stratum."""
    if list(breaks) != sorted(breaks):
        raise GridError("strata breaks must be ascending")
    names = _strata_names(len(breaks) + 1)
    counts = dict.fromkeys(names, 0)
    total = 0
    for _, _, height in chm.valid_cells():
        total += 1
        slot = len(breaks)
        for i, edge in enumerate(breaks):
            if height < edge:
                slot = i
                break
        counts[names[slot]] += 1
    if total == 0:
        return dict.fromkeys(names, 0.0)
    return {name: counts[name] / total for name in names}


def _strata_names(count: int) -> tuple[str, ...]:
    if count == len(STRATA_NAMES):
        return STRATA_NAMES
    return tuple(f"stratum_{i}" for i in range(count))

"""Individual tree detection and crown delineation from a canopy height model.

Two classical steps are implemented with no third-party dependencies:

1. **Variable-window local maxima.**  A fixed search window either merges small
   neighbouring trees or splits large crowns, because crown width grows with
   tree height.  The window radius is therefore a linear function of the height
   of the cell being tested.
2. **Marker-controlled region growing.**  Each maximum seeds a crown.  Cells are
   absorbed tallest-first and only while the surface keeps descending away from
   the apex, which places the boundary in the valley between adjacent crowns.
"""

from __future__ import annotations

import heapq
import math
from collections import deque
from dataclasses import dataclass, field, replace

from .grid import Grid, GridError, mean_of

__all__ = [
    "Crown",
    "DetectionConfig",
    "DetectionResult",
    "TreeTop",
    "crown_radius_limit",
    "detect_trees",
    "find_treetops",
    "segment_crowns",
    "window_radius",
]


@dataclass(frozen=True)
class DetectionConfig:
    """Tunable parameters for the detection pipeline.

    Attributes:
        min_height: Cells below this height cannot belong to a tree.
        smooth_radius: Radius in cells of the mean filter applied before
            searching for maxima.  ``0`` disables smoothing.
        window_intercept: Constant term of the search-window radius, in metres.
        window_slope: Growth of the search-window radius per metre of height.
        crown_intercept: Constant term of the maximum crown radius, in metres.
        crown_slope: Growth of the maximum crown radius per metre of height.
        drop_fraction: A crown cell must stay above this fraction of the apex
            height, which stops crowns bleeding into the understorey.
        min_crown_cells: Crowns smaller than this are treated as noise.
    """

    min_height: float = 2.0
    smooth_radius: int = 1
    window_intercept: float = 0.8
    window_slope: float = 0.055
    crown_intercept: float = 1.0
    crown_slope: float = 0.10
    drop_fraction: float = 0.45
    min_crown_cells: int = 3

    #: Fields that count cells and must therefore be whole numbers.
    _COUNTS = ("smooth_radius", "min_crown_cells")

    def __post_init__(self) -> None:
        for name, value in self.as_dict().items():
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                raise GridError(f"{name} must be a number, got {type(value).__name__}")
            if not math.isfinite(value):
                raise GridError(f"{name} must be a finite number, got {value!r}")
        for name in self._COUNTS:
            value = getattr(self, name)
            if value != int(value):
                raise GridError(f"{name} must be a whole number of cells, got {value!r}")
        if self.min_height < 0:
            raise GridError("min_height must be non-negative")
        if self.smooth_radius < 0:
            raise GridError("smooth_radius must be non-negative")
        if not 0.0 <= self.drop_fraction < 1.0:
            raise GridError("drop_fraction must be in [0, 1)")
        if self.window_intercept <= 0 or self.crown_intercept <= 0:
            raise GridError("window and crown intercepts must be positive")
        if self.min_crown_cells < 1:
            raise GridError("min_crown_cells must be at least 1")

    def as_dict(self) -> dict[str, float | int]:
        return {
            "min_height": self.min_height,
            "smooth_radius": self.smooth_radius,
            "window_intercept": self.window_intercept,
            "window_slope": self.window_slope,
            "crown_intercept": self.crown_intercept,
            "crown_slope": self.crown_slope,
            "drop_fraction": self.drop_fraction,
            "min_crown_cells": self.min_crown_cells,
        }


@dataclass(frozen=True)
class TreeTop:
    """A detected apex: the local maximum that seeds one crown."""

    tree_id: int
    row: int
    col: int
    x: float
    y: float
    height: float


@dataclass(frozen=True)
class Crown:
    """A delineated crown belonging to one :class:`TreeTop`."""

    tree_id: int
    apex: TreeTop
    cells: tuple[tuple[int, int], ...]
    area: float
    mean_height: float
    max_extent: float

    @property
    def height(self) -> float:
        return self.apex.height

    @property
    def x(self) -> float:
        return self.apex.x

    @property
    def y(self) -> float:
        return self.apex.y

    @property
    def radius(self) -> float:
        """Radius of a circle with the same area as the crown."""
        return math.sqrt(self.area / math.pi)

    @property
    def diameter(self) -> float:
        return 2.0 * self.radius

    def as_dict(self) -> dict[str, float | int]:
        return {
            "tree_id": self.tree_id,
            "x": round(self.x, 3),
            "y": round(self.y, 3),
            "height_m": round(self.height, 3),
            "crown_area_m2": round(self.area, 3),
            "crown_diameter_m": round(self.diameter, 3),
            "mean_height_m": round(self.mean_height, 3),
            "max_extent_m": round(self.max_extent, 3),
            "cell_count": len(self.cells),
        }


@dataclass
class DetectionResult:
    """Everything produced by one run of :func:`detect_trees`."""

    crowns: list[Crown]
    config: DetectionConfig
    label_grid: Grid
    smoothed: Grid = field(repr=False)

    def __len__(self) -> int:
        return len(self.crowns)

    def __iter__(self):
        return iter(self.crowns)

    @property
    def tops(self) -> list[TreeTop]:
        return [crown.apex for crown in self.crowns]

    def density_per_ha(self, area_ha: float | None = None) -> float:
        """Detected stems per hectare over ``area_ha`` (default: grid extent)."""
        area = self.label_grid.area_ha if area_ha is None else area_ha
        if area <= 0:
            raise GridError("area must be positive")
        return len(self.crowns) / area

    def as_dict(self) -> dict[str, object]:
        return {
            "tree_count": len(self.crowns),
            "density_per_ha": round(self.density_per_ha(), 2),
            "config": self.config.as_dict(),
            "trees": [crown.as_dict() for crown in self.crowns],
        }


def window_radius(height: float, config: DetectionConfig) -> float:
    """Search-window radius in metres for a cell of the given height."""
    return config.window_intercept + config.window_slope * max(height, 0.0)


def crown_radius_limit(height: float, config: DetectionConfig) -> float:
    """Maximum plausible crown radius in metres for a tree of this height."""
    return config.crown_intercept + config.crown_slope * max(height, 0.0)


def find_treetops(chm: Grid, config: DetectionConfig | None = None) -> list[TreeTop]:
    """Locate local maxima using a height-dependent search window.

    A cell is a candidate apex unless a **strictly higher** cell lies inside
    its own search window.  Equality never suppresses.  Ranking equal cells by
    their position in the array is not a fact about the forest and it reverses
    under reflection: on the four cells ``8 7 4 8`` the left peak was kept and
    grew a three-cell crown down the 7 and the 4, while the mirror image
    ``8 4 7 8`` kept the same *array* position, ran into the rising 7 after one
    cell, and lost the crown to the minimum-size filter — one tree became none
    on a raster that is the same forest seen from the other side.

    A connected run of equal-height candidates is still **one** treetop: a flat
    crown top is one tree, not twenty.  Its apex is the candidate nearest the
    run's centre, and where two are equally near the denser surroundings win —
    both of which mirror with the raster, unlike a row-major first.
    """
    config = config or DetectionConfig()
    candidates: list[tuple[int, int, float]] = []
    for row, col, height in chm.valid_cells():
        if height < config.min_height:
            continue
        radius_cells = max(1, round(window_radius(height, config) / chm.cellsize))
        if any(
            (other := chm.values[nrow * chm.ncols + ncol]) is not None and other > height
            for nrow, ncol in chm.window(row, col, radius_cells, circular=True)
        ):
            continue
        candidates.append((row, col, height))

    height_of = {(row, col): height for row, col, height in candidates}
    seen: set[tuple[int, int]] = set()
    tops: list[TreeTop] = []
    for row, col, height in candidates:
        if (row, col) in seen:
            continue
        # One connected run of exactly equal height is one treetop.
        plateau = [(row, col)]
        seen.add((row, col))
        queue = deque([(row, col)])
        while queue:
            here_row, here_col = queue.popleft()
            for nrow, ncol in chm.neighbors(here_row, here_col, connectivity=8):
                if (nrow, ncol) in seen or height_of.get((nrow, ncol)) != height:
                    continue
                seen.add((nrow, ncol))
                plateau.append((nrow, ncol))
                queue.append((nrow, ncol))
        tops.append(_plateau_apex(chm, plateau, height, len(tops) + 1))
    return _ranked(chm, tops)


def _ranked(chm: Grid, tops: list[TreeTop]) -> list[TreeTop]:
    """Order treetops tallest-first, breaking equal heights by outlook.

    The order is not cosmetic: crowns are seeded in it, so when two equally
    tall apexes reach for the same cell it decides which one keeps it.  Sorting
    by coordinate reverses under reflection and handed the cell to the other
    apex, which on a four-by-four raster lost a whole crown to the minimum-size
    filter.  Heights alone settle almost every pair, so the outlook — which
    costs a pass over the raster — is only computed for the apexes that tie.
    """
    by_height: dict[float, list[TreeTop]] = {}
    for top in tops:
        by_height.setdefault(top.height, []).append(top)

    ordered: list[TreeTop] = []
    for height in sorted(by_height, reverse=True):
        tied = by_height[height]
        if len(tied) > 1:
            tied.sort(key=lambda top: (_outlook(chm, (top.row, top.col)), top.row, top.col))
        ordered.extend(tied)
    return [replace(top, tree_id=index) for index, top in enumerate(ordered, start=1)]


def _outlook(chm: Grid, cell: tuple[int, int]) -> list[tuple[int, float]]:
    """How the whole raster looks from one cell: every height, by distance.

    Squared distances in whole cells and the heights at them are preserved by
    every rotation, reflection and transposition of a raster, so two cells that
    correspond under such a transform have the same outlook — and two that do
    not almost always differ somewhere in it.  That makes it a tie-break which
    travels with the plot, where the cell's own row and column reverse with it.
    """
    row, col = cell
    ncols = chm.ncols
    return sorted(
        ((other_row - row) ** 2 + (other_col - col) ** 2, value)
        for other_row in range(chm.nrows)
        for other_col in range(ncols)
        if (value := chm.values[other_row * ncols + other_col]) is not None
    )


def _plateau_apex(
    chm: Grid,
    plateau: list[tuple[int, int]],
    height: float,
    tree_id: int,
) -> TreeTop:
    """Pick the one cell of an equal-height run that represents the treetop.

    The centre of the run first, which mirrors with it.  Where several cells
    are equally central — an even-sided or lopsided run — the one with the
    lowest outlook wins, and only a run whose members are genuinely
    interchangeable falls through to row and column order.

    Centrality is measured in whole integers.  Dividing by the run length puts
    a rounded mean on both sides of the comparison, and the rounding does not
    survive a reflection: two cells equidistant from the centre of a run tied
    one way up and, by a single ulp, not the other, which chose a different
    apex and grew a different crown from it.  Scaling the offset by the square
    of the run length keeps every term a whole number, and a whole number
    mirrors exactly.
    """
    count = len(plateau)
    row_sum = sum(row for row, _ in plateau)
    col_sum = sum(col for _, col in plateau)

    def offset(cell: tuple[int, int]) -> int:
        return (cell[0] * count - row_sum) ** 2 + (cell[1] * count - col_sum) ** 2

    nearest = min(offset(cell) for cell in plateau)
    central = [cell for cell in plateau if offset(cell) == nearest]
    if len(central) == 1:
        row, col = central[0]
    else:
        row, col = min(central, key=lambda cell: (_outlook(chm, cell), cell))
    x, y = chm.cell_center(row, col)
    return TreeTop(tree_id, row, col, x, y, height)


def segment_crowns(
    chm: Grid,
    tops: list[TreeTop],
    config: DetectionConfig | None = None,
) -> tuple[list[Crown], Grid]:
    """Grow one crown per treetop and return the crowns plus a label grid.

    The label grid holds the ``tree_id`` of the owning crown for every assigned
    cell and nodata elsewhere.
    """
    config = config or DetectionConfig()
    labels: list[int] = [0] * len(chm.values)
    apex_of: dict[int, TreeTop] = {top.tree_id: top for top in tops}

    # (-height, squared cell distance from own apex, owning label, tie-breaker,
    # row, col, label, height of the recruiting cell).  Equally tall cells are
    # expanded nearest-apex first rather than in insertion order: when two
    # crowns reach for the same cell, insertion order is the order the array
    # happened to be scanned in, and it reverses under reflection.  The
    # distance is counted in whole cells, so it travels with the plot exactly —
    # differencing the map coordinates instead leaves a rounding that a
    # reflection does not reproduce, and a tie decided by an ulp is a tie
    # decided by nothing.  Where height and distance both tie, the contest is
    # settled by which crown is ranked first, which is a property of the plot;
    # only two fronts of the *same* crown fall through to insertion order,
    # where the cell is claimed either way.
    heap: list[tuple[float, int, int, int, int, int, int, float]] = []
    counter = 0
    for top in tops:
        idx = top.row * chm.ncols + top.col
        if labels[idx]:
            continue
        labels[idx] = top.tree_id
        counter += 1
        heapq.heappush(
            heap,
            (-top.height, 0, top.tree_id, counter, top.row, top.col, top.tree_id, top.height),
        )

    members: dict[int, list[tuple[int, int]]] = {top.tree_id: [] for top in tops}
    for top in tops:
        members[top.tree_id].append((top.row, top.col))

    while heap:
        _, _, _, _, row, col, label, parent_height = heapq.heappop(heap)
        apex = apex_of[label]
        limit = crown_radius_limit(apex.height, config)
        floor = max(config.min_height, config.drop_fraction * apex.height)
        for nrow, ncol in chm.neighbors(row, col, connectivity=8):
            idx = nrow * chm.ncols + ncol
            if labels[idx]:
                continue
            value = chm.values[idx]
            if value is None or value < floor:
                continue
            if value > parent_height + 1e-9:
                # Rising again: we have crossed into a neighbouring crown.
                continue
            reach = (nrow - apex.row) ** 2 + (ncol - apex.col) ** 2
            if math.sqrt(reach) * chm.cellsize > limit:
                continue
            labels[idx] = label
            members[label].append((nrow, ncol))
            counter += 1
            heapq.heappush(heap, (-value, reach, label, counter, nrow, ncol, label, value))

    crowns: list[Crown] = []
    kept_labels: set[int] = set()
    for top in tops:
        cells = members[top.tree_id]
        if len(cells) < config.min_crown_cells:
            continue
        heights = [chm.values[r * chm.ncols + c] for r, c in cells]
        valid = [h for h in heights if h is not None]
        span = max((r - top.row) ** 2 + (c - top.col) ** 2 for r, c in cells)
        max_extent = math.sqrt(span) * chm.cellsize
        crowns.append(
            Crown(
                tree_id=top.tree_id,
                apex=top,
                cells=tuple(cells),
                area=len(cells) * chm.cell_area,
                mean_height=mean_of(valid),
                max_extent=max_extent,
            )
        )
        kept_labels.add(top.tree_id)

    crowns.sort(key=lambda crown: (-crown.height, crown.tree_id))
    renumbered: list[Crown] = []
    remap: dict[int, int] = {}
    for new_id, crown in enumerate(crowns, start=1):
        remap[crown.tree_id] = new_id
        apex = crown.apex
        renumbered.append(
            Crown(
                tree_id=new_id,
                apex=TreeTop(new_id, apex.row, apex.col, apex.x, apex.y, apex.height),
                cells=crown.cells,
                area=crown.area,
                mean_height=crown.mean_height,
                max_extent=crown.max_extent,
            )
        )

    label_grid = chm.like()
    for position, label in enumerate(labels):
        if label in kept_labels:
            label_grid.values[position] = float(remap[label])
    return renumbered, label_grid


def detect_trees(chm: Grid, config: DetectionConfig | None = None) -> DetectionResult:
    """Run the full detection pipeline: smooth, find maxima, delineate crowns."""
    config = config or DetectionConfig()
    smoothed = chm.smooth_mean(config.smooth_radius) if config.smooth_radius else chm.copy()
    tops = find_treetops(smoothed, config)
    crowns, label_grid = segment_crowns(smoothed, tops, config)
    return DetectionResult(crowns=crowns, config=config, label_grid=label_grid, smoothed=smoothed)

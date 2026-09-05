"""Exact Euclidean geometry of the open ground in a raster.

The sub-canopy area of a canopy height model is a union of square cells, and
the width of an opening is the diameter of the largest circle that fits inside
that union with the canopy cells and the plot exterior as its walls.  Two
things about that circle are easy to get wrong on a grid, and both were.

*Its centre need not be a cell centre.*  A block of four open cells is a square
two cells across, and the circle that fills it is centred on the point where
the four cells meet — a lattice corner, which is as far from every cell centre
as a point can be.  Reading the distance field at cell centres alone reported a
one-cell circle for a two-cell square.

*Its radius is a Euclidean distance, not a walk.*  A chamfer transform that
steps one cell straight and ``sqrt(2)`` cells diagonally is an octile metric:
it measures a knight's move as ``1 + sqrt(2)`` where the crow flies
``sqrt(5)``, and the error does not shrink with resolution.  Seeding the walk at
the true boundary distance fixed the first step and left every later one.

So this module computes the distance field exactly, by a separable transform
against the cell squares rather than their centres, and then finds the circles
*exactly*, by enumerating where a circle can be held.  A circle inside a
rectilinear region is held in place by the walls it touches: by two parallel
walls facing each other, in which case it can slide along the corridor
between them, or by three features — wall lines and reflex corners — that
surround its centre.  Every such candidate is a small closed-form solve, and
the largest circle is one of them.

The *opening function* — for every cell centre, the radius of the largest
circle that fits and still covers it — needs one more kind of circle.  A
circle can be held by two walls and the cell centre itself: it slides along
the ridge between the two walls until the point it must cover reaches its
rim.  Such a circle is no local maximum of clearance, so no enumeration of
pinned circles finds it, yet it is the largest circle covering that cell
whenever the ridge runs away from the cell towards a bigger opening.  Those
are solved for each cell that a painted circle does not already settle.

Working in cell units keeps every wall coordinate an exact integer, and every
solve is done in offsets from the cell centre it concerns, so the same
configuration produces bit-identical numbers in every orientation of the
raster.
"""

from __future__ import annotations

import math
from collections import deque
from collections.abc import Iterable, Iterator
from dataclasses import dataclass

__all__ = [
    "Disc",
    "Segment",
    "Walls",
    "boundary_segments",
    "distance_field",
    "inscribed_circle",
    "local_circles",
    "opening_field",
    "reflex_corners",
]

#: Fraction of a cell that two lengths may differ by and still be "equal".
EPSILON = 1e-9
HALF_DIAGONAL = math.sqrt(0.5)

#: A wall feature: ``("x", level)`` or ``("y", level)`` for a wall line and
#: ``("p", x, y)`` for a reflex corner, in absolute cell coordinates.
Key = tuple


@dataclass(frozen=True)
class Segment:
    """A maximal straight piece of the boundary between open and blocked ground.

    Coordinates are in cell units: cell ``(row, col)`` occupies
    ``[col, col + 1] x [row, row + 1]``.  A horizontal segment lies on the line
    ``y = level`` and spans ``[start, stop]`` in ``x``; a vertical one lies on
    ``x = level`` and spans ``[start, stop]`` in ``y``.
    """

    horizontal: bool
    level: int
    start: int
    stop: int

    def distance(self, x: float, y: float) -> float:
        """Euclidean distance from a point to the segment."""
        return _segment_distance(self.horizontal, self.level, self.start, self.stop, x, y)

    def bounding_distance(self, x: float, y: float) -> float:
        """A cheap lower bound on :meth:`distance`, for pruning."""
        if self.horizontal:
            return max(abs(y - self.level), self.start - x, x - self.stop, 0.0)
        return max(abs(x - self.level), self.start - y, y - self.stop, 0.0)

    @property
    def key(self) -> Key:
        return ("y", self.level) if self.horizontal else ("x", self.level)

    @property
    def endpoints(self) -> tuple[tuple[float, float], tuple[float, float]]:
        if self.horizontal:
            return (float(self.start), float(self.level)), (float(self.stop), float(self.level))
        return (float(self.level), float(self.start)), (float(self.level), float(self.stop))


def _segment_distance(
    horizontal: bool, level: float, start: float, stop: float, x: float, y: float
) -> float:
    if horizontal:
        along, across = x, y - level
    else:
        along, across = y, x - level
    if along < start:
        return math.hypot(start - along, across)
    if along > stop:
        return math.hypot(along - stop, across)
    return abs(across)


@dataclass(frozen=True)
class Disc:
    """A circle that fits inside the open ground, possibly free to slide.

    ``x1 == x2`` and ``y1 == y2`` for a circle pinned at one point.  A circle
    held between two parallel walls can be centred anywhere on a segment of the
    corridor's midline, which is recorded as a horizontal or vertical run so
    that every position along it is painted rather than one representative.
    All lengths are in cell units.  ``pins`` names the wall features the circle
    touches, when they were measured.
    """

    x1: float
    y1: float
    x2: float
    y2: float
    radius: float
    pins: tuple[Key, ...] = ()

    def covers(self, x: float, y: float) -> bool:
        """Whether some position of the circle covers the point."""
        dx = max(self.x1 - x, x - self.x2, 0.0)
        dy = max(self.y1 - y, y - self.y2, 0.0)
        return math.hypot(dx, dy) <= self.radius + EPSILON


# ----------------------------------------------------------------------
# the distance field
# ----------------------------------------------------------------------
def distance_field(blocked: list[bool], nrows: int, ncols: int) -> list[float]:
    """Exact distance from every cell centre to the nearest blocked cell or edge.

    Distances are to the *square* of the blocked cell, not to its centre, and
    to the plot boundary.  In cell units, so the result for a cell whose
    neighbour straight across is blocked is ``0.5`` and for one whose only
    blocked neighbour is diagonal it is ``sqrt(2) / 2``.  Blocked cells are 0.

    The transform is separable.  Within each column the distance to the nearest
    blocked cell in the row direction is a one-dimensional scan; across a row
    the squared distance to the nearest blocked *square* is then the lower
    envelope of one parabola per column, evaluated half a cell to either side
    of the query column — because the column component of the distance to a
    square is ``|dx| - 1/2`` for any column but the cell's own, which is
    ``(x - c')^2`` at ``x = c -/+ 1/2``.  Evaluating both sides and taking the
    smaller never undercounts, since each parabola is at least its true cost on
    the far side.
    """
    inf = math.inf
    # Row-direction distance to the nearest blocked square in the same column.
    vertical = [inf] * (nrows * ncols)
    for col in range(ncols):
        nearest = inf
        for row in range(nrows):
            index = row * ncols + col
            if blocked[index]:
                nearest = 0.0
                vertical[index] = 0.0
                continue
            nearest += 1.0
            vertical[index] = nearest - 0.5 if nearest < inf else inf
        nearest = inf
        for row in range(nrows - 1, -1, -1):
            index = row * ncols + col
            if blocked[index]:
                nearest = 0.0
                continue
            nearest += 1.0
            if nearest < inf:
                vertical[index] = min(vertical[index], nearest - 0.5)

    out = [0.0] * (nrows * ncols)
    for row in range(nrows):
        base = row * ncols
        heights = [vertical[base + col] ** 2 for col in range(ncols)]
        envelope = _lower_envelope(heights)
        for col in range(ncols):
            index = base + col
            if blocked[index]:
                continue
            best = heights[col]
            for query in (col - 0.5, col + 0.5):
                best = min(best, _evaluate(envelope, heights, query))
            edge = min(col + 0.5, ncols - col - 0.5, row + 0.5, nrows - row - 0.5)
            out[index] = min(math.sqrt(best) if best < inf else inf, edge)
    return out


def _lower_envelope(heights: list[float]) -> list[tuple[int, float]]:
    """Lower envelope of the parabolas ``heights[i] + (x - i)^2``.

    Felzenszwalb and Huttenlocher's linear-time construction: a list of
    ``(vertex, from_x)`` pairs, each parabola owning the envelope from
    ``from_x`` up to the next entry's ``from_x``.  Infinite heights are skipped.
    """
    envelope: list[tuple[int, float]] = []
    for vertex, height in enumerate(heights):
        if height == math.inf:
            continue
        while envelope:
            previous, start = envelope[-1]
            crossing = _intersection(previous, heights[previous], vertex, height)
            if crossing > start:
                envelope.append((vertex, crossing))
                break
            envelope.pop()
        else:
            envelope.append((vertex, -math.inf))
    return envelope


def _intersection(left: int, left_height: float, right: int, right_height: float) -> float:
    """Where the parabola at ``right`` drops below the one at ``left``."""
    return ((right_height + right * right) - (left_height + left * left)) / (2.0 * (right - left))


def _evaluate(envelope: list[tuple[int, float]], heights: list[float], x: float) -> float:
    if not envelope:
        return math.inf
    low, high = 0, len(envelope) - 1
    while low < high:
        middle = (low + high + 1) // 2
        if envelope[middle][1] <= x:
            low = middle
        else:
            high = middle - 1
    vertex = envelope[low][0]
    return heights[vertex] + (x - vertex) ** 2


# ----------------------------------------------------------------------
# the walls
# ----------------------------------------------------------------------
def boundary_segments(blocked: list[bool], nrows: int, ncols: int) -> list[Segment]:
    """Every wall of the open ground, merged into maximal straight pieces.

    A wall is a side of an open cell whose neighbour across it is blocked or
    lies outside the plot.  Collinear consecutive walls are one segment: a
    straight canopy edge twenty cells long is one feature, not twenty.
    """
    horizontal: dict[int, list[int]] = {}
    vertical: dict[int, list[int]] = {}
    for row in range(nrows):
        for col in range(ncols):
            if blocked[row * ncols + col]:
                continue
            if row == 0 or blocked[(row - 1) * ncols + col]:
                horizontal.setdefault(row, []).append(col)
            if row == nrows - 1 or blocked[(row + 1) * ncols + col]:
                horizontal.setdefault(row + 1, []).append(col)
            if col == 0 or blocked[row * ncols + col - 1]:
                vertical.setdefault(col, []).append(row)
            if col == ncols - 1 or blocked[row * ncols + col + 1]:
                vertical.setdefault(col + 1, []).append(row)

    segments: list[Segment] = []
    for is_horizontal, table in ((True, horizontal), (False, vertical)):
        for level, cells in table.items():
            cells.sort()
            start = previous = cells[0]
            for cell in cells[1:]:
                if cell != previous + 1:
                    segments.append(Segment(is_horizontal, level, start, previous + 1))
                    start = cell
                previous = cell
            segments.append(Segment(is_horizontal, level, start, previous + 1))
    return segments


def reflex_corners(blocked: list[bool], nrows: int, ncols: int) -> list[tuple[float, float]]:
    """Lattice points where the open ground turns inward.

    A circle inside a region can be held by a corner only where the corner
    points *into* the region: three of the four cells around the lattice point
    are open, or two diagonal ones are and the region pinches to a point.  A
    convex corner — one open cell — is always farther from an interior point
    than one of its two walls, so it never pins anything and is left out,
    which keeps the enumeration of pinned circles small.  Cells beyond the plot
    count as blocked.
    """

    def open_at(row: int, col: int) -> bool:
        return 0 <= row < nrows and 0 <= col < ncols and not blocked[row * ncols + col]

    corners: list[tuple[float, float]] = []
    for y in range(nrows + 1):
        for x in range(ncols + 1):
            nw, ne = open_at(y - 1, x - 1), open_at(y - 1, x)
            sw, se = open_at(y, x - 1), open_at(y, x)
            count = nw + ne + sw + se
            if count == 3 or (count == 2 and nw == se):
                corners.append((float(x), float(y)))
    return corners


class Walls:
    """The walls of the open ground, indexed for neighbourhood queries.

    Segments and reflex corners are bucketed by the cells they border, so
    gathering the features within some distance of a point costs the area of
    that neighbourhood rather than the length of the whole boundary.  When the
    neighbourhood is larger than the boundary is long — a wide plot with few
    trees — a linear scan is used instead.
    """

    BLOCK = 4

    def __init__(self, blocked: list[bool], nrows: int, ncols: int) -> None:
        self.blocked = blocked
        self.nrows, self.ncols = nrows, ncols
        self.segments = boundary_segments(blocked, nrows, ncols)
        self.corners = reflex_corners(blocked, nrows, ncols)
        block = self.BLOCK
        self.block_cols = ncols // block + 1
        self.segment_buckets: dict[int, list[int]] = {}
        self.corner_buckets: dict[int, list[int]] = {}
        for index, segment in enumerate(self.segments):
            seen: set[int] = set()
            for along in range(segment.start, segment.stop):
                for side in (segment.level - 1, segment.level):
                    row, col = (side, along) if segment.horizontal else (along, side)
                    bucket = self._bucket(row, col)
                    if bucket not in seen:
                        seen.add(bucket)
                        self.segment_buckets.setdefault(bucket, []).append(index)
        for index, (x, y) in enumerate(self.corners):
            seen = set()
            for row in (int(y) - 1, int(y)):
                for col in (int(x) - 1, int(x)):
                    bucket = self._bucket(row, col)
                    if bucket not in seen:
                        seen.add(bucket)
                        self.corner_buckets.setdefault(bucket, []).append(index)

    def _bucket(self, row: int, col: int) -> int:
        block = self.BLOCK
        return (max(0, min(row, self.nrows - 1)) // block) * self.block_cols + max(
            0, min(col, self.ncols - 1)
        ) // block

    def open_at(self, x: float, y: float) -> bool:
        """Whether the point lies in the open ground.

        A point at positive clearance is strictly inside some cell, or on a
        line between cells that are all open; either way some cell touching
        it is open.
        """
        nrows, ncols = self.nrows, self.ncols
        for row in (math.floor(y - EPSILON), math.floor(y + EPSILON)):
            for col in (math.floor(x - EPSILON), math.floor(x + EPSILON)):
                if 0 <= row < nrows and 0 <= col < ncols and not self.blocked[row * ncols + col]:
                    return True
        return False

    def _candidates(self, x: float, y: float, reach: float) -> tuple[Iterable[int], Iterable[int]]:
        block = self.BLOCK
        row_lo = max(0, math.floor((y - reach - 1.0) / block))
        row_hi = min((self.nrows - 1) // block, math.floor((y + reach + 1.0) / block))
        col_lo = max(0, math.floor((x - reach - 1.0) / block))
        col_hi = min((self.ncols - 1) // block, math.floor((x + reach + 1.0) / block))
        buckets = (row_hi - row_lo + 1) * (col_hi - col_lo + 1)
        if buckets >= len(self.segments) + len(self.corners):
            return range(len(self.segments)), range(len(self.corners))
        segments: set[int] = set()
        corners: set[int] = set()
        for row in range(row_lo, row_hi + 1):
            base = row * self.block_cols
            for col in range(col_lo, col_hi + 1):
                segments.update(self.segment_buckets.get(base + col, ()))
                corners.update(self.corner_buckets.get(base + col, ()))
        return segments, corners

    def neighbourhood(self, anchor_x: float, anchor_y: float, reach: float) -> _Neighbourhood:
        """The features within ``reach`` of the anchor, in offsets from it."""
        return _Neighbourhood(self, anchor_x, anchor_y, reach)

    def fits(self, anchor_x: float, anchor_y: float, dx: float, dy: float, radius: float) -> bool:
        """Whether a circle centred at ``anchor + (dx, dy)`` has ``radius`` of clearance.

        Distances are measured in offsets from the anchor, with the walls
        shifted by exact amounts, so the answer does not depend on where in
        the raster the anchor lies.
        """
        x, y = anchor_x + dx, anchor_y + dy
        if not self.open_at(x, y):
            return False
        segment_ids, _ = self._candidates(x, y, radius)
        floor = radius - EPSILON
        for index in segment_ids:
            segment = self.segments[index]
            if segment.bounding_distance(x, y) >= radius:
                continue
            level = segment.level - (anchor_y if segment.horizontal else anchor_x)
            start = segment.start - (anchor_x if segment.horizontal else anchor_y)
            stop = segment.stop - (anchor_x if segment.horizontal else anchor_y)
            if _segment_distance(segment.horizontal, level, start, stop, dx, dy) < floor:
                return False
        return True


class _Neighbourhood:
    """The walls near one point, as offsets from that point.

    ``segments`` holds ``(horizontal, level, start, stop)`` in offsets from
    the anchor, with the feature key alongside; ``points`` holds reflex corner
    offsets.  All offsets are exact, being integer coordinates less a
    half-integer anchor, which is what makes the solves orientation-proof.
    """

    def __init__(self, walls: Walls, anchor_x: float, anchor_y: float, reach: float) -> None:
        self.anchor_x, self.anchor_y = anchor_x, anchor_y
        segment_ids, corner_ids = walls._candidates(anchor_x, anchor_y, reach)
        self.segments: list[tuple[bool, float, float, float]] = []
        self.keys: list[Key] = []
        for index in segment_ids:
            segment = walls.segments[index]
            if segment.bounding_distance(anchor_x, anchor_y) > reach:
                continue
            if segment.horizontal:
                relative = (
                    True,
                    segment.level - anchor_y,
                    segment.start - anchor_x,
                    segment.stop - anchor_x,
                )
            else:
                relative = (
                    False,
                    segment.level - anchor_x,
                    segment.start - anchor_y,
                    segment.stop - anchor_y,
                )
            if _segment_distance(*relative, 0.0, 0.0) <= reach:
                self.segments.append(relative)
                self.keys.append(segment.key)
        self.lines_x = sorted({s[1] for s in self.segments if not s[0]})
        self.lines_y = sorted({s[1] for s in self.segments if s[0]})
        self.points: list[tuple[float, float]] = []
        for index in corner_ids:
            cx, cy = walls.corners[index]
            dx, dy = cx - anchor_x, cy - anchor_y
            if math.hypot(dx, dy) <= reach:
                self.points.append((dx, dy))

    def clearance(self, x: float, y: float) -> float:
        return min(
            (_segment_distance(*segment, x, y) for segment in self.segments), default=math.inf
        )

    def key_of_line(self, horizontal: bool, level: float) -> Key:
        if horizontal:
            return ("y", round(level + self.anchor_y))
        return ("x", round(level + self.anchor_x))

    def key_of_point(self, point: tuple[float, float]) -> Key:
        return ("p", point[0] + self.anchor_x, point[1] + self.anchor_y)

    def pins(self, x: float, y: float, radius: float) -> tuple[Key, ...]:
        """The features the circle at ``(x, y)`` touches."""
        slack = 1e-7
        found: list[Key] = []
        for segment, key in zip(self.segments, self.keys, strict=True):
            if abs(_segment_distance(*segment, x, y) - radius) <= slack and key not in found:
                found.append(key)
        for point in self.points:
            if abs(math.hypot(point[0] - x, point[1] - y) - radius) <= slack:
                found.append(self.key_of_point(point))
        return tuple(found)


# ----------------------------------------------------------------------
# closed-form solves
# ----------------------------------------------------------------------
def _clear_runs(
    walls: _Neighbourhood,
    fixed: float,
    horizontal: bool,
    radius: float,
    lo: float,
    hi: float,
) -> list[tuple[float, float]]:
    """Sub-intervals of a corridor midline where a circle of ``radius`` fits.

    Every wall other than the two that define the corridor forbids an open
    interval of the midline; what is left is where the circle can be centred.
    Contact is allowed: a point exactly ``radius`` from a wall is in.
    """
    forbidden: list[tuple[float, float]] = []
    for is_horizontal, level, start, stop in walls.segments:
        if horizontal == is_horizontal:
            across = abs(fixed - level)
            if across >= radius - EPSILON:
                continue
            reach = math.sqrt(max(0.0, radius * radius - across * across))
            forbidden.append((start - reach, stop + reach))
        else:
            if start <= fixed <= stop:
                forbidden.append((level - radius, level + radius))
                continue
            across = min(abs(fixed - start), abs(fixed - stop))
            if across >= radius - EPSILON:
                continue
            reach = math.sqrt(max(0.0, radius * radius - across * across))
            forbidden.append((level - reach, level + reach))
    forbidden.sort()
    runs: list[tuple[float, float]] = []
    cursor = lo
    for begin, end in forbidden:
        begin += EPSILON  # the interval is open: contact at its ends is allowed
        end -= EPSILON
        if end <= cursor:
            continue
        if begin > hi:
            break
        if begin > cursor:
            runs.append((cursor, min(begin, hi)))
        cursor = max(cursor, end)
        if cursor >= hi:
            break
    if cursor <= hi:
        runs.append((cursor, hi))
    return [(a, b) for a, b in runs if a <= b + EPSILON]


def _two_lines_and_a_point(
    x_line: float,
    y_line: float,
    point: tuple[float, float],
    limit: float,
    *,
    in_cell: bool,
) -> Iterator[tuple[float, float]]:
    """Centres ``t`` from both perpendicular lines and ``t`` from the point.

    The lines are offsets from the anchor, which the circle must lie on the
    same side of the lines as — the cell centre it is centred near, or the
    point it must cover — so the diagonal to search is fixed by the lines'
    signs.  With ``in_cell`` the centre is confined to the anchor's cell and
    ``t`` is at most half a cell either way from each line's distance.
    """
    px, py = point
    u, v = x_line - px, y_line - py
    sx = -1.0 if x_line > 0.0 else 1.0
    sy = -1.0 if y_line > 0.0 else 1.0
    if in_cell:
        lo = max(abs(x_line), abs(y_line)) - 0.5 - EPSILON
        hi = min(abs(x_line), abs(y_line), limit) + 0.5 + EPSILON
    else:
        lo, hi = 0.0, limit + EPSILON
    if lo > hi:
        return
    half_sum = u * sx + v * sy
    disc = half_sum * half_sum - (u * u + v * v)
    if disc < 0.0:
        return
    root = math.sqrt(disc)
    for t in (-half_sum - root, -half_sum + root):
        if EPSILON < t <= limit + EPSILON and lo <= t <= hi:
            yield x_line + sx * t, y_line + sy * t


def _pinned_points(walls: _Neighbourhood, limit: float) -> Iterator[tuple[float, float]]:
    """Every point in the cell where a circle of radius at most ``limit`` can be held.

    Walls are lines (``x = k`` or ``y = k``) and reflex corners (points), in
    offsets from the cell centre.  Two parallel lines make a corridor and are
    handled separately; every other way of pinning a circle needs three
    features, and each combination is a small closed-form solve.  The centre
    lies in the cell, so every feature it touches is within half a diagonal
    of the same distance from the cell centre, which prunes the triples to
    features of nearly equal distance.  Spurious candidates cost nothing:
    each is measured against the real walls afterwards.
    """
    lines_x, lines_y = walls.lines_x, walls.lines_y
    tolerance = 2.0 * HALF_DIAGONAL + EPSILON
    # Corners by distance from the cell centre, for windowing.
    points = sorted(walls.points, key=lambda p: math.hypot(p[0], p[1]))
    spans = [math.hypot(p[0], p[1]) for p in points]

    def window(lo: float, hi: float) -> range:
        """Corners whose distance from the cell centre lies in ``[lo, hi]``."""
        first = 0
        while first < len(points) and spans[first] < lo:
            first += 1
        last = first
        while last < len(points) and spans[last] <= hi:
            last += 1
        return range(first, last)

    # A line of each orientation and a corner: the centre lies on a diagonal
    # through the two lines' crossing, a distance ``t`` from each, and ``t``
    # from the corner.
    for y_line in lines_y:
        for x_line in lines_x:
            if abs(abs(x_line) - abs(y_line)) > 1.0 + EPSILON:
                continue
            reach = max(abs(x_line), abs(y_line))
            for i in window(reach - 1.0 - tolerance, reach + 1.0 + tolerance):
                yield from _two_lines_and_a_point(x_line, y_line, points[i], limit, in_cell=True)

    # One line and two corners: the centre is ``r`` off the line and ``r`` from
    # both corners, so it lies on the corners' perpendicular bisector.
    for horizontal, levels in ((True, lines_y), (False, lines_x)):
        for level in levels:
            reach = abs(level)
            if reach - 0.5 > limit + EPSILON:
                continue
            candidates = window(reach - 0.5 - tolerance, reach + 0.5 + tolerance)
            for i in candidates:
                for j in candidates:
                    if j > i:
                        yield from _line_and_two_corners(horizontal, level, points[i], points[j])

    # Three corners: the circumcentre.  Two corners squeezing the circle
    # between them.  Either way all the corners touched are within a cell
    # diagonal of the same distance from the cell centre.
    for i, first in enumerate(points):
        for j in range(i + 1, len(points)):
            if spans[j] - spans[i] > tolerance:
                break
            second = points[j]
            yield (first[0] + second[0]) / 2.0, (first[1] + second[1]) / 2.0
            for k in range(j + 1, len(points)):
                if spans[k] - spans[i] > tolerance:
                    break
                centre = _circumcentre(first, second, points[k])
                if centre is not None:
                    yield centre


def _line_and_two_corners(
    horizontal: bool,
    level: float,
    first: tuple[float, float],
    second: tuple[float, float],
) -> Iterator[tuple[float, float]]:
    """Centres equidistant from one wall line and two corner points."""
    # Work in coordinates where the line is ``v = level``: ``p`` runs along
    # it and ``v`` across it.
    if horizontal:
        p1, q1, p2, q2 = first[0], first[1], second[0], second[1]
    else:
        p1, q1, p2, q2 = first[1], first[0], second[1], second[0]
    # Bisector of the two corners: A p + B v = C.
    a_coef, b_coef = 2.0 * (p1 - p2), 2.0 * (q1 - q2)
    c_coef = p1 * p1 + q1 * q1 - p2 * p2 - q2 * q2
    if abs(a_coef) < EPSILON and abs(b_coef) < EPSILON:
        return
    for sign in (-1.0, 1.0):
        # v = level + sign * r
        if abs(a_coef) > EPSILON:
            # p = alpha + beta * r
            alpha = (c_coef - b_coef * level) / a_coef
            beta = -b_coef * sign / a_coef
            g, h = alpha - p1, level - q1
            # (g + beta r)^2 + (h + sign r)^2 = r^2
            qa = beta * beta
            qb = 2.0 * (g * beta + h * sign)
            qc = g * g + h * h
            for r in _positive_roots(qa, qb, qc):
                p, v = alpha + beta * r, level + sign * r
                yield (p, v) if horizontal else (v, p)
        else:
            r = (c_coef / b_coef - level) / sign
            if r <= EPSILON:
                continue
            v = level + sign * r
            gap = r * r - (v - q1) ** 2
            if gap < 0.0:
                continue
            for p in (p1 - math.sqrt(gap), p1 + math.sqrt(gap)):
                yield (p, v) if horizontal else (v, p)


def _positive_roots(a: float, b: float, c: float) -> Iterator[float]:
    if abs(a) < EPSILON:
        if abs(b) > EPSILON:
            root = -c / b
            if root > EPSILON:
                yield root
        return
    disc = b * b - 4.0 * a * c
    if disc < 0.0:
        return
    sqrt_disc = math.sqrt(disc)
    for root in ((-b - sqrt_disc) / (2.0 * a), (-b + sqrt_disc) / (2.0 * a)):
        if root > EPSILON:
            yield root


def _circumcentre(
    a: tuple[float, float], b: tuple[float, float], c: tuple[float, float]
) -> tuple[float, float] | None:
    d = 2.0 * (a[0] * (b[1] - c[1]) + b[0] * (c[1] - a[1]) + c[0] * (a[1] - b[1]))
    if abs(d) < EPSILON:
        return None
    a2, b2, c2 = a[0] ** 2 + a[1] ** 2, b[0] ** 2 + b[1] ** 2, c[0] ** 2 + c[1] ** 2
    ux = (a2 * (b[1] - c[1]) + b2 * (c[1] - a[1]) + c2 * (a[1] - b[1])) / d
    uy = (a2 * (c[0] - b[0]) + b2 * (a[0] - c[0]) + c2 * (b[0] - a[0])) / d
    return ux, uy


def _through_point(a: Key, b: Key, limit: float) -> Iterator[tuple[float, float]]:
    """Centres of circles touching features ``a`` and ``b`` and passing through the origin.

    Features are given in offsets from the point the circle must pass
    through, which plays the part of a third feature: a corner at the origin.
    Two parallel lines are left out — a circle between them can slide, and
    the corridor it sweeps is painted directly.
    """
    origin = (0.0, 0.0)
    if a[0] == "p" and b[0] == "p":
        centre = _circumcentre((a[1], a[2]), (b[1], b[2]), origin)
        if centre is not None:
            yield centre
        return
    if a[0] == "p" or b[0] == "p":
        line, point = (b, a) if a[0] == "p" else (a, b)
        horizontal = line[0] == "y"
        yield from _line_and_two_corners(horizontal, line[1], (point[1], point[2]), origin)
        return
    if a[0] == b[0]:
        return
    x_line, y_line = (a[1], b[1]) if a[0] == "x" else (b[1], a[1])
    yield from _two_lines_and_a_point(x_line, y_line, origin, limit, in_cell=False)


# ----------------------------------------------------------------------
# the circles
# ----------------------------------------------------------------------
def local_circles(
    walls: Walls,
    distances: list[float],
    cell: int,
    *,
    allowed: Iterable[int] | None = None,
    at_least: float = 0.0,
) -> list[Disc]:
    """Every held circle centred inside ``cell`` that reaches its centre value.

    ``distances`` is the exact cell-centre field.  Only walls that can be the
    nearest to some point of the cell are gathered, which is a thin ring around
    it: nothing is closer to its centre than its own field value, and nothing
    farther than that plus two half-diagonals can be the nearest wall to any
    point inside it.  Before enumerating anything, the farthest point of the
    cell from each wall bounds what any circle inside it can reach; if that
    falls short of ``at_least``, the cell has nothing to add.  ``allowed``
    restricts the cells a centre may lie in, for measuring one gap of an
    opening the filter has divided.

    Everything is solved in offsets from the cell centre, so the circles come
    out bit-identical however the raster is turned.
    """
    ncols, nrows = walls.ncols, walls.nrows
    row, col = divmod(cell, ncols)
    centre_x, centre_y = col + 0.5, row + 0.5
    limit = distances[cell] + HALF_DIAGONAL + EPSILON
    if limit < at_least - EPSILON:
        return []
    near = walls.neighbourhood(centre_x, centre_y, distances[cell] + 2.0 * HALF_DIAGONAL)

    # Nothing is held by fewer than two features, and two hold a circle only
    # when they face each other: parallel lines or a pinch between corners.
    # Most cells of open ground have a single wall within reach and end here.
    lines_x, lines_y, points = near.lines_x, near.lines_y, near.points
    count = len(lines_x) + len(lines_y) + len(points)
    if count < 2 or (count == 2 and max(len(lines_x), len(lines_y), len(points)) < 2):
        return []

    # The farthest corner of the cell from each wall caps what fits inside it.
    box = ((-0.5, -0.5), (0.5, -0.5), (-0.5, 0.5), (0.5, 0.5))
    ceiling = min(
        (max(_segment_distance(*segment, x, y) for x, y in box) for segment in near.segments),
        default=math.inf,
    )
    if ceiling < at_least - EPSILON:
        return []
    limit = min(limit, ceiling + EPSILON)
    permitted = None if allowed is None else set(allowed)

    def inside(x: float, y: float) -> bool:
        if not (-0.5 - EPSILON <= x <= 0.5 + EPSILON and -0.5 - EPSILON <= y <= 0.5 + EPSILON):
            return False
        # A point at positive clearance is strictly inside some cell; the
        # cells touching it are all open or all blocked.
        ax, ay = centre_x + x, centre_y + y
        for r in (math.floor(ay - EPSILON), math.floor(ay + EPSILON)):
            for c in (math.floor(ax - EPSILON), math.floor(ax + EPSILON)):
                if not (0 <= r < nrows and 0 <= c < ncols):
                    continue
                index = r * ncols + c
                if not walls.blocked[index] and (permitted is None or index in permitted):
                    return True
        return False

    found: list[Disc] = []
    floor = max(distances[cell], at_least) - EPSILON

    # Corridors: two parallel walls facing each other.
    for lines, horizontal in ((near.lines_y, True), (near.lines_x, False)):
        for i, low in enumerate(lines):
            for high in lines[i + 1 :]:
                radius = (high - low) / 2.0
                if radius > limit or radius < floor:
                    continue
                mid = (low + high) / 2.0
                for start, stop in _clear_runs(near, mid, horizontal, radius, -0.5, 0.5):
                    probe = (start + stop) / 2.0
                    x, y = (probe, mid) if horizontal else (mid, probe)
                    if not inside(x, y) or near.clearance(x, y) < radius - EPSILON:
                        continue
                    pins = near.pins(
                        start if horizontal else mid, mid if horizontal else start, radius
                    )
                    pins += tuple(
                        key
                        for key in near.pins(
                            stop if horizontal else mid, mid if horizontal else stop, radius
                        )
                        if key not in pins
                    )
                    if horizontal:
                        found.append(
                            Disc(
                                centre_x + start,
                                centre_y + mid,
                                centre_x + stop,
                                centre_y + mid,
                                radius,
                                pins,
                            )
                        )
                    else:
                        found.append(
                            Disc(
                                centre_x + mid,
                                centre_y + start,
                                centre_x + mid,
                                centre_y + stop,
                                radius,
                                pins,
                            )
                        )

    # Points held by three features.
    for x, y in _pinned_points(near, limit):
        if not inside(x, y):
            continue
        radius = near.clearance(x, y)
        if radius < floor or radius > limit:
            continue
        pins = near.pins(x, y, radius)
        found.append(Disc(centre_x + x, centre_y + y, centre_x + x, centre_y + y, radius, pins))
    return found


def inscribed_circle(walls: Walls, distances: list[float], members: list[int]) -> Disc:
    """The largest circle inside the open ground centred in ``members``.

    Exact.  The circle's centre lies within ``sqrt(2) / 2`` of the centre of
    the cell it is in, and the field value there is within that of the radius,
    so only the cells whose field value is that close to the best one need to
    be searched, and the search around each is the closed-form enumeration in
    :func:`local_circles`.  Where several circles tie, the one found first is
    returned; the radius is what callers rely on.
    """
    top = max(members, key=lambda index: distances[index])
    row, col = divmod(top, walls.ncols)
    best = Disc(col + 0.5, row + 0.5, col + 0.5, row + 0.5, distances[top])
    band = [index for index in members if distances[index] >= best.radius - HALF_DIAGONAL - EPSILON]
    allowed = set(members)
    for index in band:
        for disc in local_circles(walls, distances, index, allowed=allowed):
            if disc.radius > best.radius + EPSILON:
                best = disc
    return best


# ----------------------------------------------------------------------
# the opening function
# ----------------------------------------------------------------------
def _first_unsettled(row_next: list[int], column: int) -> int:
    """First column at or after ``column`` that no disc has settled yet.

    ``row_next`` is a disjoint-set forest over one raster row: a settled column
    points forward, an unsettled one points at itself.  Paths are compressed on
    the way out, so a whole scan costs near-linear time in the number of cells
    rather than in the area of the discs it paints.
    """
    probe = column
    chain: list[int] = []
    while row_next[probe] != probe:
        chain.append(probe)
        probe = row_next[probe]
    for node in chain:
        row_next[node] = probe
    return probe


def paint(sources: list[Disc], nrows: int, ncols: int) -> list[float]:
    """Radius of the largest source covering each cell centre, in cell units.

    Sources are painted largest-first, so the first to reach a cell is the
    largest that ever will and settles it for good.  Cells already settled are
    stepped over through a per-row "next unsettled column" forest instead of
    being rewritten, which bounds the work by the size of the raster rather
    than by the summed area of the discs.  A source that is free to slide
    along a corridor midline paints the whole stadium it sweeps.
    """
    radii = [0.0] * (nrows * ncols)
    nxt = [list(range(ncols + 1)) for _ in range(nrows)]
    for disc in sorted(sources, key=lambda d: -d.radius):
        radius = disc.radius
        first_row = max(0, math.ceil(min(disc.y1, disc.y2) - radius - 0.5 - EPSILON))
        last_row = min(nrows - 1, math.floor(max(disc.y1, disc.y2) + radius - 0.5 + EPSILON))
        for row in range(first_row, last_row + 1):
            centre_y = row + 0.5
            across = max(disc.y1 - centre_y, centre_y - disc.y2, 0.0)
            half = math.sqrt(max(0.0, radius * radius - across * across))
            low = max(0, math.ceil(disc.x1 - half - 0.5 - EPSILON))
            high = min(ncols - 1, math.floor(disc.x2 + half - 0.5 + EPSILON))
            row_next = nxt[row]
            base = row * ncols
            probe = _first_unsettled(row_next, low)
            while probe <= high:
                radii[base + probe] = radius
                row_next[probe] = probe + 1
                probe = _first_unsettled(row_next, probe + 1)
    return radii


def _centre_discs(walls: Walls, distances: list[float]) -> list[Disc]:
    """The cell-centred discs that are not inside a neighbour's.

    Containment is the honest geometric test ``|p - q| + d(p) <= d(q)``,
    checked against the eight adjacent cells only.  It can never drop a disc
    that some cell still needs: whatever the dropped disc covered, the
    containing disc covers with a radius at least as large.
    """
    nrows, ncols = walls.nrows, walls.ncols
    discs: list[Disc] = []
    for index, reach in enumerate(distances):
        if walls.blocked[index] or reach <= 0.0:
            continue
        row, col = divmod(index, ncols)
        for nrow in (row - 1, row, row + 1):
            if not 0 <= nrow < nrows:
                continue
            for ncol in (col - 1, col, col + 1):
                if not 0 <= ncol < ncols or (nrow == row and ncol == col):
                    continue
                step = 1.0 if (nrow == row or ncol == col) else math.sqrt(2.0)
                if distances[nrow * ncols + ncol] >= reach + step - EPSILON:
                    break
            else:
                continue
            break
        else:
            discs.append(Disc(col + 0.5, row + 0.5, col + 0.5, row + 0.5, reach))
    return discs


def _component_maxima(walls: Walls, values: list[float]) -> list[float]:
    """For every cell, the largest ``values`` entry in its 8-connected opening."""
    nrows, ncols = walls.nrows, walls.ncols
    blocked = walls.blocked
    out = [0.0] * len(values)
    seen = [False] * len(values)
    for start in range(len(values)):
        if blocked[start] or seen[start]:
            continue
        seen[start] = True
        members = [start]
        queue: deque[int] = deque([start])
        while queue:
            index = queue.popleft()
            row, col = divmod(index, ncols)
            for nrow in (row - 1, row, row + 1):
                if not 0 <= nrow < nrows:
                    continue
                for ncol in (col - 1, col, col + 1):
                    if not 0 <= ncol < ncols:
                        continue
                    other = nrow * ncols + ncol
                    if seen[other] or blocked[other]:
                        continue
                    seen[other] = True
                    members.append(other)
                    queue.append(other)
        top = max(values[index] for index in members)
        for index in members:
            out[index] = top
    return out


def _ridge_floor(first: Key, second: Key) -> float:
    """A lower bound on the clearance anywhere along the ridge between two features.

    Two corners: half their separation, at the midpoint of their bisector.  A
    corner and a wall line: half the corner's distance from the line, at the
    apex of the parabola between them.  Two lines: the clearance grows away
    from where they cross, so the pinned circles at the ridge's ends bound it
    already, and there is nothing to lower.
    """
    if first[0] == "p" and second[0] == "p":
        return math.hypot(first[1] - second[1], first[2] - second[2]) / 2.0
    if first[0] == "p" or second[0] == "p":
        line, point = (second, first) if first[0] == "p" else (first, second)
        across = point[2] if line[0] == "y" else point[1]
        return abs(across - line[1]) / 2.0
    return math.inf


def _corner_pairs(walls: Walls) -> Iterator[tuple[Key, Key]]:
    """Wall lines that meet end to end, whose ridge starts at that corner with no clearance."""
    ends: dict[tuple[float, float], list[Key]] = {}
    for segment in walls.segments:
        for point in segment.endpoints:
            ends.setdefault(point, []).append(segment.key)
    for keys in ends.values():
        for i, first in enumerate(keys):
            for second in keys[i + 1 :]:
                if first != second:
                    yield first, second


def _at_radius(a: Key, b: Key, radius: float) -> Iterator[tuple[float, float]]:
    """Points exactly ``radius`` from features ``a`` and ``b``, in offsets.

    These are the vertices of the region that a circle of ``radius`` can be
    centred in, where two of its offset walls cross.  Parallel lines are left
    out: the corridor between them is painted directly.
    """
    if a[0] == "p" and b[0] == "p":
        # Two circles of equal radius: on the perpendicular bisector.
        mx, my = (a[1] + b[1]) / 2.0, (a[2] + b[2]) / 2.0
        ux, uy = b[1] - a[1], b[2] - a[2]
        span = math.hypot(ux, uy)
        if span < EPSILON or span > 2.0 * radius + EPSILON:
            return
        height = math.sqrt(max(0.0, radius * radius - span * span / 4.0))
        nx, ny = -uy / span, ux / span
        yield mx + nx * height, my + ny * height
        yield mx - nx * height, my - ny * height
        return
    if a[0] == "p" or b[0] == "p":
        line, point = (b, a) if a[0] == "p" else (a, b)
        horizontal = line[0] == "y"
        along, across = (point[1], point[2]) if horizontal else (point[2], point[1])
        for sign in (-1.0, 1.0):
            level = line[1] + sign * radius
            gap = radius * radius - (level - across) ** 2
            if gap < 0.0:
                continue
            for p in (along - math.sqrt(gap), along + math.sqrt(gap)):
                yield (p, level) if horizontal else (level, p)
        return
    if a[0] == b[0]:
        return
    x_line, y_line = (a[1], b[1]) if a[0] == "x" else (b[1], a[1])
    for sx in (-1.0, 1.0):
        for sy in (-1.0, 1.0):
            yield x_line + sx * radius, y_line + sy * radius


def _feet(near: _Neighbourhood, radius: float) -> Iterator[tuple[float, float]]:
    """Where a circle of ``radius`` touching one wall sits nearest the anchor.

    For every feature closer than ``radius`` to the anchor, the centre pushed
    straight away from it until the circle just touches.  Such a centre is
    within ``radius`` of the anchor by construction; whether the circle fits
    is for the caller to check.
    """
    for horizontal, levels in ((True, near.lines_y), (False, near.lines_x)):
        for level in levels:
            if abs(level) >= radius:
                continue
            away = level - math.copysign(radius, level)
            yield (0.0, away) if horizontal else (away, 0.0)
    for px, py in near.points:
        span = math.hypot(px, py)
        if span >= radius or span < EPSILON:
            continue
        scale = 1.0 - radius / span
        yield px * scale, py * scale


def opening_field(walls: Walls, distances: list[float], *, at_least: float = 0.0) -> list[float]:
    """The opening function at every cell centre, in cell units.

    Cell ``c`` gets the largest ``r`` such that some circle of radius ``r``
    fits inside the open ground and covers the centre of ``c``; blocked cells
    get 0.  Thresholding the result at ``r`` is exactly the morphological
    opening by a disc of that radius.

    The circles that can be largest for some cell are of three kinds, and all
    three are enumerated exactly.  Circles held still by the walls — pinned at
    a point by three features, or free to slide between two parallel ones —
    are found in every open cell by :func:`local_circles` and painted, along
    with the disc centred on every cell.  What remains are circles held by two
    walls *and the cell centre they must cover*: their centres lie on the
    ridge equidistant from those two walls, so only pairs of walls that share
    a ridge are tried, and those are read off the pinned circles — a pinned
    circle's features are pairwise ridge neighbours, and every ridge ends at a
    pinned circle.  Clearance along a ridge peaks at a pinned circle, so a
    ridge whose largest pinned circle is no bigger than what a cell already
    has is not tried for that cell.  A cell is spared the solve altogether
    when the painted circles already give it the most its opening can hold,
    or when no cell whose square could hold a bigger centre lies within reach.

    ``at_least`` asks only whether each cell reaches that radius: circles that
    cannot are not searched for, and the values returned are exact for cells
    at or above it and merely lower bounds below.  That question has a
    cheaper form — is any point within ``at_least`` of the cell centre at
    least ``at_least`` from every wall — whose answer is found where the
    nearest such point can be: pushed straight off one wall, or at a crossing
    of two walls' offsets.
    """
    nrows, ncols = walls.nrows, walls.ncols
    blocked = walls.blocked
    threshold = max(at_least - EPSILON, 0.0)

    centres = _centre_discs(walls, distances)
    sources = list(centres)
    # Every pair of features a pinned circle touches shares a ridge.  The
    # clearance along a ridge peaks at the pinned circles at its ends, but
    # can dip between them: the ridge between two corners is their bisector,
    # lowest at their midpoint, and a corner's ridge with a wall is a
    # parabola, lowest at its apex.  So a ridge's range is the pinned radii
    # at its ends, extended down to that geometric minimum — and to nothing
    # where two walls meet at a convex corner.
    spans: dict[tuple[Key, Key], tuple[float, float]] = {}

    def adjoin(keys: tuple[Key, ...], radius: float) -> None:
        for i, first in enumerate(keys):
            for second in keys[i + 1 :]:
                if first == second:
                    continue
                pair = (first, second) if first < second else (second, first)
                bottom, top = spans.get(pair, (math.inf, -math.inf))
                spans[pair] = (min(bottom, radius), max(top, radius))

    for index in range(len(distances)):
        if blocked[index]:
            continue
        for disc in local_circles(walls, distances, index, at_least=at_least):
            sources.append(disc)
            adjoin(disc.pins, disc.radius)
    for first, second in _corner_pairs(walls):
        pair = (first, second) if first < second else (second, first)
        if pair in spans:
            spans[pair] = (0.0, spans[pair][1])
    adjacency: dict[Key, list[tuple[float, float, Key]]] = {}
    for (first, second), (bottom, top) in spans.items():
        bottom = min(bottom, _ridge_floor(first, second))
        adjacency.setdefault(first, []).append((top, bottom, second))
        adjacency.setdefault(second, []).append((top, bottom, first))
    for ranked in adjacency.values():
        ranked.sort(key=lambda entry: -entry[0])

    values = paint(sources, nrows, ncols)
    ceiling = _component_maxima(walls, values)
    # Where a circle bigger than the painted one could be centred: within
    # ``sqrt(2) / 2`` of a cell centre the clearance is within that of the
    # cell's value, so paint each cell's value plus a half diagonal over a
    # disc that reaches every point whose square it could cover.  The
    # cell-centred discs inside a neighbour's stay inside it after both grow
    # by the same amount, so only the maximal ones need painting.
    reach = paint(
        [
            Disc(disc.x1, disc.y1, disc.x2, disc.y2, disc.radius + 2.0 * HALF_DIAGONAL)
            for disc in centres
        ],
        nrows,
        ncols,
    )

    for index in range(len(distances)):
        if blocked[index]:
            continue
        best = values[index]
        if best >= ceiling[index] - EPSILON:
            continue
        # No circle in this opening is larger than the largest one painted in
        # it, which is exact: the largest is pinned, and pinned circles are
        # painted.
        bound = min(reach[index] - HALF_DIAGONAL, ceiling[index])
        row, col = divmod(index, ncols)
        anchor = (col + 0.5, row + 0.5)
        if at_least > 0.0:
            if best >= threshold or bound < threshold:
                continue
            if _covered(walls, distances, adjacency, anchor, threshold):
                values[index] = threshold
            continue
        if bound <= best + EPSILON:
            continue
        values[index] = _largest_through(walls, distances, adjacency, anchor, best, bound)
    return values


def _offsets(near: _Neighbourhood) -> dict[Key, tuple[Key, float]]:
    """The neighbourhood's features by absolute key, as an offset feature and a bound.

    The bound is half the feature's distance from the anchor, which is the
    smallest circle through the anchor that can reach it: a circle of radius
    ``r`` centred ``r`` from the anchor comes no closer than ``2 r`` to
    anything.
    """
    offsets: dict[Key, tuple[Key, float]] = {}
    for horizontal, levels in ((True, near.lines_y), (False, near.lines_x)):
        for level in levels:
            key = ("y", level) if horizontal else ("x", level)
            offsets[near.key_of_line(horizontal, level)] = (key, abs(level) / 2.0)
    for point in near.points:
        offsets[near.key_of_point(point)] = (
            ("p", point[0], point[1]),
            math.hypot(point[0], point[1]) / 2.0,
        )
    return offsets


def _ridges(
    offsets: dict[Key, tuple[Key, float]],
    adjacency: dict[Key, list[tuple[float, float, Key]]],
    floor: float,
    cap: float,
) -> Iterator[tuple[Key, Key]]:
    """Ridge-sharing pairs among the features whose clearance crosses ``(floor, cap]``.

    Each feature's ridges are ranked by their highest clearance, so the scan
    stops at the first that cannot reach above ``floor`` — raised to what the
    anchor's own distance from each feature already demands.
    """
    for first, (relative_first, reach_first) in offsets.items():
        limit = max(floor, reach_first)
        for top, bottom, second in adjacency.get(first, ()):
            if top <= limit:
                break
            entry = offsets.get(second)
            if entry is None or second <= first or bottom > cap or top <= entry[1]:
                continue
            yield relative_first, entry[0]


def _covered(
    walls: Walls,
    distances: list[float],
    adjacency: dict[Key, list[tuple[float, float, Key]]],
    anchor: tuple[float, float],
    radius: float,
) -> bool:
    """Whether some circle of ``radius`` fits and covers the anchor."""
    anchor_x, anchor_y = anchor
    near = walls.neighbourhood(anchor_x, anchor_y, 2.0 * radius + EPSILON)
    for x, y in _feet(near, radius):
        if _fits(walls, distances, anchor_x, anchor_y, x, y, radius):
            return True
    # A crossing of two walls' offsets at ``radius`` lies on their ridge where
    # the clearance is exactly ``radius``, so only ridges whose clearance
    # passes through that value can hold one.
    for first, second in _ridges(_offsets(near), adjacency, radius - EPSILON, radius + EPSILON):
        for x, y in _at_radius(first, second, radius):
            if math.hypot(x, y) <= radius + EPSILON and _fits(
                walls, distances, anchor_x, anchor_y, x, y, radius
            ):
                return True
    return False


def _largest_through(
    walls: Walls,
    distances: list[float],
    adjacency: dict[Key, list[tuple[float, float, Key]]],
    anchor: tuple[float, float],
    best: float,
    bound: float,
) -> float:
    """The largest circle through the anchor held by two ridge-sharing walls."""
    anchor_x, anchor_y = anchor
    near = walls.neighbourhood(anchor_x, anchor_y, 2.0 * bound + EPSILON)
    for first, second in _ridges(_offsets(near), adjacency, best + EPSILON, bound + EPSILON):
        for x, y in _through_point(first, second, bound):
            radius = math.hypot(x, y)
            if radius <= best + EPSILON or radius > bound + EPSILON:
                continue
            if _fits(walls, distances, anchor_x, anchor_y, x, y, radius):
                best = radius
    return best


def _fits(
    walls: Walls,
    distances: list[float],
    anchor_x: float,
    anchor_y: float,
    dx: float,
    dy: float,
    radius: float,
) -> bool:
    """:meth:`Walls.fits`, with the distance field consulted first.

    Clearance is 1-Lipschitz, so the field value at the centre of the cell the
    point lies in settles most cases without touching a wall: the point's
    clearance is within its distance to that centre of the cell's value.
    """
    x, y = anchor_x + dx, anchor_y + dy
    row, col = math.floor(y), math.floor(x)
    if not (0 <= row < walls.nrows and 0 <= col < walls.ncols):
        return False
    index = row * walls.ncols + col
    if walls.blocked[index]:
        return walls.fits(anchor_x, anchor_y, dx, dy, radius)
    offset = math.hypot(x - (col + 0.5), y - (row + 0.5))
    if distances[index] + offset < radius - EPSILON:
        return False
    if distances[index] - offset >= radius - EPSILON:
        return True
    return walls.fits(anchor_x, anchor_y, dx, dy, radius)

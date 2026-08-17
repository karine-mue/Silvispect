"""Field inventory records, CSV interchange, and detection/field matching."""

from __future__ import annotations

import csv
import io
import math
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, replace
from pathlib import Path

from .detect import Crown

__all__ = [
    "TREE_FIELDS",
    "InventoryError",
    "MatchResult",
    "Plot",
    "Tree",
    "match_trees",
    "parse_trees",
    "read_trees",
    "trees_from_crowns",
    "trees_to_csv",
    "write_trees",
]

#: Canonical column order written by :func:`trees_to_csv`.
TREE_FIELDS: tuple[str, ...] = (
    "tree_id",
    "x",
    "y",
    "species",
    "dbh_cm",
    "height_m",
    "crown_diameter_m",
    "status",
)

#: Accepted column aliases, mapped onto the canonical field names.
COLUMN_ALIASES: dict[str, str] = {
    "id": "tree_id",
    "treeid": "tree_id",
    "tree": "tree_id",
    "stem_id": "tree_id",
    "easting": "x",
    "x_m": "x",
    "northing": "y",
    "y_m": "y",
    "sp": "species",
    "species_code": "species",
    "dbh": "dbh_cm",
    "diameter": "dbh_cm",
    "diameter_cm": "dbh_cm",
    "d": "dbh_cm",
    "height": "height_m",
    "h": "height_m",
    "ht": "height_m",
    "crown_diameter": "crown_diameter_m",
    "cd": "crown_diameter_m",
    "crown_width_m": "crown_diameter_m",
    "state": "status",
}

LIVE_STATUSES = frozenset({"live", "living", "alive", "l", "1"})


class InventoryError(ValueError):
    """Raised for malformed inventory input."""


@dataclass
class Tree:
    """One stem record.

    Coordinates are in the same map units as the rasters (metres by
    convention).  ``dbh_cm`` is the diameter at breast height and ``height_m``
    the total height; either may be missing in real field data, which the
    inspection rules report rather than silently impute.
    """

    tree_id: str
    x: float
    y: float
    species: str = ""
    dbh_cm: float | None = None
    height_m: float | None = None
    crown_diameter_m: float | None = None
    status: str = "live"

    @property
    def is_live(self) -> bool:
        return self.status.strip().lower() in LIVE_STATUSES

    @property
    def basal_area_m2(self) -> float | None:
        """Cross-sectional area at breast height in square metres."""
        if self.dbh_cm is None or self.dbh_cm <= 0:
            return None
        return math.pi * (self.dbh_cm / 200.0) ** 2

    def as_dict(self) -> dict[str, object]:
        return {
            "tree_id": self.tree_id,
            "x": self.x,
            "y": self.y,
            "species": self.species,
            "dbh_cm": self.dbh_cm,
            "height_m": self.height_m,
            "crown_diameter_m": self.crown_diameter_m,
            "status": self.status,
        }


@dataclass
class Plot:
    """A named collection of trees with a known ground area."""

    plot_id: str
    trees: list[Tree]
    area_ha: float

    def __post_init__(self) -> None:
        if self.area_ha <= 0:
            raise InventoryError("plot area must be positive")

    def live_trees(self) -> list[Tree]:
        return [tree for tree in self.trees if tree.is_live]

    def __len__(self) -> int:
        return len(self.trees)


def _normalise_key(key: str) -> str:
    cleaned = key.strip().lower().replace(" ", "_").replace("-", "_")
    cleaned = cleaned.lstrip("﻿")
    return COLUMN_ALIASES.get(cleaned, cleaned)


def _to_float(value: str | None, *, field_name: str, row_number: int) -> float | None:
    if value is None:
        return None
    text = value.strip()
    if text == "" or text.upper() in {"NA", "NAN", "NULL", "-"}:
        return None
    try:
        number = float(text)
    except ValueError as exc:
        raise InventoryError(
            f"row {row_number}: {field_name} value {value!r} is not numeric"
        ) from exc
    if not math.isfinite(number):
        raise InventoryError(
            f"row {row_number}: {field_name} value {value!r} is not a finite number"
        )
    return number


def parse_trees(text: str) -> list[Tree]:
    """Parse a CSV document into :class:`Tree` records.

    Column names are matched case-insensitively and a handful of common
    aliases (``dbh``, ``height``, ``easting`` ...) are accepted.  ``x``, ``y``
    are required; every other column is optional.
    """
    reader = csv.reader(io.StringIO(text))
    try:
        header = next(reader)
    except StopIteration:
        return []
    keys = [_normalise_key(name) for name in header]
    if "x" not in keys or "y" not in keys:
        raise InventoryError("inventory CSV must provide 'x' and 'y' columns")
    trees: list[Tree] = []
    for row_number, row in enumerate(reader, start=2):
        if not any(cell.strip() for cell in row):
            continue
        record = {key: (row[i] if i < len(row) else "") for i, key in enumerate(keys)}
        x = _to_float(record.get("x"), field_name="x", row_number=row_number)
        y = _to_float(record.get("y"), field_name="y", row_number=row_number)
        if x is None or y is None:
            raise InventoryError(f"row {row_number}: missing coordinates")
        tree_id = (record.get("tree_id") or "").strip() or str(len(trees) + 1)
        status = (record.get("status") or "live").strip() or "live"
        trees.append(
            Tree(
                tree_id=tree_id,
                x=x,
                y=y,
                species=(record.get("species") or "").strip(),
                dbh_cm=_to_float(record.get("dbh_cm"), field_name="dbh_cm", row_number=row_number),
                height_m=_to_float(
                    record.get("height_m"), field_name="height_m", row_number=row_number
                ),
                crown_diameter_m=_to_float(
                    record.get("crown_diameter_m"),
                    field_name="crown_diameter_m",
                    row_number=row_number,
                ),
                status=status,
            )
        )
    return trees


def read_trees(path: str | Path) -> list[Tree]:
    return parse_trees(Path(path).read_text(encoding="utf-8"))


def trees_to_csv(trees: Iterable[Tree], *, precision: int = 3) -> str:
    """Serialise tree records to CSV using the canonical column order."""
    buffer = io.StringIO()
    writer = csv.writer(buffer, lineterminator="\n")
    writer.writerow(TREE_FIELDS)
    for tree in trees:
        writer.writerow(
            [
                tree.tree_id,
                round(tree.x, precision),
                round(tree.y, precision),
                tree.species,
                "" if tree.dbh_cm is None else round(tree.dbh_cm, precision),
                "" if tree.height_m is None else round(tree.height_m, precision),
                "" if tree.crown_diameter_m is None else round(tree.crown_diameter_m, precision),
                tree.status,
            ]
        )
    return buffer.getvalue()


def write_trees(trees: Iterable[Tree], path: str | Path, *, precision: int = 3) -> Path:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(trees_to_csv(trees, precision=precision), encoding="utf-8")
    return target


def trees_from_crowns(
    crowns: Sequence[Crown],
    *,
    species: str = "",
    prefix: str = "D",
    dbh_model=None,
) -> list[Tree]:
    """Convert delineated crowns into inventory records.

    Args:
        crowns: Crowns produced by :func:`silvispect.detect.detect_trees`.
        species: Species label applied to every record.
        prefix: Prefix for the generated ``tree_id`` values.
        dbh_model: Optional fitted height-diameter model used to back out a
            diameter from the measured height (see :mod:`silvispect.allometry`).
    """
    trees: list[Tree] = []
    for crown in crowns:
        dbh = None if dbh_model is None else dbh_model.invert(crown.height)
        trees.append(
            Tree(
                tree_id=f"{prefix}{crown.tree_id}",
                x=crown.x,
                y=crown.y,
                species=species,
                dbh_cm=dbh,
                height_m=crown.height,
                crown_diameter_m=crown.diameter,
            )
        )
    return trees


@dataclass(frozen=True)
class MatchResult:
    """Outcome of matching detected crowns against field-measured trees."""

    matches: tuple[tuple[Crown, Tree, float], ...]
    omissions: tuple[Tree, ...]
    commissions: tuple[Crown, ...]
    tolerance: float

    @property
    def matched(self) -> int:
        return len(self.matches)

    @property
    def recall(self) -> float:
        """Share of reference trees that were detected."""
        reference = self.matched + len(self.omissions)
        return self.matched / reference if reference else 0.0

    @property
    def precision(self) -> float:
        """Share of detections that correspond to a reference tree."""
        detected = self.matched + len(self.commissions)
        return self.matched / detected if detected else 0.0

    @property
    def f1(self) -> float:
        denominator = self.recall + self.precision
        if denominator == 0:
            return 0.0
        return 2 * self.recall * self.precision / denominator

    @property
    def height_bias(self) -> float | None:
        """Mean of ``detected - reference`` height over matched pairs."""
        errors = self._height_errors()
        if not errors:
            return None
        return sum(errors) / len(errors)

    @property
    def height_rmse(self) -> float | None:
        errors = self._height_errors()
        if not errors:
            return None
        return math.sqrt(sum(e * e for e in errors) / len(errors))

    @property
    def mean_offset(self) -> float | None:
        """Mean planimetric distance between matched apexes and stems."""
        if not self.matches:
            return None
        return sum(distance for _, _, distance in self.matches) / len(self.matches)

    def _height_errors(self) -> list[float]:
        return [
            crown.height - tree.height_m
            for crown, tree, _ in self.matches
            if tree.height_m is not None
        ]

    def as_dict(self) -> dict[str, object]:
        return {
            "matched": self.matched,
            "omissions": len(self.omissions),
            "commissions": len(self.commissions),
            "tolerance_m": self.tolerance,
            "recall": round(self.recall, 4),
            "precision": round(self.precision, 4),
            "f1": round(self.f1, 4),
            "height_bias_m": _round_or_none(self.height_bias, 3),
            "height_rmse_m": _round_or_none(self.height_rmse, 3),
            "mean_offset_m": _round_or_none(self.mean_offset, 3),
        }


def match_trees(
    crowns: Sequence[Crown],
    reference: Sequence[Tree],
    *,
    tolerance: float = 2.5,
    live_only: bool = True,
) -> MatchResult:
    """Greedily pair crowns with reference trees by planimetric distance.

    Candidate pairs closer than ``tolerance`` are considered in ascending
    distance order and accepted while both partners are still unpaired, which
    yields a stable one-to-one assignment without needing a full Hungarian
    solve.
    """
    if tolerance <= 0:
        raise InventoryError("tolerance must be positive")
    targets = [tree for tree in reference if tree.is_live] if live_only else list(reference)

    candidates: list[tuple[float, int, int]] = []
    for i, crown in enumerate(crowns):
        for j, tree in enumerate(targets):
            distance = math.hypot(crown.x - tree.x, crown.y - tree.y)
            if distance <= tolerance:
                candidates.append((distance, i, j))
    candidates.sort()

    used_crowns: set[int] = set()
    used_trees: set[int] = set()
    matches: list[tuple[Crown, Tree, float]] = []
    for distance, i, j in candidates:
        if i in used_crowns or j in used_trees:
            continue
        used_crowns.add(i)
        used_trees.add(j)
        matches.append((crowns[i], targets[j], distance))

    matches.sort(key=lambda item: item[0].tree_id)
    omissions = tuple(tree for j, tree in enumerate(targets) if j not in used_trees)
    commissions = tuple(crown for i, crown in enumerate(crowns) if i not in used_crowns)
    return MatchResult(tuple(matches), omissions, commissions, tolerance)


def shift_trees(trees: Iterable[Tree], dx: float, dy: float) -> list[Tree]:
    """Return copies of ``trees`` translated by ``(dx, dy)``."""
    return [replace(tree, x=tree.x + dx, y=tree.y + dy) for tree in trees]


def _round_or_none(value: float | None, digits: int) -> float | None:
    return None if value is None else round(value, digits)

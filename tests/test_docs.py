"""The reference documentation describes the behaviour that is implemented."""

from __future__ import annotations

from pathlib import Path

DOCS = Path(__file__).resolve().parent.parent / "docs"


def _section(text: str, start: str, end: str) -> str:
    begin = text.index(start)
    return text[begin : text.index(end, begin)]


def test_the_detector_tie_break_has_one_current_description():
    """Two superseded rules were left in place beside the current one.

    A reader met the distance-height profile and its row-and-column fallback
    as the rule, then the canonical orientation as the rule, and could not
    tell which was implemented.  History may stay, labelled as history; the
    normative text names one mechanism.
    """
    text = (DOCS / "concepts.md").read_text()
    normative = _section(text, "### Variable-window local maxima", "*Superseded.*")
    assert "canonical orientation" in normative
    for stale in ("(distance, height)", "row-major", "row and column"):
        assert stale not in normative, stale
    history = _section(text, "*Superseded.*", "The CHM is mean-filtered first")
    assert "row-major" in history and "(distance, height)" in history
    assert "Neither is used" in history


def test_the_gap_width_documentation_describes_exact_geometry():
    """The chamfer walk is history, and is described as such."""
    text = (DOCS / "concepts.md").read_text()
    current = _section(
        text, "The field is the union of every circle that fits.", "Two earlier forms"
    )
    assert "exactly" in current and "lattice corner" in current
    assert "the chamfer carries on" not in text
    assert "chamfer-shortest" not in text

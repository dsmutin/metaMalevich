"""Mandatory checks for length-weighted profiles."""

from __future__ import annotations

import pytest

from metamalevich.reprofile import profile_nodes
from metamalevich.taxonomy import parse_kraken_report

pytestmark = pytest.mark.mandatory

REPORT = """\
100.00\t2\t0\tR\t1\troot
100.00\t2\t0\tS\t10\t  speciesA
"""


def test_unclassified_mass_is_counted_once() -> None:
    """Taxon 0 and the residual probability are the same unclassified bases."""
    taxonomy = parse_kraken_report(REPORT, source="toy", version="toy")
    rows, summary = profile_nodes({"n": {0: 0.25, 10: 0.75}}, {"n": 100}, taxonomy)
    assert summary["unclassified_bases"] == pytest.approx(25)
    assert summary["unclassified_fraction"] == pytest.approx(0.25)
    assert rows[0]["assigned_bases"] == pytest.approx(75)
    assert rows[0]["relative_abundance"] == pytest.approx(0.75)


def test_empty_distribution_is_unclassified() -> None:
    """A node with no colours contributes its whole length once."""
    taxonomy = parse_kraken_report(REPORT, source="toy", version="toy")
    rows, summary = profile_nodes({"n": {}}, {"n": 40}, taxonomy)
    assert rows == []
    assert summary["unclassified_bases"] == pytest.approx(40)

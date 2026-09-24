"""Mandatory checks for cutting confident disagreements before leakage."""

from __future__ import annotations

import pytest

from metamalevich.edge_cut import cut_then_leak
from metamalevich.leakage import decaying_leakage

pytestmark = pytest.mark.mandatory


def test_confident_disagreement_stays_one_hot() -> None:
    """A confident disagreeing edge is dropped, so one-hot nodes stay one-hot."""
    distributions = {"a": {1: 1.0}, "b": {2: 1.0}}
    edges = [
        {"source": "a", "target": "b"},
        {"source": "b", "target": "a"},
    ]
    result = cut_then_leak(distributions, edges, taxonomy=None, confident=0.8, decay=0.5, iterations=4)
    assert result["a"] == {1: 1.0}
    assert result["b"] == {2: 1.0}
    leaked = decaying_leakage(distributions, edges, decay=0.5, iterations=4)
    assert leaked["a"] != {1: 1.0}


def test_uncertain_edge_still_mixes() -> None:
    """An edge between nodes below the confidence cutoff still leaks."""
    distributions = {"a": {1: 0.6, 2: 0.4}, "b": {2: 0.6, 1: 0.4}}
    edges = [
        {"source": "a", "target": "b"},
        {"source": "b", "target": "a"},
    ]
    result = cut_then_leak(distributions, edges, confident=0.8, decay=0.5, iterations=4)
    expected = decaying_leakage(distributions, edges, decay=0.5, iterations=4)
    assert result == expected
    assert result["a"][1] < 0.6
    assert result["a"][2] > 0.4


def test_input_edge_list_is_unchanged() -> None:
    """Cutting copies the kept edges and leaves the caller's list in place."""
    distributions = {"a": {1: 1.0}, "b": {2: 1.0}}
    edges = [{"source": "a", "target": "b", "weight": 3.0}]
    before = [dict(edge) for edge in edges]
    cut_then_leak(distributions, edges)
    assert edges == before
    assert distributions == {"a": {1: 1.0}, "b": {2: 1.0}}

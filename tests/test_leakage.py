"""Mandatory checks for decaying label leakage."""

from __future__ import annotations

import pytest

from metamalevich.aggregate import normalize
from metamalevich.leakage import decaying_leakage

pytestmark = pytest.mark.mandatory

_LINK = (
    {"source": "confident", "target": "wrong"},
    {"source": "wrong", "target": "confident"},
)


def test_decay_zero_returns_renormalized_input() -> None:
    """A decay of 0 keeps the restart distribution and does not edit the input."""
    distributions = {"a": {1: 2.0, 2: 2.0}, "b": {3: 4.0, 4: 0.0}}
    original = {"a": {1: 2.0, 2: 2.0}, "b": {3: 4.0, 4: 0.0}}
    edges = [{"source": "a", "target": "b", "weight": 2.0}]
    result = decaying_leakage(distributions, edges, decay=0, iterations=3)
    assert result["a"] == normalize(distributions["a"])
    assert result["b"] == {3: 1.0}
    assert set(result) == {"a", "b"}
    assert distributions == original


def test_restart_keeps_confident_label() -> None:
    """Restart leaves the confident taxon in front, with leaked mass below the decay."""
    distributions = {"confident": {1: 1.0}, "wrong": {2: 1.0}}
    result = decaying_leakage(distributions, list(_LINK), decay=0.5, iterations=6)
    confident = result["confident"]
    assert confident[1] > confident[2]
    assert 0.0 < confident[2] < 0.5
    assert distributions["confident"] == {1: 1.0}
    assert distributions["wrong"] == {2: 1.0}


def test_two_iterations_move_mass_without_unbounded_leakage() -> None:
    """Two iterations leak label mass, and further iterations stay below the decay."""
    distributions = {"confident": {1: 1.0}, "wrong": {2: 1.0}}
    decay = 0.5
    two = decaying_leakage(distributions, list(_LINK), decay=decay, iterations=2)
    later = decaying_leakage(distributions, list(_LINK), decay=decay, iterations=8)
    assert two["confident"] != {1: 1.0}
    assert 0.0 < two["confident"][2] < decay
    assert 0.0 < later["confident"][2] < decay


def test_invalid_decay_raises() -> None:
    """Decay outside [0, 1] is rejected."""
    for decay in (-0.01, 1.01):
        with pytest.raises(ValueError):
            decaying_leakage({"a": {1: 1.0}}, [], decay=decay)


def test_disconnected_node_is_unchanged() -> None:
    """A node with no positive-weight neighbours keeps its renormalized evidence."""
    distributions = {
        "alone": {1: 1.0},
        "confident": {1: 1.0},
        "wrong": {2: 1.0},
    }
    edges = [
        {"source": "confident", "target": "wrong", "weight": 1.0},
        {"source": "wrong", "target": "confident", "weight": 1.0},
        {"source": "alone", "target": "wrong", "weight": 0.0},
    ]
    result = decaying_leakage(distributions, edges, decay=0.5, iterations=4)
    assert result["alone"] == {1: 1.0}
    assert result["confident"][2] > 0.0
    assert set(result) == set(distributions)

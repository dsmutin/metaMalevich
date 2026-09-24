"""Mandatory checks for locking confident calls against neighbour flips."""

from __future__ import annotations

import copy

import pytest

from metamalevich.confident_lock import confident_lock
from metamalevich.taxonomy import Taxon, Taxonomy

pytestmark = pytest.mark.mandatory


def _taxonomy() -> Taxonomy:
    """Two genera, each with two species."""
    records = [
        Taxon(1, None, "R", "root", "toy", "v"),
        Taxon(10, 1, "G", "genus-a", "toy", "v"),
        Taxon(11, 10, "S", "species-a1", "toy", "v"),
        Taxon(12, 10, "S", "species-a2", "toy", "v"),
        Taxon(20, 1, "G", "genus-b", "toy", "v"),
        Taxon(21, 20, "S", "species-b1", "toy", "v"),
    ]
    return Taxonomy(records, source="toy", version="v")


def test_confident_node_is_unchanged_when_the_only_neighbour_disagrees() -> None:
    """A call at or above the cutoff stays put even if its only neighbour disagrees."""
    taxonomy = _taxonomy()
    distributions = {"a": {11: 0.9, 12: 0.1}, "b": {12: 1.0}}
    edges = [{"source": "a", "target": "b", "weight": 1.0}]
    locked = confident_lock(distributions, edges, taxonomy)
    assert locked["a"] == {11: 0.9, 12: 0.1}
    assert locked["a"] is not distributions["a"]
    assert locked["b"] == {12: 1.0}


def test_uncertain_node_adopts_unanimous_same_genus_neighbour() -> None:
    """An uncertain node takes a strong neighbour vote that stays in its genus."""
    taxonomy = _taxonomy()
    distributions = {"a": {11: 0.6, 12: 0.4}, "b": {12: 1.0}, "c": {12: 1.0}}
    edges = [
        {"source": "a", "target": "b", "weight": 1.0},
        {"source": "a", "target": "c", "weight": 1.0},
    ]
    locked = confident_lock(distributions, edges, taxonomy)
    assert locked["a"] == {12: 1.0}
    assert locked["b"] == {12: 1.0}
    assert locked["c"] == {12: 1.0}


def test_different_genus_neighbour_is_ignored() -> None:
    """A strong vote from another genus does not replace an uncertain node."""
    taxonomy = _taxonomy()
    distributions = {"a": {11: 0.6, 12: 0.4}, "b": {21: 1.0}}
    edges = [{"source": "a", "target": "b", "weight": 1.0}]
    locked = confident_lock(distributions, edges, taxonomy)
    assert locked["a"] == {11: 0.6, 12: 0.4}


def test_inputs_are_not_mutated() -> None:
    """Distributions and edges are unchanged, and the result does not alias them."""
    taxonomy = _taxonomy()
    distributions = {"a": {11: 0.6, 12: 0.4}, "b": {12: 1.0}}
    edges = [{"source": "a", "target": "b", "weight": 1.0}]
    distributions_before = copy.deepcopy(distributions)
    edges_before = copy.deepcopy(edges)
    locked = confident_lock(distributions, edges, taxonomy)
    assert distributions == distributions_before
    assert edges == edges_before
    locked["a"][12] = 0.0
    locked["b"][12] = 0.0
    assert distributions == distributions_before
    assert edges == edges_before


def test_thresholds_outside_unit_interval_raise() -> None:
    """confident and vote_min must lie in [0, 1]. The endpoints are accepted."""
    taxonomy = _taxonomy()
    distributions = {"a": {11: 1.0}}
    for confident in (-0.01, 1.01):
        with pytest.raises(ValueError):
            confident_lock(distributions, [], taxonomy, confident=confident)
    for vote_min in (-0.01, 1.01):
        with pytest.raises(ValueError):
            confident_lock(distributions, [], taxonomy, vote_min=vote_min)
    assert confident_lock(distributions, [], taxonomy, confident=0.0, vote_min=0.0)["a"] == {11: 1.0}
    assert confident_lock(distributions, [], taxonomy, confident=1.0, vote_min=1.0)["a"] == {11: 1.0}

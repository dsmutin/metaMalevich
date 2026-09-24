"""Mandatory checks for unanimous neighbour rescue."""

from __future__ import annotations

import copy

import pytest

from metamalevich.taxonomy import Taxonomy, Taxon
from metamalevich.unanimous_rescue import unanimous_rescue

pytestmark = pytest.mark.mandatory


def _taxonomy() -> Taxonomy:
    """Two genera, each with one species pair used by the rescue checks."""
    return Taxonomy(
        [
            Taxon(1, None, "R", "root", "toy", "v"),
            Taxon(10, 1, "G", "GenusA", "toy", "v"),
            Taxon(11, 10, "S", "speciesA", "toy", "v"),
            Taxon(12, 10, "S", "speciesB", "toy", "v"),
            Taxon(20, 1, "G", "GenusB", "toy", "v"),
            Taxon(21, 20, "S", "speciesC", "toy", "v"),
        ],
        source="toy",
        version="v",
    )


def _scene() -> tuple[Taxonomy, dict[str, dict[int, float]], list[dict]]:
    """Uncertain, confident, disputed, and cross-genus nodes with two neighbours."""
    distributions = {
        "uncertain": {11: 0.5, 12: 0.25, 21: 0.25},
        "confident": {11: 0.75, 12: 0.25},
        "disputed": {11: 0.5, 12: 0.25, 21: 0.25},
        "foreign": {21: 0.5, 11: 0.25, 12: 0.25},
        "blank": {0: 0.5, 12: 0.25, 21: 0.25},
        "isolated": {11: 0.5, 12: 0.25, 21: 0.25},
        "left": {12: 0.75, 11: 0.25},
        "right": {12: 1.0},
        "other": {21: 1.0},
    }
    edges = [
        {"source": "uncertain", "target": "left", "weight": 1.0},
        {"source": "uncertain", "target": "right", "weight": 2.0},
        {"source": "confident", "target": "left", "weight": 1.0},
        {"source": "confident", "target": "right", "weight": 1.0},
        {"source": "disputed", "target": "left", "weight": 1.0},
        {"source": "disputed", "target": "other", "weight": 1.0},
        {"source": "foreign", "target": "left", "weight": 1.0},
        {"source": "foreign", "target": "right", "weight": 1.0},
        {"source": "blank", "target": "left", "weight": 1.0},
        {"source": "blank", "target": "right", "weight": 1.0},
        {"source": "isolated", "target": "left", "weight": 0.0},
        {"source": "isolated", "target": "right", "weight": -1.0},
    ]
    return _taxonomy(), distributions, edges


def test_uncertain_unanimous_neighbours_become_a_mixture() -> None:
    """An uncertain node with two neighbours on one taxon keeps that taxon above 0.5."""
    taxonomy, distributions, edges = _scene()
    rescued = unanimous_rescue(distributions, edges, taxonomy)["uncertain"]
    assert rescued[12] > 0.5
    assert rescued[11] == pytest.approx(0.25)
    assert rescued[12] == pytest.approx(0.625)
    assert rescued[21] == pytest.approx(0.125)


def test_confident_node_is_unchanged() -> None:
    """A node at or above the cutoff keeps its own distribution."""
    taxonomy, distributions, edges = _scene()
    result = unanimous_rescue(distributions, edges, taxonomy)
    assert result["confident"] == {11: 0.75, 12: 0.25}
    assert result["confident"] is not distributions["confident"]


def test_disagreeing_neighbours_leave_the_node_unchanged() -> None:
    """Neighbours that name different taxa do not change the node."""
    taxonomy, distributions, edges = _scene()
    result = unanimous_rescue(distributions, edges, taxonomy)
    assert result["disputed"] == {11: 0.5, 12: 0.25, 21: 0.25}


def test_inputs_are_not_mutated() -> None:
    """Distributions and edges are unchanged, and the result does not alias them."""
    taxonomy, distributions, edges = _scene()
    distributions_before = copy.deepcopy(distributions)
    edges_before = copy.deepcopy(edges)
    result = unanimous_rescue(distributions, edges, taxonomy)
    assert distributions == distributions_before
    assert edges == edges_before
    assert result["uncertain"] is not distributions["uncertain"]
    result["uncertain"][12] = 0.0
    assert distributions["uncertain"][12] == pytest.approx(0.25)
    assert result["confident"] is not distributions["confident"]


def test_uncertain_outside_unit_interval_raises() -> None:
    """Values outside [0, 1] are rejected."""
    taxonomy, distributions, edges = _scene()
    for cutoff in (-0.01, 1.01):
        with pytest.raises(ValueError):
            unanimous_rescue(distributions, edges, taxonomy, uncertain=cutoff)
    unanimous_rescue(distributions, edges, taxonomy, uncertain=0.0)
    unanimous_rescue(distributions, edges, taxonomy, uncertain=1.0)


def test_different_genus_is_not_adopted() -> None:
    """Unanimous neighbours in another genus leave the node unchanged."""
    taxonomy, distributions, edges = _scene()
    result = unanimous_rescue(distributions, edges, taxonomy)
    assert result["foreign"] == {21: 0.5, 11: 0.25, 12: 0.25}


def test_unclassified_node_adopts_a_neighbour_genus() -> None:
    """An uncertain unclassified node adopts a neighbour that has a genus."""
    taxonomy, distributions, edges = _scene()
    rescued = unanimous_rescue(distributions, edges, taxonomy)["blank"]
    assert rescued[12] > 0.5
    assert rescued[12] == pytest.approx(0.625)


def test_node_without_positive_outgoing_neighbour_is_unchanged() -> None:
    """Zero and negative edge weights do not count as neighbours."""
    taxonomy, distributions, edges = _scene()
    result = unanimous_rescue(distributions, edges, taxonomy)
    assert result["isolated"] == {11: 0.5, 12: 0.25, 21: 0.25}

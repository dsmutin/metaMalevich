"""Mandatory checks for the contig line-graph flip."""

from __future__ import annotations

import copy

import pytest

from metamalevich.flip import flip_graph, infer_on_flipped_graph, project_to_nodes

pytestmark = pytest.mark.mandatory


def test_flipped_node_sums_endpoint_taxa() -> None:
    """An edge node keeps every taxon on either endpoint, summed rather than min-ed."""
    edges = [
        {"edge_id": "e1", "source": "A", "target": "B", "orientation": "+"},
    ]
    nodes = {"A": {1: 0.8, 2: 0.2}, "B": {2: 1.0}}
    flipped, flipped_edges = flip_graph(edges, nodes)
    assert flipped_edges == []
    assert set(flipped) == {"e1"}
    assert set(flipped["e1"]) == {1, 2}
    assert flipped["e1"][2] > flipped["e1"][1]
    assert flipped["e1"][1] == pytest.approx(0.4)
    assert flipped["e1"][2] == pytest.approx(0.6)


def test_path_produces_one_flipped_edge() -> None:
    """A path of two contig edges shares one vertex and becomes one line edge."""
    edges = [
        {"edge_id": "e1", "source": "A", "target": "B", "orientation": "++", "weight": 0.4},
        {"edge_id": "e2", "source": "B", "target": "C", "orientation": "++", "weight": 0.9},
    ]
    nodes = {"A": {1: 1.0}, "B": {1: 1.0}, "C": {2: 1.0}}
    _flipped, flipped_edges = flip_graph(edges, nodes)
    assert len(flipped_edges) == 1
    edge = flipped_edges[0]
    assert edge == {
        "edge_id": "line:e1:e2",
        "source": "e1",
        "target": "e2",
        "orientation": "+",
        "weight": pytest.approx(0.4),
    }
    _reversed, reversed_edges = flip_graph(list(reversed(edges)), nodes)
    assert reversed_edges == flipped_edges


def test_project_to_nodes_puts_distribution_on_both_endpoints() -> None:
    """Both endpoints of an edge receive that edge's flipped distribution."""
    edges = [
        {"edge_id": "e1", "source": "A", "target": "B", "orientation": "+", "weight": 1.0},
    ]
    flipped = {"e1": {1: 0.4, 2: 0.6}}
    fallback = {"A": {1: 1.0}, "B": {2: 1.0}}
    projected = project_to_nodes(edges, flipped, fallback)
    assert set(projected) == {"A", "B"}
    assert projected["A"][1] == pytest.approx(0.4)
    assert projected["A"][2] == pytest.approx(0.6)
    assert projected["B"][1] == pytest.approx(0.4)
    assert projected["B"][2] == pytest.approx(0.6)


def test_infer_identity_resolver_returns_original_node_keys() -> None:
    """A resolver that echoes its input still labels the original contigs."""
    edges = [
        {"edge_id": "e1", "source": "A", "target": "B", "orientation": "+", "weight": 1.0},
    ]
    nodes = {"A": {1: 0.8, 2: 0.2}, "B": {2: 1.0}, "Z": {0: 1.0}}
    calls = {"n": 0}

    def resolver(flipped_distributions, flipped_edges):
        calls["n"] += 1
        assert set(flipped_distributions) == {"e1"}
        assert flipped_edges == []
        return flipped_distributions

    result = infer_on_flipped_graph(edges, nodes, resolver)
    assert calls["n"] == 1
    assert set(result) == {"A", "B", "Z"}
    assert result["A"][1] == pytest.approx(0.4)
    assert result["A"][2] == pytest.approx(0.6)
    assert result["B"][1] == pytest.approx(0.4)
    assert result["B"][2] == pytest.approx(0.6)
    assert result["Z"] == {0: 1.0}


def test_inputs_are_not_mutated() -> None:
    """Graph edges and distributions are unchanged after flip, project, and infer."""
    edges = [
        {"edge_id": "e1", "source": "A", "target": "B", "orientation": "+", "weight": 1.0},
        {"edge_id": "e2", "source": "B", "target": "C", "orientation": "+", "weight": 0.5},
    ]
    nodes = {"A": {1: 0.8, 2: 0.2}, "B": {2: 1.0}, "C": {1: 0.3, 3: 0.7}}
    flipped = {"e1": {1: 0.4, 2: 0.6}, "e2": {2: 1.0}}
    edges_before = copy.deepcopy(edges)
    nodes_before = copy.deepcopy(nodes)
    flipped_before = copy.deepcopy(flipped)

    flip_graph(edges, nodes)
    projected = project_to_nodes(edges, flipped, nodes)
    infer_on_flipped_graph(edges, nodes, lambda distributions, _flipped_edges: distributions)

    assert edges == edges_before
    assert nodes == nodes_before
    assert flipped == flipped_before
    projected["A"][1] = 0.0
    assert nodes["A"][1] == pytest.approx(0.8)

    isolated = {"Z": {0: 1.0}}
    isolated_before = copy.deepcopy(isolated)
    kept = project_to_nodes([], isolated, isolated)
    kept["Z"][0] = 0.0
    assert isolated == isolated_before
    assert kept["Z"] is not isolated["Z"]

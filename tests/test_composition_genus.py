"""Mandatory checks for genus transfer on the 4-mer graph."""

from __future__ import annotations

import pytest

from metamalevich.composition_genus import composition_genus_graph

pytestmark = pytest.mark.mandatory


def test_unlabelled_node_copies_the_identical_donor() -> None:
    """A contig with the same sequence as a labelled donor receives that genus."""
    sequence = "ACGTACGTACGT"
    labels = {"donor": 10, "query": 0}
    sequences = {"donor": sequence, "query": sequence}
    transferred = composition_genus_graph(sequences, labels)
    assert transferred["query"] == 10
    assert transferred["donor"] == 10
    assert labels["query"] == 0


def test_labelled_node_is_not_replaced() -> None:
    """A positive label stays even when another donor is present."""
    transferred = composition_genus_graph(
        {"a": "AAAAAAAAAAAA", "b": "CCCCCCCCCCCC"},
        {"a": 1, "b": 2},
    )
    assert transferred == {"a": 1, "b": 2}


def test_no_donor_leaves_the_query_unlabelled() -> None:
    """Nothing is invented when every label is missing."""
    transferred = composition_genus_graph({"a": "ACGTACGTACGT"}, {"a": 0})
    assert transferred["a"] == 0


def test_min_sim_can_refuse_a_weak_neighbour() -> None:
    """A similarity floor leaves a distant query unlabelled."""
    transferred = composition_genus_graph(
        {"donor": "AAAAAAAAAAAAAAAA", "query": "CCCCCCCCCCCCCCCC"},
        {"donor": 5, "query": 0},
        min_sim=0.99,
    )
    assert transferred["query"] == 0


def test_min_sim_outside_unit_interval_raises() -> None:
    """The cosine floor is a similarity, not an arbitrary weight."""
    with pytest.raises(ValueError):
        composition_genus_graph({"a": "ACGT"}, {"a": 1}, min_sim=1.5)

"""Mandatory scores for classification and abundance."""

from __future__ import annotations

from pathlib import Path

import pytest

from metamalevich.evaluate import abundance_scores, classification_scores, neighbour_agreement

pytestmark = pytest.mark.mandatory


def test_neighbour_agreement_counts_matching_endpoints() -> None:
    """Agreement is the share of edges whose classified endpoints match."""
    edges = [{"source": "a", "target": "b"}, {"source": "b", "target": "c"}]
    predicted = {"a": {1: 1.0}, "b": {1: 1.0}, "c": {2: 1.0}}
    assert neighbour_agreement(edges, predicted) == pytest.approx(0.5)


def test_unpredicted_truth_taxon_scores_zero() -> None:
    """A species that is never predicted pulls the macro F1 down to include a zero."""
    predicted = {"a": {1: 1.0}, "b": {1: 1.0}}
    truth = {"a": 1, "b": 2}
    scores = classification_scores(predicted, truth)
    assert scores["precision"] == pytest.approx(0.25)
    assert scores["recall"] == pytest.approx(0.5)
    assert scores["f1"] == pytest.approx(1.0 / 3.0)
    assert scores["accuracy"] == pytest.approx(0.5)


def test_beats_initial_requires_both_l1_and_f1() -> None:
    """Lower L1 with a worse macro F1 does not count as beating the hard colouring."""
    from metamalevich.bench import ABUNDANCE_CHART_HYPOTHESES, compare_to_initial

    rows = [
        {"dataset": "d", "hypothesis": "initial_colouring", "l1": 0.2, "f1": 0.9},
        {"dataset": "d", "hypothesis": "probability_sum", "l1": 0.1, "f1": 0.95},
        {"dataset": "d", "hypothesis": "gated_neighbour", "l1": 0.05, "f1": 0.8},
    ]
    by_name = {row["hypothesis"]: row for row in compare_to_initial(rows)}
    assert by_name["probability_sum"]["beats_initial"] is True
    assert by_name["gated_neighbour"]["beats_initial"] is False
    assert by_name["gated_neighbour"]["beats_initial_l1"] is True
    assert "probability_sum" in ABUNDANCE_CHART_HYPOTHESES


def test_graph_provenance_names_the_4mer_cosine() -> None:
    """Manifest text names the graph that was actually built."""
    from metamalevich.bench import GRAPH_CONSTRUCTION_METHOD, INPUT_GRAPH
    from metamalevich.bridge import _display_path

    assert "4-mer" in INPUT_GRAPH
    assert "Jaccard" not in INPUT_GRAPH
    assert "4-mer" in GRAPH_CONSTRUCTION_METHOD
    assert _display_path(Path("/repo/intermediate/cfa"), Path("/repo")) == "intermediate/cfa"


def test_abundance_keeps_unclassified_mass() -> None:
    """Predicted taxon 0 stays in the L1 total instead of being rescaled away."""
    scores = abundance_scores({1: 0.5, 0: 0.5}, {1: 1.0})
    assert scores["l1"] == pytest.approx(1.0)

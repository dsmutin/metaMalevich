"""Mandatory scores for classification and abundance."""

from __future__ import annotations

from pathlib import Path

import pytest

from metamalevich.evaluate import abundance_scores, classification_scores

pytestmark = pytest.mark.mandatory


def test_unpredicted_truth_taxon_scores_zero() -> None:
    """A species that is never predicted pulls the macro F1 down to include a zero."""
    predicted = {"a": {1: 1.0}, "b": {1: 1.0}}
    truth = {"a": 1, "b": 2}
    scores = classification_scores(predicted, truth)
    assert scores["precision"] == pytest.approx(0.25)
    assert scores["recall"] == pytest.approx(0.5)
    assert scores["f1"] == pytest.approx(1.0 / 3.0)
    assert scores["accuracy"] == pytest.approx(0.5)


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

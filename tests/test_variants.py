"""Mandatory checks that new graph hypotheses stay on the original nodes."""

from __future__ import annotations

import pytest

from metamalevich.bench import HYPOTHESES, decision_variants, graph_variants
from metamalevich.resolve import argmax_taxon
from metamalevich.taxonomy import Taxon, Taxonomy

pytestmark = pytest.mark.mandatory

EDGES = [
    {"edge_id": "e1", "source": "a", "target": "b", "orientation": "+", "weight": 1.0},
    {"edge_id": "e2", "source": "b", "target": "c", "orientation": "+", "weight": 1.0},
]
EVIDENCE = {
    "a": {1: 1.0},
    "b": {1: 0.6, 2: 0.4},
    "c": {2: 1.0},
}


def test_graph_variants_cover_flip_leakage_and_drops() -> None:
    """Every new hypothesis is scored, and each one returns the contig ids."""
    expected = {
        "edge_union",
        "leakage",
        "leakage_flipped",
        "label_drop",
        "label_drop_flipped",
        "logistic_drop",
        "logistic_drop_flipped",
    }
    assert expected <= set(HYPOTHESES)
    variants = graph_variants(EVIDENCE, EDGES)
    assert set(variants) == expected
    for name, distributions in variants.items():
        assert set(distributions) == {"a", "b", "c"}, name


def test_decision_variants_stay_on_original_nodes() -> None:
    """The five error-split resolvers return every original contig id."""
    taxonomy = Taxonomy(
        [
            Taxon(1, None, "G", "Genus", "test", "v"),
            Taxon(2, 1, "S", "Species A", "test", "v"),
            Taxon(3, 1, "S", "Species B", "test", "v"),
        ],
        source="test",
        version="v",
    )
    expected = {
        "confident_lock",
        "unanimous_rescue",
        "cut_then_leak",
        "genus_plurality",
        "margin_mixture",
    }
    assert expected <= set(HYPOTHESES)
    variants = decision_variants(
        EVIDENCE,
        EDGES,
        taxonomy,
        {"a": 1, "b": 1, "c": 2},
        {"a": 10, "b": 10, "c": 10},
    )
    assert set(variants) == expected
    for name, distributions in variants.items():
        assert set(distributions) == {"a", "b", "c"}, name


def test_edge_union_keeps_taxa_from_both_endpoints() -> None:
    """Projecting the edge union onto contig a admits b's taxon."""
    projected = graph_variants(EVIDENCE, EDGES)["edge_union"]
    assert 2 in projected["a"]
    assert argmax_taxon(projected["a"])[0] == 1

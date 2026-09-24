"""Mandatory checks for the genus-plurality prior."""

from __future__ import annotations

import pytest

from metamalevich.genus_plurality import genus_plurality
from metamalevich.taxonomy import Taxonomy, Taxon

pytestmark = pytest.mark.mandatory


def _taxonomy() -> Taxonomy:
    """Species parents point at genus. Ranks are Kraken codes G and S."""
    return Taxonomy(
        [
            Taxon(10, None, "G", "GenusA", "toy", "v"),
            Taxon(11, 10, "S", "SpeciesA1", "toy", "v"),
            Taxon(12, 10, "S", "SpeciesA2", "toy", "v"),
            Taxon(20, None, "G", "GenusB", "toy", "v"),
            Taxon(21, 20, "S", "SpeciesB1", "toy", "v"),
        ],
        source="toy",
        version="v",
    )


def test_uncertain_node_adopts_the_single_confident_donor() -> None:
    """An uncertain same-genus node becomes a one-hot of the confident donor."""
    taxonomy = _taxonomy()
    distributions = {"donor": {11: 0.95, 12: 0.05}, "uncertain": {12: 0.4, 11: 0.2}}
    edges = [{"source": "uncertain", "target": "donor", "weight": 1.0}]
    resolved = genus_plurality(distributions, edges, taxonomy)
    assert resolved["uncertain"] == {11: 1.0}
    assert resolved["donor"] == {11: 0.95, 12: 0.05}


def test_confident_node_is_not_replaced() -> None:
    """A confident node keeps its distribution when another species has more length."""
    taxonomy = _taxonomy()
    distributions = {"short": {12: 0.99}, "long": {11: 0.99}}
    resolved = genus_plurality(distributions, [], taxonomy, lengths={"short": 1, "long": 50})
    assert resolved["short"] == {12: 0.99}
    assert resolved["long"] == {11: 0.99}


def test_genus_without_a_donor_is_unchanged() -> None:
    """A genus with only an uncertain node contributes no plurality."""
    taxonomy = _taxonomy()
    distributions = {"alone": {21: 0.4}, "donor": {11: 1.0}}
    resolved = genus_plurality(distributions, [], taxonomy)
    assert resolved["alone"] == {21: 0.4}


def test_taxon_zero_uses_calls_to_find_the_genus() -> None:
    """An unclassified node takes its genus from the raw call, including genus rank."""
    taxonomy = _taxonomy()
    distributions = {"donor": {11: 0.95}, "blank": {0: 1.0}}
    resolved = genus_plurality(distributions, [], taxonomy, calls={"blank": 10})
    assert resolved["blank"] == {11: 1.0}


def test_inputs_are_not_mutated() -> None:
    """Distributions, edges, calls, and lengths are unchanged after resolution."""
    taxonomy = _taxonomy()
    distributions = {"donor": {11: 0.95, 12: 0.05}, "uncertain": {12: 0.2}, "blank": {0: 1.0}}
    edges = [{"source": "uncertain", "target": "blank"}]
    calls = {"blank": 12}
    lengths = {"donor": 4}
    original_distributions = {node: dict(weights) for node, weights in distributions.items()}
    original_edges = [dict(edge) for edge in edges]
    original_calls = dict(calls)
    original_lengths = dict(lengths)
    resolved = genus_plurality(distributions, edges, taxonomy, calls=calls, lengths=lengths)
    assert distributions == original_distributions
    assert edges == original_edges
    assert calls == original_calls
    assert lengths == original_lengths
    assert resolved["uncertain"] == {11: 1.0}
    assert resolved["donor"] is not distributions["donor"]


def test_plurality_sums_lengths_and_breaks_ties_by_smaller_taxon_id() -> None:
    """Summed donor length picks the species. An equal sum takes the smaller id."""
    taxonomy = _taxonomy()
    by_length = {
        "a": {11: 1.0},
        "b": {11: 1.0},
        "c": {12: 1.0},
        "uncertain": {12: 0.1},
    }
    resolved = genus_plurality(by_length, [], taxonomy, lengths={"a": 3, "b": 3, "c": 5})
    assert resolved["uncertain"] == {11: 1.0}
    tied = {"left": {12: 1.0}, "right": {11: 1.0}, "uncertain": {0: 0.4, 12: 0.3}}
    tied_resolved = genus_plurality(tied, [], taxonomy, calls={"uncertain": 10}, lengths={"left": 5, "right": 5})
    assert tied_resolved["uncertain"] == {11: 1.0}


def test_node_without_a_genus_stays_unchanged() -> None:
    """A taxon that does not walk to genus rank is copied."""
    taxonomy = _taxonomy()
    distributions = {"orphan": {999: 0.1}, "donor": {11: 1.0}}
    resolved = genus_plurality(distributions, [], taxonomy)
    assert resolved["orphan"] == {999: 0.1}


def test_thresholds_outside_unit_interval() -> None:
    """confident and uncertain outside [0, 1] are rejected. Endpoints are accepted."""
    taxonomy = _taxonomy()
    with pytest.raises(ValueError):
        genus_plurality({}, [], taxonomy, confident=1.1)
    with pytest.raises(ValueError):
        genus_plurality({}, [], taxonomy, uncertain=-0.01)
    assert genus_plurality({}, [], taxonomy, confident=0.0, uncertain=1.0) == {}

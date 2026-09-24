"""Mandatory checks for keeping a close same-genus pair."""

from __future__ import annotations

import pytest

from metamalevich.resolve import argmax_taxon
from metamalevich.margin_mixture import margin_mixture
from metamalevich.taxonomy import Taxonomy, Taxon

pytestmark = pytest.mark.mandatory


def _taxonomy() -> Taxonomy:
    """Two genera, each with two species, linked by parent ids."""
    return Taxonomy(
        [
            Taxon(10, None, "G", "GenusA", "toy", "v"),
            Taxon(11, 10, "S", "speciesA1", "toy", "v"),
            Taxon(12, 10, "S", "speciesA2", "toy", "v"),
            Taxon(20, None, "G", "GenusB", "toy", "v"),
            Taxon(21, 20, "S", "speciesB1", "toy", "v"),
        ],
        source="toy",
        version="v",
    )


def test_close_same_genus_pair_drops_the_third_taxon() -> None:
    """A 0.55/0.45 same-genus pair is kept and a third taxon is dropped."""
    taxonomy = _taxonomy()
    distributions = {"n": {11: 0.55, 12: 0.45, 21: 0.1}}
    original = {"n": {11: 0.55, 12: 0.45, 21: 0.1}}
    mixed = margin_mixture(distributions, [{"source": "n", "target": "m"}], taxonomy)
    assert set(mixed["n"]) == {11, 12}
    assert mixed["n"][11] == pytest.approx(0.55)
    assert mixed["n"][12] == pytest.approx(0.45)
    assert argmax_taxon(mixed["n"])[0] == 11
    assert distributions == original


def test_wide_same_genus_gap_is_copied_whole() -> None:
    """A 0.9/0.1 same-genus pair stays in the full distribution."""
    taxonomy = _taxonomy()
    distributions = {"n": {11: 0.9, 12: 0.1, 21: 0.05}}
    mixed = margin_mixture(distributions, [], taxonomy)
    assert mixed["n"] == {11: 0.9, 12: 0.1, 21: 0.05}
    assert mixed["n"] is not distributions["n"]


def test_different_genus_top_two_are_copied_whole() -> None:
    """Top labels from two genera are not reduced to a pair."""
    taxonomy = _taxonomy()
    distributions = {"n": {11: 0.55, 21: 0.45, 12: 0.1}}
    mixed = margin_mixture(distributions, [], taxonomy)
    assert mixed["n"] == {11: 0.55, 21: 0.45, 12: 0.1}


def test_argmax_stays_the_leader_and_inputs_are_not_mutated() -> None:
    """The mixture leader is the argmax, and the input maps are unchanged."""
    taxonomy = _taxonomy()
    distributions = {"n": {12: 0.5, 11: 0.5, 21: 0.2}}
    edges = [{"source": "n", "target": "m", "weight": 1.0}]
    original_weights = dict(distributions["n"])
    original_edges = [dict(edge) for edge in edges]
    mixed = margin_mixture(distributions, edges, taxonomy, margin=0.4)
    assert argmax_taxon(mixed["n"])[0] == argmax_taxon(original_weights)[0] == 11
    assert set(mixed["n"]) == {11, 12}
    assert distributions["n"] == original_weights
    assert edges == original_edges
    assert mixed["n"] is not distributions["n"]


def test_taxon_zero_is_not_a_same_genus_partner() -> None:
    """Unclassified mass blocks the pair even when the gap is small."""
    taxonomy = _taxonomy()
    distributions = {"n": {11: 0.55, 0: 0.45, 12: 0.1}}
    mixed = margin_mixture(distributions, [], taxonomy)
    assert mixed["n"] == {11: 0.55, 0: 0.45, 12: 0.1}


def test_margin_outside_unit_interval() -> None:
    """A margin outside [0, 1] is rejected. The endpoints are accepted."""
    taxonomy = _taxonomy()
    with pytest.raises(ValueError):
        margin_mixture({}, [], taxonomy, margin=1.1)
    with pytest.raises(ValueError):
        margin_mixture({}, [], taxonomy, margin=-0.01)
    assert margin_mixture({}, [], taxonomy, margin=0.0) == {}
    assert margin_mixture({}, [], taxonomy, margin=1.0) == {}

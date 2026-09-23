"""Mandatory contracts for taxonomy, colour aggregation, and the toy reprofile."""

from __future__ import annotations

import pytest

from metamalevich.aggregate import aggregate
from metamalevich.evidence import EvidenceGraph, colours_from_weights, make_layer
from metamalevich.taxonomy import parse_kraken_report
from metamalevich.toy import run_toy

pytestmark = pytest.mark.mandatory

REPORT = """\
100.00\t4\t0\tR\t1\troot
100.00\t4\t0\tD\t2\t  Bacteria
50.00\t2\t0\tG\t3\t    GenusA
25.00\t1\t1\tS\t10\t      speciesA
25.00\t1\t1\tS\t11\t      speciesB
"""


def test_lca_and_rank_rollup() -> None:
    """LCA stays at the genus when two species are both supported."""
    taxonomy = parse_kraken_report(REPORT, source="kraken2", version="unit")
    assert taxonomy.lca([10, 11]) == 3
    assert taxonomy.ancestor_at_rank(10, "S") == 10
    assert taxonomy.ancestor_at_rank(3, "S") is None
    assert taxonomy.parent(10) == 3


def test_majority_is_explicit_and_multilabel_survives() -> None:
    """A tie keeps both taxa, and replace is the only way to drop a layer."""
    taxonomy = parse_kraken_report(REPORT, source="kraken2", version="unit")
    tied = aggregate({10: 4, 11: 4}, "majority", taxonomy)
    assert tied == {10: 0.5, 11: 0.5}
    graph = EvidenceGraph()
    layer = make_layer(
        layer_id="first",
        source="kraken2",
        method="count",
        database="kraken2",
        database_version="unit",
        parameters={},
        taxonomy_version="unit",
        operation="replace",
    )
    graph.apply("node", {"n1": colours_from_weights({10: 2, 11: 2}, evidence_type="kmer", source="kraken2")}, layer)
    assert set(graph.distribution("node", "n1")) == {10, 11}
    with pytest.raises(ValueError):
        graph.apply(
            "node",
            {"n1": colours_from_weights({10: 1}, evidence_type="kmer", source="kraken2")},
            make_layer(
                layer_id="",
                source="kraken2",
                method="count",
                database="kraken2",
                database_version="unit",
                parameters={},
                taxonomy_version="unit",
                operation="replace",
            ),
        )


def test_toy_reprofile_beats_initial_colouring() -> None:
    """Graph resolution on the toy community lowers abundance L1."""
    result = run_toy()
    assert result["beats_initial"] is True
    assert result["resolved_l1"] < result["initial_l1"]
    assert result["accuracy"] == 1.0

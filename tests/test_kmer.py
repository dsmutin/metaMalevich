"""Mandatory check that the C++ composition graph is symmetrized."""

from __future__ import annotations

from pathlib import Path

import pytest

from metamalevich.native import kmer_graph
from metamalevich.tables import read_tsv

pytestmark = pytest.mark.mandatory


def test_kmer_graph_writes_both_directions(tmp_path: Path) -> None:
    """Identical compositions become a pair of directed edges."""
    fasta = tmp_path / "contigs.fasta"
    fasta.write_text(
        ">a\n" + ("ACGT" * 20) + "\n>b\n" + ("ACGT" * 20) + "\n>c\n" + ("TTTT" * 20) + "\n",
        encoding="utf-8",
    )
    edges_path = tmp_path / "edges.tsv"
    kmer_graph(fasta, edges_path, top_k=2, min_sim=0.15)
    pairs = {(row["source"], row["target"]) for row in read_tsv(edges_path)}
    assert ("a", "b") in pairs
    assert ("b", "a") in pairs
    lengths = {row["node_id"] for row in read_tsv(Path(str(edges_path) + ".lengths.tsv"))}
    assert lengths == {"a", "b", "c"}

"""Mandatory checks for the C++ composition graph."""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from metamalevich.native import compile_tool, kmer_graph
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


def test_kmer_graph_refuses_an_all_pairs_scan_past_the_limit(tmp_path: Path) -> None:
    """More sequences than --max-nodes stops before the quadratic comparison."""
    fasta = tmp_path / "contigs.fasta"
    fasta.write_text(">a\nACGT\n>b\nACGT\n>c\nTTTT\n", encoding="utf-8")
    binary = compile_tool("kmer_knn")
    result = subprocess.run(
        [str(binary), str(fasta), str(tmp_path / "edges.tsv"), "--max-nodes", "2"],
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode != 0
    assert "all-pairs" in result.stderr

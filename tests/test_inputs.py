"""Mandatory checks that disagreeing inputs fail before a profile is written."""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from metamalevich.bench import require_consistent_ids
from metamalevich.native import kraken_counts
from metamalevich.tables import read_tsv

pytestmark = pytest.mark.mandatory


def test_id_mismatch_names_the_dataset() -> None:
    """FASTA, Kraken, truth, and lengths must describe the same contigs."""
    with pytest.raises(ValueError, match="toy: contig ids"):
        require_consistent_ids(
            dataset="toy",
            sequences={"a"},
            calls={"a"},
            truth={"b"},
            lengths={"a"},
        )


def test_kraken_count_sums_kmer_pairs(tmp_path: Path) -> None:
    """Numeric taxid:count pairs are summed, and the classifier call is stored separately."""
    source = tmp_path / "kraken.output"
    source.write_text("C\tn1\t10\t4\t10:3 0:1\nU\tn2\t0\t4\t0:4\n", encoding="utf-8")
    counts = tmp_path / "counts.tsv"
    calls = tmp_path / "calls.tsv"
    kraken_counts(source, counts, calls)
    summed = {(row["seq_id"], row["taxon_id"]): int(row["count"]) for row in read_tsv(counts)}
    assert summed[("n1", "10")] == 3
    assert summed[("n1", "0")] == 1
    call_rows = {row["seq_id"]: row["taxon_id"] for row in read_tsv(calls)}
    assert call_rows == {"n1": "10", "n2": "0"}


def test_kraken_count_rejects_a_short_line(tmp_path: Path) -> None:
    """A truncated Kraken row is an error, not a skipped contig."""
    source = tmp_path / "kraken.output"
    source.write_text("C\tonly-three\t10\n", encoding="utf-8")
    with pytest.raises(subprocess.CalledProcessError):
        kraken_counts(source, tmp_path / "counts.tsv", tmp_path / "calls.tsv")


def test_kraken_count_rejects_duplicate_sequence_ids(tmp_path: Path) -> None:
    """The same sequence id must not be summed twice or keep only the last call."""
    line = "C\tn1\t10\t4\t10:3 0:1\n"
    source = tmp_path / "kraken.output"
    source.write_text(line + line, encoding="utf-8")
    with pytest.raises(subprocess.CalledProcessError):
        kraken_counts(source, tmp_path / "counts.tsv", tmp_path / "calls.tsv")

"""Build the C++ helpers used for Kraken aggregation and the k-mer graph."""

from __future__ import annotations

import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CPP = ROOT / "cpp"


def compile_tool(name: str) -> Path:
    """Compile ``cpp/<name>.cpp`` when the binary is missing or older than the source."""
    source = CPP / f"{name}.cpp"
    binary = CPP / name
    if not source.is_file():
        raise FileNotFoundError(source)
    if binary.is_file() and binary.stat().st_mtime >= source.stat().st_mtime:
        return binary
    subprocess.run(
        ["g++", "-O3", "-std=c++17", "-o", str(binary), str(source)],
        check=True,
    )
    return binary


def kraken_counts(kraken_output: Path, destination: Path, calls: Path | None = None) -> Path:
    """Write aggregated ``seq_id, taxon_id, count`` rows and optional classifier calls."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    binary = compile_tool("kraken_count")
    command = [str(binary), str(kraken_output), str(destination)]
    if calls is not None:
        calls.parent.mkdir(parents=True, exist_ok=True)
        command.append(str(calls))
    subprocess.run(command, check=True)
    return destination


def kmer_graph(fasta: Path, destination: Path, *, top_k: int, min_sim: float) -> Path:
    """Write a symmetrized minimizer-Jaccard graph and a sibling lengths file."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    binary = compile_tool("kmer_knn")
    subprocess.run(
        [
            str(binary),
            str(fasta),
            str(destination),
            "--top-k",
            str(top_k),
            "--min-sim",
            str(min_sim),
        ],
        check=True,
    )
    return destination

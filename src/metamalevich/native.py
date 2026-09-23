"""Build the C++ helpers used for Kraken aggregation and the k-mer graph."""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path


def _source_dirs() -> list[Path]:
    """Directories that may hold the C++ sources, first match wins."""
    found = []
    if os.environ.get("METAMALEVICH_CPP"):
        found.append(Path(os.environ["METAMALEVICH_CPP"]))
    found.append(Path(__file__).resolve().parent / "cpp")
    found.append(Path(__file__).resolve().parents[2] / "cpp")
    prefix = os.environ.get("CONDA_PREFIX")
    if prefix:
        found.append(Path(prefix) / "share" / "metamalevich" / "cpp")
    return found


def tool_source(name: str) -> Path:
    """Return the C++ source for a native tool."""
    filename = f"{name}.cpp"
    for directory in _source_dirs():
        path = directory / filename
        if path.is_file():
            return path
    raise FileNotFoundError(
        f"missing C++ source {filename}. Set METAMALEVICH_CPP or install the conda package."
    )


def _compiler() -> str:
    """Use ``g++`` from ``PATH``, including the conda compiler name."""
    for name in ("g++", "x86_64-conda-linux-gnu-g++"):
        found = shutil.which(name)
        if found:
            return found
    raise FileNotFoundError("g++ was not found. Install gxx_linux-64 from conda-forge.")


def _binary_dir(source: Path) -> Path:
    """Write the binary beside a writable checkout, otherwise into a cache directory."""
    repo_cpp = Path(__file__).resolve().parents[2] / "cpp"
    if source.parent.resolve() == repo_cpp.resolve() and os.access(repo_cpp, os.W_OK):
        return repo_cpp
    cache = Path(os.environ.get("METAMALEVICH_CACHE", Path.home() / ".cache" / "metamalevich"))
    cache.mkdir(parents=True, exist_ok=True)
    return cache


def compile_tool(name: str) -> Path:
    """Compile ``cpp/<name>.cpp`` when the binary is missing or older than the source."""
    source = tool_source(name)
    binary = _binary_dir(source) / name
    if binary.is_file() and binary.stat().st_mtime >= source.stat().st_mtime:
        return binary
    subprocess.run(
        [_compiler(), "-O3", "-std=c++17", "-o", str(binary), str(source)],
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
    """Write a symmetrized canonical 4-mer cosine graph and a sibling lengths file."""
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

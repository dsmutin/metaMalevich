"""Load a MetaMetro ToCUMG and keep the colour namespaces this tool asked for."""

from __future__ import annotations

from pathlib import Path
from typing import Sequence

from metametro.bench import load_bench_cdbg
from metametro.bench.colourings import colouring_namespace
from metametro.formats.cdbg.model import Cdbg


def load_selected(bench_dir: Path, colourings: Sequence[str] | None = None) -> Cdbg:
    """Return the CDBG from a ``benchbuild`` directory, filtered by colouring names.

    Each token is a MetaMetro colouring name (``kraken2``, ``decaying``,
    ``composition_kmeans``) or a colour-dictionary namespace (``taxon``,
    ``sample``). An empty list keeps every layer.
    """
    names = [str(name) for name in colourings or ()]
    if not names:
        return load_bench_cdbg(Path(bench_dir))
    namespaces: list[str] = []
    for name in names:
        mapped = colouring_namespace(name)
        token = mapped if mapped else name
        if token not in namespaces:
            namespaces.append(token)
    return load_bench_cdbg(Path(bench_dir), namespaces=namespaces)

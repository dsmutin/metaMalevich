"""Build a MetaMetro bench and run one resolver on its colour layers.

The resolver reads the colour namespaces MetaMetro wrote. It does not read
``ground_truth/`` and it does not score a simulation label.
"""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import Sequence

from metamalevich.resolve import gated_neighbour
from metamalevich.tocumg import load_selected

SOLVE_COLOURINGS = ("composition_kmeans", "decaying")
RESOLVER = "gated_neighbour"
_PHAGE_GENOMES = ("T1", "T3", "T4", "T5", "T7")
_PHAGE_PROGRAMS = ("samovar", "megahit", "megahit_toolkit")


def solve_benches(
    names: Sequence[str],
    workdir: Path,
    colourings: Sequence[str] | None = None,
    *,
    execute: bool = False,
) -> dict:
    """Build each bench name and resolve the ones that have a ToCUMG.

    ``colourings`` defaults to composition k-means and decaying neighbour
    signal. Those layers are MetaMetro colourings. An external bench such as
    ``phage_10`` is assembled only when ``execute`` is set, or when the phage
    programs and genome FASTA files are already on this machine. Otherwise the
    bench records its contract and ``solved`` is false.
    """
    from metametro.bench.build import build
    from metametro.bench.paths import repo_root
    from metametro.bench.registry import resolve

    requested = [str(name) for name in names]
    if not requested:
        raise ValueError("solve requires at least one bench name")
    layers = tuple(str(name) for name in (colourings or SOLVE_COLOURINGS))
    if not layers:
        raise ValueError("solve requires at least one colouring")
    root = Path(workdir)
    root.mkdir(parents=True, exist_ok=True)
    rows: list[dict] = []
    for name in requested:
        spec = resolve(name)
        destination = root / spec.name
        run_execute = bool(execute) or _phage_ready(spec.name, repo_root())
        result = build(
            name,
            outdir=destination,
            colourings=layers,
            execute=run_execute,
        )
        row = {
            "bench": name,
            "canonical_name": result.spec.name,
            "build_status": result.status,
            "solved": False,
            "colourings": list(result.colourings),
        }
        if result.status == "built":
            row.update(_resolve_built(result.outdir, layers))
            applied = set(result.colourings)
            skipped = [name for name in layers if name not in applied]
            if skipped:
                row["colourings_skipped"] = skipped
        elif result.status == "contract":
            row["missing"] = _phage_missing(result.spec.name, repo_root())
        else:
            raise ValueError(f"{name} finished with status {result.status}")
        rows.append(row)
    solved = [row for row in rows if row["solved"]]
    unresolved_graph = [
        row["bench"] for row in rows if row["build_status"] == "built" and not row["solved"]
    ]
    ok = bool(solved) and not unresolved_graph
    return {
        "status": "ok" if ok else "incomplete",
        "ok": ok,
        "resolver": RESOLVER,
        "colourings": list(layers),
        "benches": rows,
    }


def _resolve_built(outdir: Path, colourings: Sequence[str]) -> dict:
    """Run gated neighbour smoothing on one built ToCUMG."""
    graph = load_selected(outdir, colourings)
    namespaces = sorted({str(row["namespace"]) for row in graph.colors or []})
    distributions = _colour_distributions(graph)
    coloured = sum(1 for weights in distributions.values() if weights)
    if coloured == 0:
        raise ValueError(
            f"{outdir} has no unitig colours in namespaces {', '.join(namespaces) or '(none)'}"
        )
    edges = [
        {
            "edge_id": link.link_id,
            "source": link.source,
            "target": link.target,
            "weight": 1.0,
        }
        for link in graph.links
    ]
    gated_neighbour(distributions, edges)
    return {
        "solved": True,
        "n_unitigs": len(graph.unitigs),
        "n_links": len(graph.links),
        "n_coloured_unitigs": coloured,
        "namespaces": namespaces,
    }


def _colour_distributions(graph) -> dict[str, dict[int, float]]:
    """Uniform mass on the colour ids MetaMetro left on each unitig."""
    known: set[int] = set()
    for row in graph.colors or []:
        token = str(row.get("color_id", "")).strip()
        if token.lstrip("-").isdigit():
            known.add(int(token))
    distributions: dict[str, dict[int, float]] = {}
    for unitig in graph.unitigs:
        colour_ids = []
        for colour_id in unitig.color_ids:
            value = int(colour_id)
            if value in known and value not in colour_ids:
                colour_ids.append(value)
        if not colour_ids:
            distributions[unitig.unitig_id] = {}
            continue
        share = 1.0 / len(colour_ids)
        distributions[unitig.unitig_id] = {colour_id: share for colour_id in colour_ids}
    return distributions


def _phage_ready(canonical_name: str, root: Path) -> bool:
    """True when this machine can assemble the five-phage bench."""
    return canonical_name.startswith("phage_species_5") and not _phage_missing(canonical_name, root)


def _phage_missing(canonical_name: str, root: Path) -> list[str]:
    """Programs and FASTA files the phage bench still needs. Empty for other benches."""
    if not canonical_name.startswith("phage_species_5"):
        return []
    missing = [program for program in _PHAGE_PROGRAMS if shutil.which(program) is None]
    genomes = root / "data" / "raw" / "genomes"
    for name in _PHAGE_GENOMES:
        path = genomes / f"{name}.fna"
        if not path.is_file() or path.stat().st_size == 0:
            missing.append(f"data/raw/genomes/{name}.fna")
    return missing

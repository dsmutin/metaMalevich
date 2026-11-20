"""Build a MetaMetro bench and run one resolver on its taxonomic colour layer.

The resolver reads taxa only from the ``kraken2`` namespace, whose colour
values are NCBI taxon ids, checked against an NCBI TAXDUMP. It does not read
``ground_truth/`` and it does not score a simulation label.
"""

from __future__ import annotations

import os
import shutil
from pathlib import Path
from typing import Sequence

from metamalevich.resolve import gated_neighbour
from metamalevich.taxonomy import Taxonomy, parse_taxdump
from metamalevich.tocumg import load_selected

TAXONOMIC_COLOURING = "kraken2"
TAXONOMIC_NAMESPACE = "kraken2"
SOLVE_COLOURINGS = (TAXONOMIC_COLOURING, "composition_kmeans", "decaying")
RESOLVER = "gated_neighbour"
NO_KRAKEN_REASON = (
    "no kraken2 input: pass kraken_output or kraken_db (--kraken-output, --kraken-db). "
    "Pinned kraken2 outputs for the MetaMetro benches are not released yet (MetaMetro #1163)"
)
_PHAGE_GENOMES = ("T1", "T3", "T4", "T5", "T7")
_PHAGE_PROGRAMS = ("samovar", "megahit", "megahit_toolkit")


def solve_benches(
    names: Sequence[str],
    workdir: Path,
    colourings: Sequence[str] | None = None,
    *,
    execute: bool = False,
    kraken_db: Path | None = None,
    kraken_output: Path | None = None,
    kaiju_db: Path | None = None,
    kaiju_output: Path | None = None,
    taxdump: Path | None = None,
) -> dict:
    """Build each bench name and resolve the ones that have a ToCUMG.

    ``colourings`` defaults to ``kraken2``, composition k-means, and the
    decaying neighbour signal. ``kraken2`` is the taxonomic layer and is
    required. The classifier paths go to MetaMetro ``build()`` unchanged, so
    one ``kraken_output`` serves every bench in the call. ``taxdump`` is an
    NCBI TAXDUMP directory (``nodes.dmp``, ``names.dmp``); it defaults to the
    ``TAXDUMP`` environment variable and is required once a bench has a
    kraken2 input.

    Without ``kraken_output`` or ``kraken_db`` an in-process bench, or an
    external bench that would be assembled, is skipped with a reason. An
    external bench such as ``phage_10`` is assembled only when ``execute`` is
    set, or when the phage programs and genome FASTA files are already on this
    machine. Otherwise the bench records its contract and ``solved`` is false.
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
    if TAXONOMIC_COLOURING not in layers:
        raise ValueError(
            f"solve requires the taxonomic colouring {TAXONOMIC_COLOURING}; "
            f"got {', '.join(layers)}"
        )
    has_kraken = kraken_output is not None or kraken_db is not None
    taxonomy: Taxonomy | None = None
    if has_kraken:
        taxonomy = _load_taxdump(taxdump)
    classifiers = {
        "kraken_db": kraken_db,
        "kraken_output": kraken_output,
        "kaiju_db": kaiju_db,
        "kaiju_output": kaiju_output,
    }
    root = Path(workdir)
    root.mkdir(parents=True, exist_ok=True)
    rows: list[dict] = []
    for name in requested:
        spec = resolve(name)
        destination = root / spec.name
        run_execute = bool(execute) or _phage_ready(spec.name, repo_root())
        if not has_kraken and (spec.kind == "inprocess" or run_execute):
            row = {
                "bench": name,
                "canonical_name": spec.name,
                "build_status": "skipped",
                "solved": False,
                "skipped": NO_KRAKEN_REASON,
            }
            missing = _phage_missing(spec.name, repo_root())
            if missing:
                row["missing"] = missing
            rows.append(row)
            continue
        result = build(
            name,
            outdir=destination,
            colourings=layers,
            execute=run_execute,
            **classifiers,
        )
        row = {
            "bench": name,
            "canonical_name": result.spec.name,
            "build_status": result.status,
            "solved": False,
            "colourings": list(result.colourings),
        }
        if result.status == "built":
            if taxonomy is None:
                raise ValueError(f"{name} was built without a kraken2 input")
            row.update(_resolve_built(result.outdir, layers, taxonomy))
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
    skipped_rows = [row for row in rows if row["build_status"] == "skipped"]
    unresolved_graph = [
        row["bench"] for row in rows if row["build_status"] == "built" and not row["solved"]
    ]
    ok = not unresolved_graph and bool(solved or skipped_rows)
    if not ok:
        status = "incomplete"
    elif solved:
        status = "ok"
    else:
        status = "skipped"
    summary = {
        "status": status,
        "ok": ok,
        "resolver": RESOLVER,
        "taxonomic_layer": TAXONOMIC_COLOURING,
        "colourings": list(layers),
        "benches": rows,
    }
    if taxonomy is not None:
        summary["taxdump"] = taxonomy.version
    return summary


def _load_taxdump(taxdump: Path | None) -> Taxonomy:
    """NCBI TAXDUMP from the argument or ``TAXDUMP``. Stop if neither is set."""
    path = taxdump if taxdump is not None else os.environ.get("TAXDUMP")
    if not path:
        raise ValueError(
            "solve needs an NCBI TAXDUMP directory for the kraken2 layer: "
            "pass taxdump (--taxdump) or set TAXDUMP"
        )
    return parse_taxdump(Path(path))


def _resolve_built(outdir: Path, colourings: Sequence[str], taxonomy: Taxonomy) -> dict:
    """Run gated neighbour smoothing on the kraken2 taxa of one built ToCUMG."""
    graph = load_selected(outdir, colourings)
    namespaces = sorted({str(row["namespace"]) for row in graph.colors or []})
    if TAXONOMIC_NAMESPACE not in namespaces:
        raise ValueError(
            f"{outdir} has no {TAXONOMIC_NAMESPACE} colours; namespaces: {', '.join(namespaces) or '(none)'}"
        )
    distributions = _colour_distributions(graph)
    taxa = sorted({taxon for weights in distributions.values() for taxon in weights})
    unknown = [taxon for taxon in taxa if taxon not in taxonomy.by_id]
    if unknown:
        shown = ", ".join(str(taxon) for taxon in unknown[:10])
        raise ValueError(
            f"{outdir}: {len(unknown)} kraken2 taxon ids are not in TAXDUMP {taxonomy.version}: {shown}"
        )
    coloured = sum(1 for weights in distributions.values() if weights)
    if coloured == 0:
        raise ValueError(f"{outdir}: the {TAXONOMIC_NAMESPACE} layer colours no unitig")
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
        "n_taxa": len(taxa),
        "namespaces": namespaces,
    }


def _colour_distributions(graph) -> dict[str, dict[int, float]]:
    """Uniform mass on the kraken2 taxon ids MetaMetro left on each unitig.

    Only the ``kraken2`` namespace holds taxa; its colour value is the NCBI
    taxon id. Cluster ids from other namespaces are not taxa.
    """
    taxon_of: dict[int, int] = {}
    for row in graph.colors or []:
        if str(row.get("namespace")) != TAXONOMIC_NAMESPACE:
            continue
        colour = str(row.get("color_id", "")).strip()
        value = str(row.get("value", "")).strip()
        if not value.isdigit():
            raise ValueError(f"{TAXONOMIC_NAMESPACE} colour {colour} has a non-taxid value {value!r}")
        taxon_of[int(colour)] = int(value)
    distributions: dict[str, dict[int, float]] = {}
    for unitig in graph.unitigs:
        taxa: list[int] = []
        for colour_id in unitig.color_ids:
            taxon = taxon_of.get(int(colour_id))
            if taxon is not None and taxon != 0 and taxon not in taxa:
                taxa.append(taxon)
        if not taxa:
            distributions[unitig.unitig_id] = {}
            continue
        share = 1.0 / len(taxa)
        distributions[unitig.unitig_id] = {taxon: share for taxon in taxa}
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

"""Build a MetaMetro bench and run one resolver on its taxonomic colour layer.

The resolver reads taxa only from the ``kraken2`` namespace, whose colour
values are NCBI taxon ids, checked against an NCBI TAXDUMP. Each solved bench
gets per-unitig assignments and a length-weighted profile. Gold labels are
read from ``ground_truth/`` with MetaMetro ``load_bench_gold`` only after the
resolver has run, and only to score it.
"""

from __future__ import annotations

import csv
import os
from pathlib import Path
from typing import Sequence

from metamalevich.evaluate import classification_scores
from metamalevich.genus_plurality import genus_plurality
from metamalevich.resolve import argmax_taxon, resolve
from metamalevich.taxonomy import Taxonomy, parse_taxdump
from metamalevich.tocumg import load_selected

TAXONOMIC_COLOURING = "kraken2"
TAXONOMIC_NAMESPACE = "kraken2"
SOLVE_COLOURINGS = (TAXONOMIC_COLOURING, "composition_kmeans", "decaying")
RESOLVERS = ("genus_plurality", "gated_neighbour", "probability_sum", "hard")
RESOLVER = RESOLVERS[0]
GOLD_NAMESPACE = "taxon"
NO_KRAKEN_REASON = (
    "no kraken2 input: no kraken_output or kraken_db (--kraken-output, --kraken-db) was given, "
    "and MetaMetro pins no Kraken2 output for this bench"
)
_PHAGE_PREFIX = "phage_species_5"
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
    resolver: str = RESOLVER,
) -> dict:
    """Build each bench name and resolve the ones that have a ToCUMG.

    ``colourings`` defaults to ``kraken2``, composition k-means, and the
    decaying neighbour signal. ``kraken2`` is the taxonomic layer and is
    required. The classifier paths go to MetaMetro ``build()`` unchanged, so
    one ``kraken_output`` serves every bench in the call. ``taxdump`` is an
    NCBI TAXDUMP directory (``nodes.dmp``, ``names.dmp``); it defaults to the
    ``TAXDUMP`` environment variable. Without either, a bench whose kraken2
    layer comes from a MetaMetro pin uses the trimmed TAXDUMP pinned beside it
    (``metametro.bench.kraken_pin.pinned_taxdump``); any other bench with a
    kraken2 input stops before it is built.

    A bench has a kraken2 input when ``kraken_output`` or ``kraken_db`` is
    given, or when MetaMetro pins a Kraken2 output for it
    (``metametro.bench.kraken_pin.kraken_pin``, as for ``phage_10``); MetaMetro
    then paints the layer from the pin itself. Without either, an in-process
    bench, or an external bench that would be assembled, is skipped with a
    reason. An
    external bench such as ``phage_10`` is assembled when ``execute`` is set,
    or when MetaMetro's own checks find its pinned genomes and programs.
    Otherwise the bench records its contract, ``missing`` lists the errors
    MetaMetro raised (with the ``metametro benchfetch`` command for missing
    genomes), and ``solved`` is false.

    ``resolver`` is one of ``RESOLVERS`` (default ``genus_plurality``). A
    solved bench writes ``assignments.tsv`` and ``profile.tsv`` under
    ``<workdir>/solve/<bench>/``. ``metrics`` scores the resolved argmax
    against the bench gold ``taxon`` labels when those labels are NCBI taxon
    ids; otherwise ``metrics`` is null and ``metrics_reason`` says why.
    """
    from metametro.bench.build import build
    from metametro.bench.registry import resolve

    if resolver not in RESOLVERS:
        raise ValueError(f"unknown resolver {resolver}; choose one of {', '.join(RESOLVERS)}")
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
    chosen = taxdump if taxdump is not None else os.environ.get("TAXDUMP") or None
    taxonomy: Taxonomy | None = None
    if has_kraken or chosen is not None:
        taxonomy = _load_taxdump(chosen)
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
        missing = _phage_preflight(spec)
        run_execute = bool(execute) or (spec.name.startswith(_PHAGE_PREFIX) and not missing)
        pinned = not has_kraken and _has_kraken_pin(spec.name)
        will_build = spec.kind == "inprocess" or run_execute
        if not (has_kraken or pinned) and will_build:
            row = {
                "bench": name,
                "canonical_name": spec.name,
                "build_status": "skipped",
                "solved": False,
                "skipped": NO_KRAKEN_REASON,
            }
            if missing:
                row["missing"] = missing
            rows.append(row)
            continue
        bench_taxonomy = taxonomy
        if will_build and bench_taxonomy is None:
            bench_taxonomy = _pinned_taxonomy(spec.name)
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
        if pinned:
            row["kraken2_input"] = "metametro_pin"
        if result.status == "built":
            if bench_taxonomy is None:
                raise ValueError(f"{name} was built without a TAXDUMP")
            row["taxdump"] = bench_taxonomy.version
            row.update(
                _resolve_built(
                    result.outdir, layers, bench_taxonomy, resolver, root / "solve" / spec.name
                )
            )
            applied = set(result.colourings)
            skipped = [name for name in layers if name not in applied]
            if skipped:
                row["colourings_skipped"] = skipped
        elif result.status == "contract":
            row["missing"] = missing
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
        "resolver": resolver,
        "taxonomic_layer": TAXONOMIC_COLOURING,
        "colourings": list(layers),
        "benches": rows,
    }
    if taxonomy is not None:
        summary["taxdump"] = taxonomy.version
    return summary


def _load_taxdump(path: Path | str | None) -> Taxonomy:
    """NCBI TAXDUMP from ``--taxdump`` or ``TAXDUMP``. Stop if neither is set."""
    if not path:
        raise ValueError(
            "solve needs an NCBI TAXDUMP directory for the kraken2 layer: "
            "pass taxdump (--taxdump) or set TAXDUMP"
        )
    return parse_taxdump(Path(path))


def _pinned_taxonomy(bench: str) -> Taxonomy:
    """The trimmed TAXDUMP MetaMetro pins beside ``bench``'s Kraken2 output.

    Used only when no TAXDUMP was given and the bench's kraken2 layer comes
    from that pin. The version is ``<Kraken2 index name>-<bench>``.
    """
    from metametro.bench.kraken_pin import kraken_pin, pinned_taxdump

    found = kraken_pin(bench)
    directory = pinned_taxdump(bench) if found is not None else None
    if directory is None:
        raise ValueError(
            f"{bench} has no TAXDUMP: pass taxdump (--taxdump) or set TAXDUMP; "
            "MetaMetro pins no TAXDUMP for this bench"
        )
    index = (found[1].get("index") or {}).get("name") or "kraken2_pin"
    return parse_taxdump(Path(directory), version=f"{index}-{found[1].get('bench', bench)}")


def _has_kraken_pin(bench: str) -> bool:
    """True when MetaMetro pins a Kraken2 output for ``bench`` (checked against its SHA-256)."""
    from metametro.bench.kraken_pin import kraken_pin

    return kraken_pin(bench) is not None


def _resolve_built(
    outdir: Path, colourings: Sequence[str], taxonomy: Taxonomy, resolver: str, destination: Path
) -> dict:
    """Run ``resolver`` on the kraken2 taxa of one built ToCUMG and write its outputs."""
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
    lengths = {unitig.unitig_id: len(unitig.sequence) for unitig in graph.unitigs}
    if resolver == "genus_plurality":
        resolved = genus_plurality(distributions, edges, taxonomy, lengths=lengths)
    else:
        resolved = resolve(distributions, edges, resolver)
    profile = _length_profile(resolved, lengths)
    destination.mkdir(parents=True, exist_ok=True)
    assignments_path = destination / "assignments.tsv"
    profile_path = destination / "profile.tsv"
    _write_assignments(assignments_path, resolved, lengths, taxonomy)
    _write_profile(profile_path, profile, taxonomy)
    metrics, reason = _gold_metrics(outdir, graph, resolved)
    row = {
        "solved": True,
        "n_unitigs": len(graph.unitigs),
        "n_links": len(graph.links),
        "n_coloured_unitigs": coloured,
        "n_taxa": len(taxa),
        "namespaces": namespaces,
        "assignments": str(assignments_path),
        "profile": str(profile_path),
        "profile_fractions": {str(taxon): fraction for taxon, fraction in profile.items()},
        "metrics": metrics,
    }
    if metrics is None:
        row["metrics_reason"] = reason
    return row


def _length_profile(
    resolved: dict[str, dict[int, float]], lengths: dict[str, int]
) -> dict[int, float]:
    """Length-weighted sum of resolved distributions. Unresolved length is taxon 0."""
    mass: dict[int, float] = {}
    for unitig_id, length in lengths.items():
        weights = resolved.get(unitig_id) or {}
        total = sum(value for value in weights.values() if value > 0)
        if total <= 0:
            mass[0] = mass.get(0, 0.0) + length
            continue
        for taxon, value in weights.items():
            if value > 0:
                mass[taxon] = mass.get(taxon, 0.0) + length * value / total
    grand = sum(mass.values())
    if grand <= 0:
        return {}
    return {taxon: mass[taxon] / grand for taxon in sorted(mass, key=lambda key: (-mass[key], key))}


def _write_assignments(
    path: Path, resolved: dict[str, dict[int, float]], lengths: dict[str, int], taxonomy: Taxonomy
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["unitig_id", "length", "taxon_id", "name", "rank", "probability", "distribution"])
        for unitig_id, length in lengths.items():
            weights = resolved.get(unitig_id) or {}
            taxon, probability = argmax_taxon(weights)
            distribution = ",".join(
                f"{key}:{weights[key]:.6g}" for key in sorted(weights, key=lambda key: (-weights[key], key))
            )
            writer.writerow(
                [unitig_id, length, taxon, taxonomy.name(taxon), taxonomy.rank(taxon), f"{probability:.6g}", distribution]
            )


def _write_profile(path: Path, profile: dict[int, float], taxonomy: Taxonomy) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["taxon_id", "name", "rank", "fraction"])
        for taxon, fraction in profile.items():
            writer.writerow([taxon, taxonomy.name(taxon), taxonomy.rank(taxon), f"{fraction:.6g}"])


def _gold_metrics(
    outdir: Path, graph, resolved: dict[str, dict[int, float]]
) -> tuple[dict | None, str | None]:
    """Score resolved argmax taxa against gold taxon ids, or say why not.

    Gold is MetaMetro ``ground_truth/taxon_nodes.tsv``. A unitig's gold is the
    union over its CFA members. Gold values must be NCBI taxon ids; bench
    label names such as ``taxon_1`` are not mapped to taxa.
    """
    from metametro.bench import load_bench_gold
    from metametro.errors import ContractError

    try:
        gold = load_bench_gold(outdir, GOLD_NAMESPACE)
    except ContractError as exc:
        return None, f"no gold {GOLD_NAMESPACE} labels on this bench: {exc}"
    names = [value for value in gold.values if not value.isdigit()]
    if names:
        shown = ", ".join(names[:5])
        return None, (
            f"gold {GOLD_NAMESPACE} labels are bench names ({shown}), not NCBI taxon ids, "
            "and ground_truth has no mapping from them to taxon ids"
        )
    by_unitig: dict[str, set[int]] = {}
    for unitig in graph.unitigs:
        labels = {int(value) for member in unitig.members for value in gold.nodes.get(member, ())}
        if labels:
            by_unitig[unitig.unitig_id] = labels
    if not by_unitig:
        return None, f"no unitig carries a gold {GOLD_NAMESPACE} label"
    single = {unitig_id: next(iter(labels)) for unitig_id, labels in by_unitig.items() if len(labels) == 1}
    hits = sum(
        1 for unitig_id, labels in by_unitig.items() if argmax_taxon(resolved.get(unitig_id) or {})[0] in labels
    )
    scores = classification_scores(resolved, single) if single else {}
    return {
        "gold_namespace": GOLD_NAMESPACE,
        "n_labelled_unitigs": len(by_unitig),
        "n_single_label_unitigs": len(single),
        "accuracy_any": hits / len(by_unitig),
        "accuracy": scores.get("accuracy"),
        "precision": scores.get("precision"),
        "recall": scores.get("recall"),
        "f1": scores.get("f1"),
    }, None


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


def _phage_preflight(spec) -> list[str]:
    """Errors MetaMetro raises before it would assemble a phage bench. Empty when ready.

    These are the checks ``build(..., execute=True)`` runs before writing:
    the pinned genome FASTA files (``check_genomes``) and the programs on PATH.
    Other benches return an empty list.
    """
    if not spec.name.startswith(_PHAGE_PREFIX):
        return []
    from metametro.bench.data.universal.catalog import require_programs
    from metametro.bench.genomes import check_genomes
    from metametro.bench.paths import repo_root
    from metametro.errors import ContractError

    errors: list[str] = []
    try:
        check_genomes(spec.name, repo_root() / "data" / "raw" / "genomes")
    except ContractError as exc:
        errors.extend(exc.errors)
    try:
        require_programs(spec, _PHAGE_PROGRAMS)
    except ContractError as exc:
        errors.extend(exc.errors)
    return errors

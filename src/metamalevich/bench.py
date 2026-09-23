"""Run colouring and resolution hypotheses on a closed community.

The k-mer graph is rebuilt from contig sequences with the ``top_k`` and
``min_sim`` recorded in ``graph/graph_meta.txt``. The original edge list was
not in the bundle, so edge ids here are new and are not joined to
``strategy_b`` rows.
"""

from __future__ import annotations

from pathlib import Path

import yaml

from metamalevich import __version__
from metamalevich.aggregate import aggregate, normalize, roll_counts
from metamalevich.bridge import export_tocumg
from metamalevich.evaluate import abundance_scores, classification_scores, path_consistency
from metamalevich.evidence import colours_from_weights, make_layer
from metamalevich.evidence import EvidenceGraph
from metamalevich.native import kmer_graph, kraken_counts, tool_source
from metamalevich.plots import write_abundance_chart, write_summary_charts
from metamalevich.reprofile import profile_nodes
from metamalevich.resolve import edge_distributions, hard_assignment, resolve
from metamalevich.tables import read_csv, read_fasta, read_tsv, write_tsv
from metamalevich.taxonomy import Taxonomy, parse_kraken_report

DATASETS = {
    "samovar10": {
        "fasta": "samovar10/assembly/assembly.fasta",
        "kraken": "samovar10/kraken2/whole_contig.output",
        "report": "samovar10/kraken2/whole_contig.report",
        "truth": "samovar10/ground_truth/ground_truth.csv",
        "abundance": "samovar10/ground_truth/genome_to_taxonomy.csv",
        "top_k": 8,
        "min_sim": 0.15,
        "graph_note": "samovar10/graph/graph_meta.txt records top_k=8 and min_sim=0.15",
    },
    "samovar10_ont1b": {
        "fasta": "samovar10_ont1b/assembly/assembly.fasta",
        "kraken": "samovar10_ont1b/kraken2/whole_contig.output",
        "report": "samovar10_ont1b/kraken2/whole_contig.report",
        "truth": "samovar10_ont1b/ground_truth/ground_truth.csv",
        "abundance": "samovar10_ont1b/ground_truth/genome_to_taxonomy.csv",
        "top_k": 8,
        "min_sim": 0.15,
        "graph_note": "samovar10_ont1b/graph/graph_meta.txt records top_k=8 and min_sim=0.15",
    },
}

INPUT_GRAPH = "canonical 4-mer cosine kNN rebuilt from assembly contigs"
GRAPH_CONSTRUCTION_METHOD = "symmetrized top-k cosine on canonical 4-mer frequencies"

ABUNDANCE_CHART_HYPOTHESES = (
    "initial_colouring",
    "probability_sum",
    "gated_neighbour",
    "bayesian_edge",
)

HYPOTHESES = (
    "initial_colouring",
    "probability_sum",
    "lca",
    "gated_neighbour",
    "bayesian_edge",
)

PROFILE_COLUMNS = [
    "taxon_id",
    "rank",
    "name",
    "estimated_abundance",
    "relative_abundance",
    "assigned_bases",
    "assigned_reads",
    "unique_bases",
    "shared_bases",
    "ambiguous_bases",
    "confidence",
]


def cache_signature(payload: dict) -> str:
    """Stable sidecar text for one cached table."""
    return yaml.safe_dump(payload, sort_keys=True)


def cache_is_fresh(path: Path, signature: str) -> bool:
    """True when ``path`` exists and its sidecar still matches ``signature``."""
    meta = Path(str(path) + ".meta")
    if not path.is_file() or not meta.is_file():
        return False
    return meta.read_text(encoding="utf-8") == signature


def write_cache_signature(path: Path, signature: str) -> None:
    """Record the inputs that produced ``path``."""
    Path(str(path) + ".meta").write_text(signature, encoding="utf-8")


def _group_counts(rows: list[dict[str, str]]) -> dict[str, dict[int, float]]:
    grouped: dict[str, dict[int, float]] = {}
    for row in rows:
        seq_id = row["seq_id"]
        taxon_id = int(row["taxon_id"])
        table = grouped.setdefault(seq_id, {})
        table[taxon_id] = table.get(taxon_id, 0.0) + float(row["count"])
    return grouped


def _species_truth(taxonomy: Taxonomy, taxon_id: int) -> int:
    rolled = taxonomy.ancestor_at_rank(taxon_id, "S")
    return taxon_id if rolled is None else rolled


def _one_hot_call(call: int, taxonomy: Taxonomy) -> dict[int, float]:
    if call == 0:
        return {0: 1.0}
    species = taxonomy.ancestor_at_rank(call, "S")
    return {0 if species is None else species: 1.0}


def _lca_distribution(counts: dict[int, float], taxonomy: Taxonomy) -> dict[int, float]:
    weights = aggregate(counts, "lca", taxonomy, min_fraction=0.05)
    if not weights:
        return {0: 1.0}
    taxon_id = next(iter(weights))
    if taxon_id == 0:
        return {0: 1.0}
    species = taxonomy.ancestor_at_rank(taxon_id, "S")
    return {0 if species is None else species: 1.0}


def _evidence_distributions(
    counts: dict[str, dict[int, float]],
    taxonomy: Taxonomy,
    calls: dict[str, int],
) -> dict[str, dict[int, float]]:
    """Species-level k-mer distributions.

    K-mers placed above species, including unclassified, are not spread across
    child species. A node with no species-level k-mers keeps the Kraken call
    when that call itself rolls to species.
    """
    distributions = {}
    for node_id in set(counts) | set(calls):
        rolled = roll_counts(counts.get(node_id, {}), taxonomy, "S")
        rolled.pop(0, None)
        if rolled:
            distributions[node_id] = normalize(rolled)
        else:
            distributions[node_id] = _one_hot_call(calls.get(node_id, 0), taxonomy)
    return distributions


def require_consistent_ids(
    *,
    dataset: str,
    sequences: set[str],
    calls: set[str],
    truth: set[str],
    lengths: set[str],
) -> None:
    """Fail when contig ids disagree across FASTA, Kraken calls, truth, and lengths."""
    problems = []
    if sequences != calls:
        problems.append(f"fasta/kraken symmetric difference {sorted(sequences ^ calls)[:8]}")
    if sequences != truth:
        problems.append(f"fasta/truth symmetric difference {sorted(sequences ^ truth)[:8]}")
    if sequences != lengths:
        problems.append(f"fasta/lengths symmetric difference {sorted(sequences ^ lengths)[:8]}")
    if not sequences:
        problems.append("FASTA has no sequences")
    if problems:
        raise ValueError(f"{dataset}: contig ids do not match: " + "; ".join(problems))


def _load_edges(path: Path) -> list[dict]:
    edges = []
    for row in read_tsv(path):
        edges.append(
            {
                "edge_id": row["edge_id"],
                "source": row["source"],
                "target": row["target"],
                "orientation": row.get("orientation") or "++",
                "weight": float(row["weight"]),
            }
        )
    return edges


def _relative_from_profile(rows: list[dict], unclassified_fraction: float) -> dict[int, float]:
    predicted = {int(row["taxon_id"]): float(row["relative_abundance"]) for row in rows}
    if unclassified_fraction > 0:
        predicted[0] = predicted.get(0, 0.0) + unclassified_fraction
    return predicted


def run_dataset(data_root: Path, name: str, *, intermediate: Path, benchmark: Path) -> list[dict]:
    """Colour one dataset, resolve it, and write one folder per hypothesis."""
    spec = DATASETS[name]
    fasta = data_root / spec["fasta"]
    for label, path in (
        ("assembly", fasta),
        ("kraken", data_root / spec["kraken"]),
        ("report", data_root / spec["report"]),
        ("truth", data_root / spec["truth"]),
        ("abundance", data_root / spec["abundance"]),
    ):
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(f"{name}: missing or empty {label} file: {path}")
    kraken = data_root / spec["kraken"]
    report_path = data_root / spec["report"]
    work = intermediate / name
    work.mkdir(parents=True, exist_ok=True)
    counts_path = work / "kraken_counts.tsv"
    calls_path = work / "kraken_calls.tsv"
    edges_path = work / "knn_edges.tsv"
    kraken_signature = cache_signature(
        {
            "tool": "kraken_count",
            "source_mtime_ns": tool_source("kraken_count").stat().st_mtime_ns,
            "input_size": kraken.stat().st_size,
            "input_mtime_ns": kraken.stat().st_mtime_ns,
        }
    )
    if not (cache_is_fresh(counts_path, kraken_signature) and calls_path.is_file()):
        kraken_counts(kraken, counts_path, calls_path)
        write_cache_signature(counts_path, kraken_signature)
    graph_signature = cache_signature(
        {
            "tool": "kmer_knn",
            "source_mtime_ns": tool_source("kmer_knn").stat().st_mtime_ns,
            "fasta_size": fasta.stat().st_size,
            "fasta_mtime_ns": fasta.stat().st_mtime_ns,
            "top_k": int(spec["top_k"]),
            "min_sim": float(spec["min_sim"]),
        }
    )
    lengths_path = Path(str(edges_path) + ".lengths.tsv")
    if not (cache_is_fresh(edges_path, graph_signature) and lengths_path.is_file()):
        kmer_graph(fasta, edges_path, top_k=int(spec["top_k"]), min_sim=float(spec["min_sim"]))
        write_cache_signature(edges_path, graph_signature)
    taxonomy = parse_kraken_report(
        report_path.read_text(encoding="utf-8"),
        source="kraken2",
        version=spec["report"],
    )
    counts = _group_counts(read_tsv(counts_path))
    calls = {row["seq_id"]: int(row["taxon_id"]) for row in read_tsv(calls_path)}
    edges = _load_edges(edges_path)
    lengths = {row["node_id"]: int(row["length"]) for row in read_tsv(Path(str(edges_path) + ".lengths.tsv"))}
    truth_rows = read_csv(data_root / spec["truth"])
    node_truth = {}
    for row in truth_rows:
        node_truth[row["contig_id"]] = _species_truth(taxonomy, int(row["taxon_id"]))
        if row.get("length"):
            lengths[row["contig_id"]] = int(float(row["length"]))
    sequence_ids = {node_id for node_id, _sequence in read_fasta(fasta)}
    require_consistent_ids(
        dataset=name,
        sequences=sequence_ids,
        calls=set(calls),
        truth=set(node_truth),
        lengths=set(lengths),
    )
    if any(lengths[node_id] <= 0 for node_id in sequence_ids):
        raise ValueError(f"{name}: every contig needs a positive length")
    abundance_rows = read_csv(data_root / spec["abundance"])
    truth_abundance: dict[int, float] = {}
    for row in abundance_rows:
        taxon_id = _species_truth(taxonomy, int(row["taxon_id"]))
        truth_abundance[taxon_id] = truth_abundance.get(taxon_id, 0.0) + float(row["genome_abundance"])

    evidence = _evidence_distributions(counts, taxonomy, calls)
    calls_dist = {node_id: _one_hot_call(calls.get(node_id, 0), taxonomy) for node_id in evidence}
    lca_dist = {node_id: _lca_distribution(counts.get(node_id, {}), taxonomy) for node_id in evidence}
    edge_dist = edge_distributions(edges, evidence)
    graph_methods = {
        "gated_neighbour": resolve(evidence, edges, "gated_neighbour", edge_dist),
        "bayesian_edge": resolve(evidence, edges, "bayesian_edge", edge_dist),
    }
    resolved = {
        "initial_colouring": hard_assignment(calls_dist),
        "probability_sum": evidence,
        "lca": lca_dist,
        **graph_methods,
    }
    _write_evidence_tables(work, evidence, edge_dist, taxonomy)
    sequences = {node_id: sequence for node_id, sequence in read_fasta(fasta)}
    node_taxa = {node_id: [taxon_id for taxon_id in dist if taxon_id != 0] for node_id, dist in evidence.items()}
    edge_taxa = {edge_id: [taxon_id for taxon_id in dist if taxon_id != 0] for edge_id, dist in edge_dist.items()}
    tocumg = export_tocumg(
        root=Path(__file__).resolve().parents[2],
        graph_id=name,
        sequences=sequences,
        edges=edges,
        node_taxa=node_taxa,
        edge_taxa=edge_taxa,
        taxonomy_names={taxon_id: taxonomy.name(taxon_id) for taxon_id in taxonomy.by_id},
        destination=work / "tocumg",
        relative_to=Path.cwd(),
    )
    _store_colour_layer(work, taxonomy, spec)

    metric_rows = []
    for hypothesis in HYPOTHESES:
        metric_rows.append(
            _write_hypothesis(
                hypothesis=hypothesis,
                dataset=name,
                distributions=resolved[hypothesis],
                lengths=lengths,
                taxonomy=taxonomy,
                edges=edges,
                node_truth=node_truth,
                truth_abundance=truth_abundance,
                benchmark=benchmark,
                spec=spec,
                tocumg=tocumg,
                report_path=report_path,
            )
        )
    return metric_rows


def _store_colour_layer(work: Path, taxonomy: Taxonomy, spec: dict) -> None:
    layer = make_layer(
        layer_id="tca-kmer-probability-sum",
        source="kraken2",
        method="probability_sum",
        database="kraken2",
        database_version=taxonomy.version,
        parameters={"aggregation": "probability_sum", "rank": "S", "top_k": spec["top_k"], "min_sim": spec["min_sim"]},
        taxonomy_version=taxonomy.version,
        operation="replace",
    )
    payload = {
        "layer_id": layer.layer_id,
        "source": layer.source,
        "method": layer.method,
        "database": layer.database,
        "database_version": layer.database_version,
        "parameters": layer.parameters,
        "taxonomy_version": layer.taxonomy_version,
        "timestamp": layer.timestamp,
        "operation": layer.operation,
    }
    (work / "colour_layer.yaml").write_text(yaml.safe_dump(payload, sort_keys=False), encoding="utf-8")


def _write_evidence_tables(
    work: Path,
    evidence: dict[str, dict[int, float]],
    edge_dist: dict[str, dict[int, float]],
    taxonomy: Taxonomy,
) -> None:
    node_rows = []
    graph = EvidenceGraph()
    incoming = {}
    for node_id, weights in evidence.items():
        incoming[node_id] = colours_from_weights(weights, evidence_type="kmer", source="kraken2", counts=weights)
        for taxon_id, weight in sorted(weights.items()):
            node_rows.append({"node_id": node_id, "taxon_id": taxon_id, "weight": f"{weight:.6f}", "rank": taxonomy.rank(taxon_id)})
    graph.apply(
        "node",
        incoming,
        make_layer(
            layer_id="observed-nodes",
            source="kraken2",
            method="probability_sum",
            database="kraken2",
            database_version=taxonomy.version,
            parameters={"rank": "S"},
            taxonomy_version=taxonomy.version,
            operation="replace",
        ),
    )
    edge_rows = []
    edge_incoming = {}
    for edge_id, weights in edge_dist.items():
        edge_incoming[edge_id] = colours_from_weights(weights, evidence_type="shared_kmer_mass", source="endpoint_min", counts=weights)
        for taxon_id, weight in sorted(weights.items()):
            edge_rows.append({"edge_id": edge_id, "taxon_id": taxon_id, "weight": f"{weight:.6f}"})
    graph.apply(
        "edge",
        edge_incoming,
        make_layer(
            layer_id="observed-edges",
            source="endpoint_min",
            method="probability_sum",
            database="kraken2",
            database_version=taxonomy.version,
            parameters={"edge_colour": "min_endpoint_mass"},
            taxonomy_version=taxonomy.version,
            operation="replace",
        ),
    )
    write_tsv(work / "node_evidence.tsv", node_rows, ["node_id", "taxon_id", "weight", "rank"])
    write_tsv(work / "edge_evidence.tsv", edge_rows, ["edge_id", "taxon_id", "weight"])


def _write_hypothesis(
    *,
    hypothesis: str,
    dataset: str,
    distributions: dict[str, dict[int, float]],
    lengths: dict[str, int],
    taxonomy: Taxonomy,
    edges: list[dict],
    node_truth: dict[str, int],
    truth_abundance: dict[int, float],
    benchmark: Path,
    spec: dict,
    tocumg: dict,
    report_path: Path,
) -> dict:
    rows, summary = profile_nodes(distributions, lengths, taxonomy)
    scored_nodes = {node_id: dist for node_id, dist in distributions.items() if node_id in node_truth}
    scored_truth = {node_id: node_truth[node_id] for node_id in scored_nodes}
    classification = classification_scores(scored_nodes, scored_truth)
    predicted = _relative_from_profile(rows, float(summary["unclassified_fraction"]))
    abundance = abundance_scores(predicted, truth_abundance)
    metrics = {
        "dataset": dataset,
        "hypothesis": hypothesis,
        "l1": abundance["l1"],
        "bray_curtis": abundance["bray_curtis"],
        "pearson": abundance["pearson"],
        "spearman": abundance["spearman"],
        "precision": classification["precision"],
        "recall": classification["recall"],
        "f1": classification["f1"],
        "accuracy": classification["accuracy"],
        "path_consistency": path_consistency(edges, distributions),
        "unclassified_fraction": summary["unclassified_fraction"],
        "unclassified_bases": summary["unclassified_bases"],
        "n_nodes": classification["n_nodes"],
    }
    folder = benchmark / hypothesis / dataset
    folder.mkdir(parents=True, exist_ok=True)
    write_tsv(folder / "profile.tsv", rows, PROFILE_COLUMNS)
    assignments = []
    for node_id, dist in sorted(distributions.items()):
        if not dist:
            taxon_id, prob = 0, 0.0
        else:
            taxon_id = min(dist, key=lambda key: (-dist[key], key))
            prob = dist[taxon_id]
        assignments.append(
            {
                "node_id": node_id,
                "taxon_id": taxon_id,
                "probability": f"{prob:.6f}",
                "truth_taxon_id": node_truth.get(node_id, ""),
                "length": lengths.get(node_id, 0),
            }
        )
    write_tsv(
        folder / "node_assignments.tsv",
        assignments,
        ["node_id", "taxon_id", "probability", "truth_taxon_id", "length"],
    )
    manifest = {
        "input_graph": INPUT_GRAPH,
        "graph_construction_method": GRAPH_CONSTRUCTION_METHOD,
        "k": 4,
        "assembler": "contigs already present in the dataset bundle",
        "tca_method": "kraken2 k-mer counts rolled to species, aggregation=probability_sum",
        "taxonomic_database": "kraken2 report shipped with the dataset",
        "taxonomy_database_version": spec["report"],
        "classifier": "kraken2",
        "resolver": hypothesis,
        "resolver_parameters": _resolver_parameters(hypothesis),
        "profiling_method": "length-weighted posterior; assigned_reads is 0 because no read-to-graph map is in the bundle",
        "software_version": __version__,
        "graph_note": spec["graph_note"],
        "tocumg": tocumg,
        "metrics": metrics,
        "unclassified_reads": 0,
    }
    (folder / "metrics.yaml").write_text(yaml.safe_dump(manifest, sort_keys=False), encoding="utf-8")
    (folder / "summary.tsv").write_text(
        "key\tvalue\n" + "\n".join(f"{key}\t{value}" for key, value in metrics.items()) + "\n",
        encoding="utf-8",
    )
    return metrics


def _resolver_parameters(hypothesis: str) -> dict:
    if hypothesis == "initial_colouring":
        return {"source": "kraken2 classified taxid", "rollup": "S", "assignment": "hard"}
    if hypothesis == "probability_sum":
        return {"aggregation": "probability_sum", "rollup": "S", "graph": False}
    if hypothesis == "lca":
        return {"aggregation": "lca", "min_fraction": 0.05, "rollup": "S"}
    if hypothesis == "gated_neighbour":
        return {"iterations": 3, "confident": 0.8, "consensus": 0.75, "uncertain_mix": 0.7, "conflict_mix": 0.5}
    if hypothesis == "bayesian_edge":
        return {"iterations": 3, "neighbour_weight": 1.0, "edge_weight": 1.0, "floor": 1e-6}
    return {}


def run_benchmark(
    data_root: Path,
    datasets: list[str],
    output: Path,
    intermediate: Path | None = None,
) -> dict:
    """Run every requested dataset and write summary charts under ``output/summary``.

    Reusable tables go to ``intermediate``, which defaults to ``./intermediate``
    in the current directory rather than a path inside the source tree.
    """
    if intermediate is None:
        intermediate = Path("intermediate")
    intermediate = Path(intermediate)
    rows: list[dict] = []
    for name in datasets:
        if name not in DATASETS:
            raise KeyError(f"unknown dataset: {name}")
        rows.extend(run_dataset(data_root, name, intermediate=intermediate, benchmark=output))
    summary = output / "summary"
    summary.mkdir(parents=True, exist_ok=True)
    columns = [
        "dataset",
        "hypothesis",
        "l1",
        "bray_curtis",
        "pearson",
        "spearman",
        "precision",
        "recall",
        "f1",
        "accuracy",
        "path_consistency",
        "unclassified_fraction",
    ]
    write_tsv(summary / "metrics.tsv", rows, columns)
    write_summary_charts(rows, summary)
    compare = _abundance_compare(data_root, output, datasets)
    if compare:
        write_tsv(summary / "abundance_compare.tsv", compare, ["dataset", "hypothesis", "taxon_id", "name", "truth", "estimated"])
        write_abundance_chart(compare, summary)
    beats = compare_to_initial(rows)
    (summary / "versus_initial.yaml").write_text(yaml.safe_dump(beats, sort_keys=False), encoding="utf-8")
    return {"ok": True, "n_rows": len(rows), "beats": beats}


def compare_to_initial(rows: list[dict]) -> list[dict]:
    """A hypothesis beats the hard colouring only when both L1 and macro F1 improve."""
    initial_l1 = {row["dataset"]: row["l1"] for row in rows if row["hypothesis"] == "initial_colouring"}
    initial_f1 = {row["dataset"]: row["f1"] for row in rows if row["hypothesis"] == "initial_colouring"}
    beats = []
    for row in rows:
        if row["hypothesis"] == "initial_colouring":
            continue
        dataset = row["dataset"]
        better_l1 = row["l1"] < initial_l1[dataset]
        old_f1 = initial_f1[dataset]
        new_f1 = row["f1"]
        better_f1 = old_f1 is not None and new_f1 is not None and new_f1 > old_f1
        beats.append(
            {
                "dataset": dataset,
                "hypothesis": row["hypothesis"],
                "l1": row["l1"],
                "f1": new_f1,
                "initial_l1": initial_l1[dataset],
                "initial_f1": old_f1,
                "beats_initial_l1": better_l1,
                "beats_initial_f1": better_f1,
                "beats_initial": better_l1 and better_f1,
            }
        )
    return beats


def _abundance_compare(data_root: Path, output: Path, datasets: list[str]) -> list[dict]:
    rows = []
    for dataset in datasets:
        abundance: dict[int, float] = {}
        for row in read_csv(data_root / DATASETS[dataset]["abundance"]):
            taxon_id = int(row["taxon_id"])
            abundance[taxon_id] = abundance.get(taxon_id, 0.0) + float(row["genome_abundance"])
        total = sum(abundance.values()) or 1.0
        abundance = {taxon_id: value / total for taxon_id, value in abundance.items()}
        for hypothesis in ABUNDANCE_CHART_HYPOTHESES:
            profile = output / hypothesis / dataset / "profile.tsv"
            if not profile.is_file():
                continue
            for row in read_tsv(profile):
                taxon_id = int(row["taxon_id"])
                if taxon_id == 0:
                    continue
                rows.append(
                    {
                        "dataset": dataset,
                        "hypothesis": hypothesis,
                        "taxon_id": taxon_id,
                        "name": row["name"],
                        "truth": abundance.get(taxon_id, 0.0),
                        "estimated": float(row["relative_abundance"]),
                    }
                )
    return rows

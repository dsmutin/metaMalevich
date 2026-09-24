"""Error properties of nodes, edges, and taxa on scored communities.

Writes TSV summaries and Altair charts under this example's ``ds/`` directory,
and the half-strain charts under ``examples/half_strains/ds/``. Scores are
read from the benchmark and from each example's assignment tables. The script
does not refit a resolver.
"""

from __future__ import annotations

import math
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "examples" / "heldout_genera"))

from community import assembly_graph_from_fastg  # noqa: E402
from metamalevich.tables import read_tsv, write_tsv  # noqa: E402
from metamalevich.taxonomy import parse_kraken_report  # noqa: E402

FOCUS = (
    "initial_colouring",
    "probability_sum",
    "lca",
    "gated_neighbour",
    "bayesian_edge",
    "edge_union",
    "leakage",
    "leakage_flipped",
    "label_drop",
    "logistic_drop",
)


def _out(dataset: str) -> Path:
    if dataset == "half_strains":
        return ROOT / "examples" / "half_strains" / "ds"
    return ROOT / "examples" / "heldout_genera" / "ds"


def _load_taxonomy(dataset: str):
    if dataset == "samovar10":
        report = ROOT / "samovar10" / "kraken2" / "whole_contig.report"
    else:
        report = ROOT / "examples" / dataset / "work" / "reprofile" / "contigs.report"
    return parse_kraken_report(report.read_text(encoding="utf-8"), source="kraken2", version=report.name)


def _genus(taxonomy, taxon_id: int) -> int | None:
    if taxon_id == 0:
        return None
    return taxonomy.ancestor_at_rank(taxon_id, "G")


def _assignments(dataset: str, hypothesis: str) -> list[dict[str, str]]:
    if dataset == "samovar10":
        path = ROOT / "benchmark" / hypothesis / "samovar10" / "node_assignments.tsv"
    else:
        path = ROOT / "examples" / dataset / "work" / "reprofile" / "assignments" / f"{hypothesis}.tsv"
    if not path.is_file():
        raise FileNotFoundError(path)
    return read_tsv(path)


def _edges(dataset: str) -> list[dict]:
    if dataset == "samovar10":
        rows = read_tsv(ROOT / "intermediate" / "samovar10" / "knn_edges.tsv")
        return [
            {
                "edge_id": row["edge_id"],
                "source": row["source"],
                "target": row["target"],
                "weight": float(row["weight"]),
            }
            for row in rows
        ]
    fastg = ROOT / "examples" / dataset / "work" / "megahit" / "assembly.fastg"
    edges, _sequences = assembly_graph_from_fastg(fastg.read_text(encoding="utf-8", errors="replace"))
    return [
        {
            "edge_id": edge["edge_id"],
            "source": edge["source"],
            "target": edge["target"],
            "weight": float(edge.get("weight") or 1.0),
        }
        for edge in edges
    ]


def _evidence(dataset: str) -> dict[str, dict[str, float]]:
    """Per-node entropy, label count, and k-mer mass from the evidence table or counts."""
    if dataset == "samovar10":
        path = ROOT / "intermediate" / "samovar10" / "node_evidence.tsv"
        grouped: dict[str, list[tuple[float, float]]] = defaultdict(list)
        for row in read_tsv(path):
            grouped[row["node_id"]].append((float(row["weight"]), float(row["evidence_count"])))
    else:
        from metamalevich.bench import _evidence_distributions, _group_counts
        from metamalevich.native import kraken_counts

        work = ROOT / "examples" / dataset / "work" / "reprofile"
        counts = _group_counts(read_tsv(work / "kraken_counts.tsv"))
        calls = {row["seq_id"]: int(row["taxon_id"]) for row in read_tsv(work / "kraken_calls.tsv")}
        taxonomy = _load_taxonomy(dataset)
        distributions, raw = _evidence_distributions(counts, taxonomy, calls)
        grouped = {}
        for node_id, weights in distributions.items():
            grouped[node_id] = [
                (probability, raw.get(node_id, {}).get(taxon_id, 0.0))
                for taxon_id, probability in weights.items()
            ]
    features = {}
    for node_id, pairs in grouped.items():
        weights = [weight for weight, _count in pairs if weight > 0]
        entropy = 0.0
        for weight in weights:
            entropy -= weight * math.log(weight)
        features[node_id] = {
            "n_labels": len(weights),
            "entropy": entropy,
            "kmer_mass": sum(count for _weight, count in pairs),
        }
    return features


def _confidence_bin(probability: float) -> str:
    if probability < 0.5:
        return "0-0.5"
    if probability < 0.8:
        return "0.5-0.8"
    if probability < 0.95:
        return "0.8-0.95"
    return "0.95-1"


def _length_bin(length: int, edges: list[int]) -> str:
    if length <= edges[0]:
        return f"Q1:<={edges[0]}"
    if length <= edges[1]:
        return f"Q2:<={edges[1]}"
    if length <= edges[2]:
        return f"Q3:<={edges[2]}"
    return f"Q4:>{edges[2]}"


def _degree_bin(degree: int) -> str:
    if degree == 0:
        return "0"
    if degree <= 2:
        return "1-2"
    if degree <= 8:
        return "3-8"
    return "9+"


def analyze(dataset: str) -> dict[str, list[dict]]:
    """Summarise one community and write its charts. Returns the summary tables."""
    taxonomy = _load_taxonomy(dataset)
    evidence = _evidence(dataset)
    edges = _edges(dataset)
    degree: dict[str, int] = defaultdict(int)
    weight_sum: dict[str, float] = defaultdict(float)
    for edge in edges:
        if edge["weight"] <= 0:
            continue
        degree[edge["source"]] += 1
        weight_sum[edge["source"]] += edge["weight"]
    soft = {
        row["node_id"]: row
        for row in _assignments(dataset, "probability_sum")
        if row.get("truth_taxon_id") not in ("", None)
    }
    lengths = sorted(int(row["length"]) for row in soft.values())
    quartiles = [lengths[len(lengths) // 4], lengths[len(lengths) // 2], lengths[3 * len(lengths) // 4]]
    neighbour_argmax: dict[str, Counter] = defaultdict(Counter)
    for edge in edges:
        target = soft.get(edge["target"])
        if target is None or edge["weight"] <= 0:
            continue
        neighbour_argmax[edge["source"]][int(target["taxon_id"])] += edge["weight"]

    calls = {}
    for hypothesis in FOCUS:
        calls[hypothesis] = {
            row["node_id"]: int(row["taxon_id"])
            for row in _assignments(dataset, hypothesis)
            if row["node_id"] in soft
        }

    node_rows = []
    for node_id, row in soft.items():
        truth = int(row["truth_taxon_id"])
        own = int(row["taxon_id"])
        votes = neighbour_argmax.get(node_id, Counter())
        vote = votes.most_common(1)[0][0] if votes else None
        truth_genus = _genus(taxonomy, truth)
        own_genus = _genus(taxonomy, own)
        if degree[node_id] == 0:
            situation = "isolated"
        elif vote == truth and own != truth:
            situation = "vote_matches_truth"
        elif vote == own and own != truth:
            situation = "vote_matches_wrong_call"
        elif vote == own:
            situation = "vote_matches_call"
        else:
            situation = "vote_other"
        if own == truth:
            kind = "correct"
        elif own == 0:
            kind = "unclassified"
        elif truth_genus is not None and own_genus == truth_genus:
            kind = "same_genus"
        else:
            kind = "other_genus"
        record = {
            "dataset": dataset,
            "node_id": node_id,
            "length": int(row["length"]),
            "length_bin": _length_bin(int(row["length"]), quartiles),
            "degree": degree[node_id],
            "degree_bin": _degree_bin(degree[node_id]),
            "mean_weight": (weight_sum[node_id] / degree[node_id]) if degree[node_id] else 0.0,
            "max_p": float(row["probability"]),
            "confidence_bin": _confidence_bin(float(row["probability"])),
            "n_labels": evidence.get(node_id, {}).get("n_labels", 0),
            "entropy": evidence.get(node_id, {}).get("entropy", 0.0),
            "kmer_mass": evidence.get(node_id, {}).get("kmer_mass", 0.0),
            "situation": situation,
            "error_kind": kind,
            "truth": truth,
            "truth_name": taxonomy.name(truth),
            "call": own,
            "call_name": taxonomy.name(own),
        }
        for hypothesis, table in calls.items():
            guess = table.get(node_id, 0)
            record[hypothesis] = int(guess != truth)
        node_rows.append(record)

    by_confidence = _rate_table(node_rows, "confidence_bin", dataset)
    by_degree = _rate_table(node_rows, "degree_bin", dataset)
    by_length = _rate_table(node_rows, "length_bin", dataset)
    by_situation = _situation_table(node_rows, dataset)
    confusion = _confusion(node_rows, dataset)
    rescue = _rescue(node_rows, dataset)
    edges_summary = _edge_summary(edges, calls["probability_sum"], soft, dataset)
    taxa = _taxon_table(dataset, taxonomy)

    destination = _out(dataset)
    destination.mkdir(parents=True, exist_ok=True)
    _write(destination, dataset, "by_confidence", by_confidence)
    _write(destination, dataset, "by_degree", by_degree)
    _write(destination, dataset, "by_length", by_length)
    _write(destination, dataset, "by_situation", by_situation)
    _write(destination, dataset, "confusion", confusion)
    _write(destination, dataset, "rescue", rescue)
    _write(destination, dataset, "edges", edges_summary)
    _write(destination, dataset, "taxa", taxa)
    _charts(
        destination,
        dataset,
        by_confidence,
        by_degree,
        by_situation,
        confusion,
        rescue,
        taxa,
    )
    tables = {
        "by_confidence": by_confidence,
        "by_degree": by_degree,
        "by_length": by_length,
        "by_situation": by_situation,
        "confusion": confusion,
        "rescue": rescue,
        "edges": edges_summary,
        "taxa": taxa,
        "nodes": node_rows,
    }
    (destination / f"{dataset}_findings.md").write_text(_findings(dataset, tables), encoding="utf-8")
    return tables


def _rate_table(nodes: list[dict], column: str, dataset: str) -> list[dict]:
    rows = []
    groups = sorted({node[column] for node in nodes})
    for hypothesis in FOCUS:
        for group in groups:
            chosen = [node for node in nodes if node[column] == group]
            bases = sum(node["length"] for node in chosen) or 1
            wrong_bases = sum(node["length"] for node in chosen if node[hypothesis])
            rows.append(
                {
                    "dataset": dataset,
                    "hypothesis": hypothesis,
                    "bin": group,
                    "n": len(chosen),
                    "node_error_rate": sum(node[hypothesis] for node in chosen) / len(chosen) if chosen else 0.0,
                    "base_error_rate": wrong_bases / bases,
                }
            )
    return rows


def _situation_table(nodes: list[dict], dataset: str) -> list[dict]:
    rows = []
    for situation in sorted({node["situation"] for node in nodes}):
        chosen = [node for node in nodes if node["situation"] == situation]
        bases = sum(node["length"] for node in chosen) or 1
        wrong = [node for node in chosen if node["error_kind"] != "correct"]
        rows.append(
            {
                "dataset": dataset,
                "situation": situation,
                "n": len(chosen),
                "bases": sum(node["length"] for node in chosen),
                "node_error_rate": len(wrong) / len(chosen) if chosen else 0.0,
                "base_error_rate": sum(node["length"] for node in wrong) / bases,
                "same_genus_bases": sum(node["length"] for node in wrong if node["error_kind"] == "same_genus"),
                "other_genus_bases": sum(node["length"] for node in wrong if node["error_kind"] == "other_genus"),
                "unclassified_bases": sum(node["length"] for node in wrong if node["error_kind"] == "unclassified"),
            }
        )
    return rows


def _confusion(nodes: list[dict], dataset: str) -> list[dict]:
    counter: Counter[tuple[str, str]] = Counter()
    for node in nodes:
        if node["error_kind"] == "correct":
            continue
        counter[(node["truth_name"], node["call_name"])] += node["length"]
    rows = []
    for (truth_name, call_name), bases in counter.most_common(12):
        rows.append(
            {
                "dataset": dataset,
                "truth_name": truth_name,
                "call_name": call_name,
                "bases": bases,
            }
        )
    return rows


def _rescue(nodes: list[dict], dataset: str) -> list[dict]:
    rows = []
    for hypothesis in FOCUS:
        if hypothesis == "probability_sum":
            continue
        rescued = harmed = 0
        rescued_bases = harmed_bases = 0
        for node in nodes:
            before = node["probability_sum"]
            after = node[hypothesis]
            if before and not after:
                rescued += 1
                rescued_bases += node["length"]
            elif after and not before:
                harmed += 1
                harmed_bases += node["length"]
        rows.append(
            {
                "dataset": dataset,
                "hypothesis": hypothesis,
                "rescued_nodes": rescued,
                "harmed_nodes": harmed,
                "rescued_bases": rescued_bases,
                "harmed_bases": harmed_bases,
                "net_bases": rescued_bases - harmed_bases,
            }
        )
    return rows


def _edge_summary(edges: list[dict], calls: dict[str, int], soft: dict[str, dict], dataset: str) -> list[dict]:
    buckets = {
        "same_truth": {"n": 0, "pred_disagree": 0},
        "different_truth": {"n": 0, "pred_disagree": 0},
    }
    for edge in edges:
        left = soft.get(edge["source"])
        right = soft.get(edge["target"])
        if left is None or right is None:
            continue
        key = "same_truth" if int(left["truth_taxon_id"]) == int(right["truth_taxon_id"]) else "different_truth"
        buckets[key]["n"] += 1
        if calls.get(edge["source"]) != calls.get(edge["target"]):
            buckets[key]["pred_disagree"] += 1
    rows = []
    for key, bucket in buckets.items():
        rows.append(
            {
                "dataset": dataset,
                "endpoints": key,
                "n_edges": bucket["n"],
                "prediction_disagree": bucket["pred_disagree"],
                "disagree_fraction": (bucket["pred_disagree"] / bucket["n"]) if bucket["n"] else 0.0,
            }
        )
    return rows


def _taxon_table(dataset: str, taxonomy) -> list[dict]:
    if dataset == "samovar10":
        path = ROOT / "benchmark" / "summary" / "abundance_compare.tsv"
        rows = [row for row in read_tsv(path) if row["dataset"] == "samovar10"]
    else:
        path = ROOT / "examples" / dataset / "work" / "comparison.tsv"
        rows = read_tsv(path)
        return _taxon_from_comparison(rows, dataset)
    out = []
    for row in rows:
        if row["hypothesis"] not in ("initial_colouring", "probability_sum", "logistic_drop", "leakage"):
            continue
        truth = float(row["truth"])
        estimated = float(row["estimated"])
        out.append(
            {
                "dataset": dataset,
                "hypothesis": row["hypothesis"],
                "name": row["name"],
                "truth": truth,
                "estimated": estimated,
                "error": estimated - truth,
            }
        )
    return out


def _taxon_from_comparison(rows: list[dict], dataset: str) -> list[dict]:
    out = []
    keep = {"initial_colouring", "probability_sum", "logistic_drop", "leakage", "kraken2"}
    for row in rows:
        if row["method"] not in keep:
            continue
        truth = float(row["truth"])
        estimated = float(row["estimated"])
        if truth == 0 and estimated == 0:
            continue
        out.append(
            {
                "dataset": dataset,
                "hypothesis": row["method"],
                "name": row["name"],
                "truth": truth,
                "estimated": estimated,
                "error": estimated - truth,
            }
        )
    return out


def _write(destination: Path, dataset: str, name: str, rows: list[dict]) -> None:
    if not rows:
        return
    write_tsv(destination / f"{dataset}_{name}.tsv", rows, list(rows[0].keys()))


def _charts(destination, dataset, by_confidence, by_degree, by_situation, confusion, rescue, taxa) -> None:
    import altair as alt
    import pandas as pd

    confidence = pd.DataFrame(by_confidence)
    degree = pd.DataFrame(by_degree)
    situation = pd.DataFrame(by_situation)
    pairs = pd.DataFrame(confusion)
    harm = pd.DataFrame(rescue)
    species = pd.DataFrame(taxa)
    charts = {
        "confidence": alt.Chart(confidence)
        .mark_bar()
        .encode(
            x=alt.X("bin:N", title="Soft argmax probability", sort=["0-0.5", "0.5-0.8", "0.8-0.95", "0.95-1"]),
            y=alt.Y("base_error_rate:Q", title="Base-weighted error rate"),
            color=alt.Color("hypothesis:N", title="Hypothesis"),
            xOffset="hypothesis:N",
            tooltip=["hypothesis", "bin", "n", "base_error_rate", "node_error_rate"],
        )
        .properties(title=f"{dataset}: error rate by call confidence", width=420, height=260),
        "degree": alt.Chart(degree)
        .mark_bar()
        .encode(
            x=alt.X("bin:N", title="Outgoing degree", sort=["0", "1-2", "3-8", "9+"]),
            y=alt.Y("base_error_rate:Q", title="Base-weighted error rate"),
            color=alt.Color("hypothesis:N", title="Hypothesis"),
            xOffset="hypothesis:N",
            tooltip=["hypothesis", "bin", "n", "base_error_rate"],
        )
        .properties(title=f"{dataset}: error rate by node degree", width=420, height=260),
        "situation": alt.Chart(situation)
        .mark_bar()
        .encode(
            x=alt.X("situation:N", title="Neighbour vote vs soft call"),
            y=alt.Y("bases:Q", title="Bases"),
            color=alt.Color("base_error_rate:Q", title="Base error rate"),
            tooltip=["situation", "n", "bases", "base_error_rate", "same_genus_bases", "unclassified_bases"],
        )
        .properties(title=f"{dataset}: where the neighbour vote sits", width=420, height=260),
    }
    if not pairs.empty:
        pairs["pair"] = pairs["truth_name"] + " → " + pairs["call_name"]
        charts["confusion"] = (
            alt.Chart(pairs)
            .mark_bar()
            .encode(
                x=alt.X("bases:Q", title="Wrong bases"),
                y=alt.Y("pair:N", sort="-x", title="Truth → probability_sum call"),
                tooltip=["truth_name", "call_name", "bases"],
            )
            .properties(title=f"{dataset}: largest species confusions", width=420, height=280)
        )
    if not harm.empty:
        charts["rescue"] = (
            alt.Chart(harm)
            .mark_bar()
            .encode(
                x=alt.X("hypothesis:N", title="Hypothesis", axis=alt.Axis(labelAngle=-40)),
                y=alt.Y("net_bases:Q", title="Rescued minus harmed bases vs probability_sum"),
                tooltip=["hypothesis", "rescued_nodes", "harmed_nodes", "rescued_bases", "harmed_bases"],
            )
            .properties(title=f"{dataset}: graph methods versus the soft profile", width=420, height=260)
        )
    if not species.empty:
        focus = species[species["hypothesis"].isin(["initial_colouring", "probability_sum", "logistic_drop"])]
        charts["taxa"] = (
            alt.Chart(focus)
            .mark_bar()
            .encode(
                x=alt.X("name:N", title="Taxon", axis=alt.Axis(labelAngle=-40)),
                y=alt.Y("error:Q", title="Estimated minus truth"),
                color=alt.Color("hypothesis:N", title="Hypothesis"),
                xOffset="hypothesis:N",
                tooltip=["hypothesis", "name", "truth", "estimated", "error"],
            )
            .properties(title=f"{dataset}: taxon abundance error", width=520, height=280)
        )
    for name, chart in charts.items():
        html = destination / f"{dataset}_{name}.html"
        chart.save(str(html))
        (destination / f"{dataset}_{name}.json").write_text(chart.to_json(), encoding="utf-8")


def _findings(dataset: str, tables: dict) -> str:
    nodes = tables["nodes"]
    n = len(nodes)
    bases = sum(node["length"] for node in nodes) or 1
    wrong = [node for node in nodes if node["error_kind"] != "correct"]
    lines = [
        f"# {dataset} error properties",
        "",
        "Soft calls are `probability_sum`. Base rates weight contig length.",
        f"Scored nodes: {n}. Wrong nodes: {len(wrong)} "
        f"({len(wrong) / n:.3f}). Wrong bases: {sum(node['length'] for node in wrong) / bases:.3f}.",
        "",
        "## Neighbour situation",
        "",
    ]
    for row in tables["by_situation"]:
        lines.append(
            f"- {row['situation']}: n={row['n']}, bases={row['bases']}, "
            f"base error={row['base_error_rate']:.3f}, same-genus wrong bases={row['same_genus_bases']}, "
            f"unclassified bases={row['unclassified_bases']}"
        )
    lines.extend(["", "## Edges", ""])
    for row in tables["edges"]:
        lines.append(
            f"- {row['endpoints']}: {row['n_edges']} edges, "
            f"prediction disagree {row['disagree_fraction']:.3f}"
        )
    lines.extend(["", "## Rescue versus probability_sum", ""])
    for row in sorted(tables["rescue"], key=lambda item: item["net_bases"]):
        lines.append(
            f"- {row['hypothesis']}: rescued {row['rescued_nodes']} ({row['rescued_bases']} bp), "
            f"harmed {row['harmed_nodes']} ({row['harmed_bases']} bp), net {row['net_bases']} bp"
        )
    lines.extend(["", "## Largest confusions", ""])
    for row in tables["confusion"][:8]:
        lines.append(f"- {row['truth_name']} → {row['call_name']}: {row['bases']} bp")
    lines.append("")
    return "\n".join(lines)


def main() -> None:
    """Analyse samovar10, the held-out genera, and the half-strain subset."""
    analyze("samovar10")
    analyze("heldout_genera")
    analyze("half_strains")


if __name__ == "__main__":
    main()

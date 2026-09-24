"""Genus-rank scores for low75 and low75half.

Kraken2 and Kaiju are the read baselines. The graph method labels each
unclassified contig with the nearest Kraken-labelled contig in canonical
4-mer cosine, then assigns that genus to unclassified reads that map onto
the contig. Species rank is not used: the database species is not the
simulated species.

Outputs land in this directory and in ``examples/low75half/ds/``.
"""

from __future__ import annotations

import csv
import importlib.util
import shutil
import subprocess
import sys
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))

from metamalevich.composition_genus import composition_genus_graph  # noqa: E402
from metamalevich.tables import write_tsv  # noqa: E402

TAXDUMP = Path("/mnt/tank/scratch/partition-metagenomics/databases/taxdump/nodes.dmp")
ENGINE = Path("/mnt/tank/scratch/dsmutin/tools/my/samovar/samovar/src/samovar/taxonomy_engine.py")


def _parser():
    """Samovar NCBI walker. A missing rank returns None."""
    if not TAXDUMP.is_file() or TAXDUMP.stat().st_size == 0:
        raise SystemExit(f"missing NCBI nodes.dmp: {TAXDUMP}")
    if not ENGINE.is_file():
        raise SystemExit(f"missing Samovar taxonomy engine: {ENGINE}")
    spec = importlib.util.spec_from_file_location("taxonomy_engine", ENGINE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.NCBITaxonomyParser(str(TAXDUMP))


def _genus(parser, taxon_id: int) -> int:
    if taxon_id <= 0:
        return 0
    return parser.get_ancestor_by_rank(taxon_id, "genus") or 0


def _read_fasta(path: Path) -> dict[str, str]:
    sequences: dict[str, str] = {}
    name = ""
    chunks: list[str] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(">"):
            if name:
                sequences[name] = "".join(chunks).upper()
            name = line[1:].split()[0]
            chunks = []
        else:
            chunks.append(line.strip())
    if name:
        sequences[name] = "".join(chunks).upper()
    return sequences


def _classifier(path: Path) -> dict[str, int]:
    calls = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        parts = line.split("\t")
        if len(parts) < 3:
            raise SystemExit(f"classifier row has fewer than 3 columns: {path}")
        calls[parts[1].split("/")[0]] = 0 if parts[0] == "U" else int(parts[2])
    if not calls:
        raise SystemExit(f"classifier output is empty: {path}")
    return calls


def _map_reads(fasta: Path, fastq: Path, destination: Path) -> dict[str, str]:
    """Primary minimap2 short-read placement. Read name without /1 or /2."""
    minimap = shutil.which("minimap2")
    if not minimap:
        raise SystemExit("required program is not on PATH: minimap2")
    completed = subprocess.run(
        [minimap, "-x", "sr", "-t", "4", str(fasta), str(fastq)],
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        raise SystemExit(completed.stderr.strip() or "minimap2 failed")
    destination.write_text(completed.stdout, encoding="utf-8")
    best: dict[str, tuple[int, str]] = {}
    for line in completed.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) < 10:
            continue
        read_id = parts[0].split("/")[0]
        matches = int(parts[9])
        target = parts[5]
        previous = best.get(read_id)
        if previous is None or matches > previous[0]:
            best[read_id] = (matches, target)
    return {read_id: target for read_id, (_matches, target) in best.items()}


def _l1(predicted: dict[int, float], truth: dict[int, float]) -> float:
    keys = set(predicted) | set(truth)
    return float(sum(abs(predicted.get(key, 0.0) - truth.get(key, 0.0)) for key in keys))


def _fraction(counts: Counter) -> dict[int, float]:
    total = sum(counts.values())
    if total <= 0:
        raise SystemExit("profile has no mass")
    return {taxon_id: count / total for taxon_id, count in counts.items()}


def _truth(example: Path, parser) -> dict[int, float]:
    accessions = {
        row["accession"]: row
        for row in csv.DictReader((example / "accessions.tsv").open(encoding="utf-8"), delimiter="\t")
    }
    counts: Counter = Counter()
    abundance = example / "work" / "abundance.csv"
    for line in abundance.read_text(encoding="utf-8").splitlines()[1:]:
        if not line.strip():
            continue
        accession, raw_count = line.split(",")
        row = accessions.get(accession)
        if row is None:
            raise SystemExit(f"abundance accession is not in the pin table: {accession}")
        counts[int(row["genus_taxid"])] += int(raw_count)
    return _fraction(counts)


def score_example(example_name: str, parser) -> tuple[list[dict], list[dict]]:
    """Return metric rows and the contig error breakdown for one community."""
    example = ROOT / "examples" / example_name
    work = example / "work"
    fasta = work / "megahit" / "graph_nodes.fa"
    fastq = work / "iss" / "initial" / "sample_full_R1.fastq"
    for path in (fasta, fastq, work / "classify" / "kraken2.output", work / "classify" / "kaiju.output"):
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(f"missing input: {path}")
    sequences = _read_fasta(fasta)
    contig_calls = {
        row["seq_id"]: int(row["taxon_id"])
        for row in csv.DictReader((work / "reprofile" / "kraken_calls.tsv").open(encoding="utf-8"), delimiter="\t")
    }
    assignments = [
        row
        for row in csv.DictReader(
            (work / "reprofile" / "assignments" / "initial_colouring.tsv").open(encoding="utf-8"),
            delimiter="\t",
        )
        if row["truth_taxon_id"]
    ]
    labels = {node_id: _genus(parser, contig_calls.get(node_id, 0)) for node_id in sequences}
    transferred = composition_genus_graph(sequences, labels)
    placement = _map_reads(fasta, fastq, work / "reprofile" / "reads_to_contigs.paf")
    truth = _truth(example, parser)
    kraken = _classifier(work / "classify" / "kraken2.output")
    kaiju = _classifier(work / "classify" / "kaiju.output")

    def read_profile(calls: dict[str, int], fill: bool) -> dict[int, float]:
        counts: Counter = Counter()
        for read_id, taxon_id in calls.items():
            genus_id = _genus(parser, taxon_id)
            if fill and genus_id == 0:
                contig_id = placement.get(read_id)
                if contig_id is not None:
                    genus_id = transferred.get(contig_id, 0)
            counts[genus_id] += 1
        return _fraction(counts)

    rows = []
    for method, profile in (
        ("kraken2", read_profile(kraken, False)),
        ("kaiju", read_profile(kaiju, False)),
        ("kraken2_4mer_graph", read_profile(kraken, True)),
        ("kaiju_4mer_graph", read_profile(kaiju, True)),
    ):
        rows.append(
            {
                "dataset": example_name,
                "method": method,
                "rank": "genus",
                "weight": "read_count",
                "l1": f"{_l1(profile, truth):.6g}",
                "unclassified": f"{profile.get(0, 0.0):.6g}",
            }
        )
    error_rows = _contig_errors(example_name, assignments, labels, transferred, parser)
    return rows, error_rows


def _contig_errors(dataset, assignments, labels, transferred, parser) -> list[dict]:
    buckets = {
        "kraken_contig": Counter(),
        "4mer_graph": Counter(),
    }
    for row in assignments:
        length = int(row["length"])
        truth = _genus(parser, int(row["truth_taxon_id"]))
        for method, table in (("kraken_contig", labels), ("4mer_graph", transferred)):
            predicted = table.get(row["node_id"], 0)
            if predicted and predicted == truth:
                kind = "correct"
            elif predicted == 0:
                kind = "unclassified"
            else:
                kind = "wrong_genus"
            buckets[method][kind] += length
    rows = []
    for method, counts in buckets.items():
        total = sum(counts.values()) or 1
        for kind in ("correct", "unclassified", "wrong_genus"):
            rows.append(
                {
                    "dataset": dataset,
                    "method": method,
                    "kind": kind,
                    "bases": counts[kind],
                    "base_share": f"{counts[kind] / total:.6g}",
                }
            )
    return rows


def _charts(metrics: list[dict], errors: list[dict], destination: Path) -> None:
    import altair as alt
    import pandas as pd

    destination.mkdir(parents=True, exist_ok=True)
    frame = pd.DataFrame(metrics)
    frame["l1"] = frame["l1"].astype(float)
    chart = (
        alt.Chart(frame)
        .mark_bar()
        .encode(
            x=alt.X("method:N", title="Method", axis=alt.Axis(labelAngle=-30)),
            y=alt.Y("l1:Q", title="Genus L1 against read truth"),
            color=alt.Color("dataset:N", title="Community"),
            xOffset="dataset:N",
            tooltip=["dataset", "method", "l1", "unclassified"],
        )
        .properties(title="Genus abundance error: read baselines and the 4-mer graph", width=420, height=260)
    )
    chart.save(str(destination / "genus_l1.html"))
    (destination / "genus_l1.json").write_text(chart.to_json(), encoding="utf-8")
    error = pd.DataFrame(errors)
    error["base_share"] = error["base_share"].astype(float)
    share = (
        alt.Chart(error)
        .mark_bar()
        .encode(
            x=alt.X("method:N", title="Contig labelling"),
            y=alt.Y("base_share:Q", title="Share of assembled bases"),
            color=alt.Color("kind:N", title="Genus call"),
            tooltip=["dataset", "method", "kind", "base_share", "bases"],
        )
        .properties(title="Where contig genus errors sit", width=320, height=240)
        .facet(facet=alt.Facet("dataset:N", title="Community"), columns=2)
    )
    share.save(str(destination / "genus_error_share.html"))
    (destination / "genus_error_share.json").write_text(share.to_json(), encoding="utf-8")


def main() -> None:
    """Score both communities and write tables plus Altair charts."""
    parser = _parser()
    metrics: list[dict] = []
    errors: list[dict] = []
    for name in ("low75", "low75half"):
        metric_rows, error_rows = score_example(name, parser)
        metrics.extend(metric_rows)
        errors.extend(error_rows)
    for folder in (ROOT / "examples" / "low75" / "ds", ROOT / "examples" / "low75half" / "ds"):
        write_tsv(folder / "genus_metrics.tsv", metrics, list(metrics[0].keys()))
        write_tsv(folder / "genus_contig_errors.tsv", errors, list(errors[0].keys()))
        _charts(metrics, errors, folder)
    for row in metrics:
        print(f"{row['dataset']}\t{row['method']}\t{row['l1']}")


if __name__ == "__main__":
    main()

"""low75 at ten times the read depth: 1000000 paired fragments.

Genera, accessions, and the Kraken2 and Kaiju indexes are the low75 pin.
Reads, the assembly, and scores are written under this example's ``work/``.
"""

from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PARENT = ROOT / "examples" / "low75"
HELD = ROOT / "examples" / "heldout_genera"
sys.path.insert(0, str(HELD))

import run_example as heldout  # noqa: E402
from community import load_pairs  # noqa: E402

heldout.EXAMPLE = HERE
heldout.WORK = HERE / "work"
heldout.GRAPH_ID = "low75_enriched"
heldout.TOTAL_READS = 1_000_000
heldout.REPORT_PATH = ROOT / "data" / "raw" / "low75_assembly_data_report.jsonl"


def stage_download() -> None:
    """Require the low75 FASTA pin. This example does not download again."""
    pinned = load_pairs(HERE / "accessions.tsv")
    parent = load_pairs(PARENT / "accessions.tsv")
    if pinned != parent:
        raise SystemExit("examples/low75_enriched/accessions.tsv is not the low75 pin")
    if len({row["pair_id"] for row in pinned}) != 75:
        raise SystemExit("low75_enriched accession table must contain 75 genera")
    missing = [
        row["accession"]
        for row in pinned
        if not (ROOT / "data" / "raw" / "fasta" / f"{row['accession']}.fna").is_file()
    ]
    if missing:
        raise SystemExit("missing FASTA for low75_enriched: " + ", ".join(missing))


def stage_databases() -> None:
    """Reuse the low75 indexes. The genomes are the same pin."""
    for name in ("kraken_db", "kaiju_db"):
        dest = HERE / "work" / name
        source = PARENT / "work" / name
        if dest.exists():
            continue
        if not source.is_dir():
            raise SystemExit(f"missing parent index {source}")
        dest.symlink_to(source)
    heldout_databases()


heldout_databases = heldout.stage_databases
heldout.stage_download = stage_download
heldout.stage_databases = stage_databases


if __name__ == "__main__":
    heldout.main()

"""Every other high100 family, still scored at family rank.

One hundred families is even, so this keeps 50. FASTA files come from the
high100 download.
"""

from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
HELD = ROOT / "examples" / "heldout_genera"
FULL = ROOT / "examples" / "high100"
sys.path.insert(0, str(HELD))

import run_example as heldout  # noqa: E402
from community import every_other_strain, load_pairs  # noqa: E402

heldout.EXAMPLE = HERE
heldout.WORK = HERE / "work"
heldout.GRAPH_ID = "half100half"
heldout.SCORE_RANK = "family"
heldout.REPORT_PATH = ROOT / "data" / "raw" / "high100_assembly_data_report.jsonl"


def stage_download() -> None:
    """Require the selected FASTA files from the high100 download."""
    pinned = load_pairs(HERE / "accessions.tsv")
    expected = every_other_strain(load_pairs(FULL / "accessions.tsv"))
    if pinned != expected:
        raise SystemExit("examples/half100half/accessions.tsv is not every other high100 family")
    missing = []
    for row in pinned:
        path = ROOT / "data" / "raw" / "fasta" / f"{row['accession']}.fna"
        if not path.is_file() or path.stat().st_size == 0:
            missing.append(row["accession"])
    if missing:
        raise SystemExit(
            "missing FASTA for half100half; run examples/high100/run_example.py --stage download first: "
            + ", ".join(missing)
        )


heldout.stage_download = stage_download


if __name__ == "__main__":
    heldout.main()

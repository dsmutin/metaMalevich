"""Same pipeline as low75 on every other genus in that accession table.

Seventy-five genera is an odd count, so the even positions keep 38 genera and
drop 37. FASTA files come from the low75 download. This script does not
download again and does not rewrite that manifest.
"""

from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
HELD = ROOT / "examples" / "heldout_genera"
LOW75 = ROOT / "examples" / "low75"
sys.path.insert(0, str(HELD))

import run_example as heldout  # noqa: E402
from community import every_other_strain, load_pairs  # noqa: E402

heldout.EXAMPLE = HERE
heldout.WORK = HERE / "work"
heldout.GRAPH_ID = "low75half"
heldout.REPORT_PATH = ROOT / "data" / "raw" / "low75_assembly_data_report.jsonl"


def stage_download() -> None:
    """Require the selected FASTA files from the low75 download."""
    pinned = load_pairs(HERE / "accessions.tsv")
    expected = every_other_strain(load_pairs(LOW75 / "accessions.tsv"))
    if pinned != expected:
        raise SystemExit("examples/low75half/accessions.tsv is not every other low75 genus")
    missing = []
    for row in pinned:
        path = ROOT / "data" / "raw" / "fasta" / f"{row['accession']}.fna"
        if not path.is_file() or path.stat().st_size == 0:
            missing.append(row["accession"])
    if missing:
        raise SystemExit(
            "missing FASTA for low75half; run examples/low75/run_example.py --stage download first: "
            + ", ".join(missing)
        )


heldout.stage_download = stage_download


if __name__ == "__main__":
    heldout.main()

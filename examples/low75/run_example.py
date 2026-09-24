"""75 genera: one species simulated, a different species in the Kraken2 and Kaiju databases.

Bacteria, archaea, and viruses contribute 25 genera each. Accessions are pinned
in ``accessions.tsv``. FASTA files land in ``data/raw/fasta`` and are not
committed. This download uses ``data/raw/ncbi_low75.zip`` and does not rewrite
the held-out genera manifest.
"""

from __future__ import annotations

import sys
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
HELD = ROOT / "examples" / "heldout_genera"
sys.path.insert(0, str(HELD))

import run_example as heldout  # noqa: E402
from community import load_pairs, sequence_digest  # noqa: E402

heldout.EXAMPLE = HERE
heldout.WORK = HERE / "work"
heldout.GRAPH_ID = "low75"
heldout.REPORT_PATH = ROOT / "data" / "raw" / "low75_assembly_data_report.jsonl"
ZIP_PATH = ROOT / "data" / "raw" / "ncbi_low75.zip"


def stage_download() -> None:
    """Download pinned assemblies if the zip is absent, then validate this table only."""
    rows = load_pairs(HERE / "accessions.tsv")
    if len({row["pair_id"] for row in rows}) != 75:
        raise SystemExit("low75 accession table must contain 75 genera")
    raw = ROOT / "data" / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    log_path = ROOT / "data" / "logs" / "download.log"
    if not ZIP_PATH.is_file() or not zipfile.is_zipfile(ZIP_PATH):
        datasets = heldout._tool("datasets")
        listing = raw / "low75_accessions.txt"
        listing.write_text("".join(f"{row['accession']}\n" for row in rows), encoding="utf-8")
        heldout._run(
            [
                datasets,
                "download",
                "genome",
                "accession",
                "--inputfile",
                str(listing),
                "--include",
                "genome",
                "--filename",
                str(ZIP_PATH),
            ],
            log_path,
        )
    fasta_dir = raw / "fasta"
    fasta_dir.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(ZIP_PATH) as archive:
        report_name = next((name for name in archive.namelist() if name.endswith("assembly_data_report.jsonl")), "")
        if not report_name:
            raise SystemExit("NCBI zip has no assembly_data_report.jsonl")
        report_path = raw / "low75_assembly_data_report.jsonl"
        report_path.write_bytes(archive.read(report_name))
        for row in rows:
            dest = fasta_dir / f"{row['accession']}.fna"
            if dest.is_file() and dest.stat().st_size > 0:
                continue
            genomic = [name for name in archive.namelist() if row["accession"] in name and name.endswith(".fna")]
            member = next((name for name in genomic if name.endswith("_genomic.fna")), genomic[0] if genomic else "")
            if not member:
                raise SystemExit(f"NCBI zip has no FASTA for {row['accession']}")
            dest.write_bytes(archive.read(member))
    _validate(rows, fasta_dir, report_path)


def _validate(rows: list[dict[str, str]], fasta_dir: Path, report_path: Path) -> None:
    reports = {}
    for line in report_path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        import json

        record = json.loads(line)
        accession, _name, _tax = heldout._report_identity(record)
        reports[accession] = record
    manifest = [
        "accession\trepository\tdownload_url\tlocal_path\tfilename_bytes\tsha256\t"
        "download_status\tvalidation_status\tbiological_status\tnotes\tdownload_date"
    ]
    checksums = [f"{heldout._sha256(ZIP_PATH)}  {ZIP_PATH.relative_to(ROOT)}"]
    digests: dict[str, dict[str, str]] = {}
    species: dict[str, dict[str, str]] = {}
    for row in rows:
        record = reports.get(row["accession"])
        if record is None:
            raise SystemExit(f"assembly report has no record for {row['accession']}")
        _accession, name, tax_id = heldout._report_identity(record)
        if not name.startswith(row["genus"]):
            raise SystemExit(f"{row['accession']} organism {name!r} is not genus {row['genus']}")
        if row["species_taxid"] == species.setdefault(row["pair_id"], {}).get({"sim": "db", "db": "sim"}[row["role"]]):
            raise SystemExit(f"{row['pair_id']} uses one species for both roles")
        species[row["pair_id"]][row["role"]] = row["species_taxid"]
        path = fasta_dir / f"{row['accession']}.fna"
        digest = sequence_digest(path)
        digests.setdefault(row["pair_id"], {})[row["role"]] = digest
        file_sha = heldout._sha256(path)
        checksums.append(f"{file_sha}  {path.relative_to(ROOT)}")
        manifest.append(
            "\t".join(
                [
                    row["accession"],
                    "NCBI Assembly",
                    "",
                    str(path.relative_to(ROOT)),
                    str(path.stat().st_size),
                    file_sha,
                    "success",
                    "pass",
                    "pass",
                    f"organism {name}; tax_id {tax_id}; species_taxid {row['species_taxid']}",
                    heldout.date.today().isoformat(),
                ]
            )
        )
    for pair_id, pair_digests in digests.items():
        if pair_digests.get("sim") == pair_digests.get("db"):
            raise SystemExit(f"{pair_id} sim and db assemblies have the same sequence digest")
    manifest_path = ROOT / "data" / "manifests" / "low75_download_manifest.tsv"
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text("\n".join(manifest) + "\n", encoding="utf-8")
    checksum_path = ROOT / "data" / "checksums" / "low75_checksums.txt"
    checksum_path.parent.mkdir(parents=True, exist_ok=True)
    checksum_path.write_text("\n".join(checksums) + "\n", encoding="utf-8")


heldout.stage_download = stage_download


if __name__ == "__main__":
    heldout.main()

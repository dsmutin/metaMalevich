"""Small TSV and CSV readers used by the benchmark."""

from __future__ import annotations

import csv
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    """Read a headered TSV file."""
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], columns: list[str]) -> None:
    """Write rows with a fixed column order."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def read_csv(path: Path) -> list[dict[str, str]]:
    """Read a headered comma-separated file."""
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def read_fasta(path: Path) -> list[tuple[str, str]]:
    """Return ``(record_id, sequence)`` in file order."""
    records: list[tuple[str, str]] = []
    name = ""
    chunks: list[str] = []

    def flush() -> None:
        nonlocal name, chunks
        if name:
            records.append((name, "".join(chunks).upper()))
        name = ""
        chunks = []

    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if line == "":
                continue
            if line.startswith(">"):
                flush()
                name = line[1:].split()[0]
                continue
            chunks.append(line)
    flush()
    return records


def sanitize_sequence(sequence: str) -> str:
    """Map the sequence onto the CFA alphabet."""
    return "".join(base if base in "ACGTN" else "N" for base in sequence.upper()) or "N"

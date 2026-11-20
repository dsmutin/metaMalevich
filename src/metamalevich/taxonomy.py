"""Taxon records parsed from a Kraken-style report or an NCBI TAXDUMP.

Identifiers are integers. Names are display text. A parent link is taken from
report indentation, or from ``nodes.dmp``. The caller names the taxonomy
version.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


RANK_LADDER = ("R", "D", "K", "P", "C", "O", "F", "G", "S")


def rank_key(rank: str) -> tuple[int, int]:
    """Order ranks from coarse to fine. ``S1`` is finer than ``S``."""
    text = rank.strip()
    base = text.rstrip("0123456789") or text
    suffix = text[len(base) :]
    level = RANK_LADDER.index(base) if base in RANK_LADDER else -1
    extra = int(suffix) if suffix.isdigit() else 0
    return (level, extra)


@dataclass(frozen=True)
class Taxon:
    """One taxon under an explicit taxonomy source and version."""

    taxon_id: int
    parent_taxon_id: int | None
    rank: str
    name: str
    taxonomy_source: str
    taxonomy_version: str


class Taxonomy:
    """Hierarchy used by LCA, rank rollup, and profile names."""

    def __init__(self, records: list[Taxon], *, source: str, version: str) -> None:
        self.source = source
        self.version = version
        self.by_id = {record.taxon_id: record for record in records}
        if len(self.by_id) != len(records):
            raise ValueError("duplicate taxon_id in taxonomy")

    def require(self, taxon_id: int) -> Taxon:
        """Return a taxon or raise if the identifier is absent."""
        if taxon_id not in self.by_id:
            raise KeyError(f"taxon_id {taxon_id} is not in taxonomy {self.version}")
        return self.by_id[taxon_id]

    def parent(self, taxon_id: int) -> int | None:
        """Parent identifier, or None at the root."""
        return self.require(taxon_id).parent_taxon_id

    def ancestors(self, taxon_id: int) -> list[int]:
        """Self first, then parents, stopping at the root."""
        chain = []
        seen: set[int] = set()
        current: int | None = taxon_id
        while current is not None:
            if current in seen:
                raise ValueError(f"taxonomy cycle at {current}")
            if current not in self.by_id:
                break
            seen.add(current)
            chain.append(current)
            current = self.by_id[current].parent_taxon_id
        return chain

    def lca(self, taxon_ids: list[int]) -> int:
        """Lowest common ancestor of taxa that exist in this taxonomy.

        The chosen ancestor is the deepest shared node in the parent tree.
        Rank codes are not used, so a rank outside the standard ladder cannot
        hide a finer ancestor.
        """
        present = [taxon_id for taxon_id in taxon_ids if taxon_id in self.by_id and taxon_id != 0]
        if not present:
            raise ValueError("LCA requires at least one known taxon")
        shared = set(self.ancestors(present[0]))
        for taxon_id in present[1:]:
            shared &= set(self.ancestors(taxon_id))
        if not shared:
            raise ValueError("taxa do not share an ancestor")
        return max(shared, key=lambda taxon_id: len(self.ancestors(taxon_id)))

    def ancestor_at_rank(self, taxon_id: int, rank: str) -> int | None:
        """Walk to ``rank`` (for example ``S``).

        A finer rank such as ``S1`` walks up to ``S``. A coarser rank returns
        None so genus-only evidence is not forced onto a species.
        """
        target = rank_key(rank)
        if taxon_id == 0 or taxon_id not in self.by_id:
            return None
        for current in self.ancestors(taxon_id):
            key = rank_key(self.by_id[current].rank)
            if key == target:
                return current
            if key < target:
                return None
        return None

    def name(self, taxon_id: int) -> str:
        """Scientific name, or ``unclassified`` for taxon 0."""
        if taxon_id == 0:
            return "unclassified"
        if taxon_id not in self.by_id:
            return f"taxon:{taxon_id}"
        return self.by_id[taxon_id].name

    def rank(self, taxon_id: int) -> str:
        """Rank code. Unclassified is ``U``."""
        if taxon_id == 0:
            return "U"
        if taxon_id not in self.by_id:
            return "U"
        return self.by_id[taxon_id].rank


def parse_kraken_report(text: str, *, source: str = "kraken2", version: str) -> Taxonomy:
    """Parse a Kraken2 report into parent-linked taxon records.

    ``version`` is required so later colours record which report they used.
    """
    if not version.strip():
        raise ValueError("taxonomy version is required")
    records: list[Taxon] = []
    stack: list[tuple[int, int]] = []
    for line_number, raw in enumerate(text.splitlines(), start=1):
        line = raw.rstrip("\n")
        if line.strip() == "" or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 6:
            raise ValueError(f"report line {line_number} does not have 6 columns")
        rank = parts[3].strip()
        try:
            taxon_id = int(parts[4])
        except ValueError as exc:
            raise ValueError(f"report line {line_number} has a non-integer taxon id") from exc
        name_field = parts[5]
        depth = len(name_field) - len(name_field.lstrip(" "))
        name = name_field.strip()
        while stack and stack[-1][0] >= depth:
            stack.pop()
        parent = stack[-1][1] if stack else None
        records.append(
            Taxon(
                taxon_id=taxon_id,
                parent_taxon_id=parent,
                rank=rank,
                name=name,
                taxonomy_source=source,
                taxonomy_version=version,
            )
        )
        stack.append((depth, taxon_id))
    if not records:
        raise ValueError("taxonomy report is empty")
    return Taxonomy(records, source=source, version=version)


NCBI_RANK_CODES = {
    "superkingdom": "D",
    "domain": "D",
    "realm": "D",
    "kingdom": "K",
    "phylum": "P",
    "class": "C",
    "order": "O",
    "family": "F",
    "genus": "G",
    "species": "S",
}


def _dmp_fields(line: str) -> list[str]:
    """Split one ``*.dmp`` row on its tab, pipe, tab separator."""
    text = line.rstrip("\n")
    if text.endswith("\t|"):
        text = text[:-2]
    return [field.strip() for field in text.split("\t|\t")]


def parse_taxdump(directory: Path, *, version: str | None = None) -> Taxonomy:
    """Parse NCBI ``nodes.dmp`` and ``names.dmp`` from a TAXDUMP directory.

    Ranks get Kraken codes. A rank outside the standard ladder takes the code
    of its nearest standard ancestor plus the number of steps to it (a strain
    under a species is ``S1``), as in a Kraken2 report. Taxon 1 is ``R``.
    ``version`` defaults to the directory name.
    """
    root = Path(directory)
    nodes_path = root / "nodes.dmp"
    names_path = root / "names.dmp"
    missing = [str(path) for path in (nodes_path, names_path) if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"TAXDUMP is missing {', '.join(missing)}")
    label = version if version is not None else root.resolve().name
    if not label.strip():
        raise ValueError("taxonomy version is required")
    parents: dict[int, int | None] = {}
    ranks: dict[int, str] = {}
    with nodes_path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            fields = _dmp_fields(line)
            if len(fields) < 3:
                raise ValueError(f"nodes.dmp line {line_number} has fewer than 3 fields")
            try:
                taxon_id, parent_id = int(fields[0]), int(fields[1])
            except ValueError as exc:
                raise ValueError(f"nodes.dmp line {line_number} has a non-integer taxon id") from exc
            if taxon_id in parents:
                raise ValueError(f"nodes.dmp repeats taxon_id {taxon_id}")
            parents[taxon_id] = None if parent_id == taxon_id else parent_id
            ranks[taxon_id] = fields[2]
    if not parents:
        raise ValueError(f"{nodes_path} is empty")
    names: dict[int, str] = {}
    with names_path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            fields = _dmp_fields(line)
            if len(fields) < 4:
                raise ValueError(f"names.dmp line {line_number} has fewer than 4 fields")
            if fields[3] == "scientific name":
                names[int(fields[0])] = fields[1]
    codes: dict[int, str] = {}
    for start in parents:
        chain: list[int] = []
        current: int | None = start
        while current is not None and current not in codes:
            if current in chain:
                raise ValueError(f"taxonomy cycle at {current}")
            chain.append(current)
            current = parents.get(current)
        for taxon_id in reversed(chain):
            codes[taxon_id] = _rank_code(taxon_id, ranks[taxon_id], codes.get(parents[taxon_id] or -1))
    records = [
        Taxon(
            taxon_id=taxon_id,
            parent_taxon_id=parent,
            rank=codes[taxon_id],
            name=names.get(taxon_id, f"taxon:{taxon_id}"),
            taxonomy_source="ncbi_taxdump",
            taxonomy_version=label,
        )
        for taxon_id, parent in parents.items()
    ]
    return Taxonomy(records, source="ncbi_taxdump", version=label)


def _rank_code(taxon_id: int, rank: str, parent_code: str | None) -> str:
    """Kraken rank code for one node, given its parent's code."""
    if taxon_id == 1:
        return "R"
    if rank in NCBI_RANK_CODES:
        return NCBI_RANK_CODES[rank]
    if parent_code is None:
        return "R1"
    base = parent_code.rstrip("0123456789") or parent_code
    suffix = parent_code[len(base) :]
    return f"{base}{int(suffix) + 1 if suffix else 1}"

"""Taxon records parsed from a Kraken-style report.

Identifiers are integers. Names are display text. A parent link is taken from
report indentation, and the report file name is the taxonomy version unless
the caller sets one.
"""

from __future__ import annotations

from dataclasses import dataclass


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

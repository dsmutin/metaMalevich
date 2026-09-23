"""Multi-label taxonomic colours.

A colour is evidence, not a resolved class. Layers record the method that
produced them. Writing a layer requires an explicit operation.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime, timezone


@dataclass(frozen=True)
class Colour:
    """One taxonomic hypothesis on one graph element."""

    taxon_id: int
    evidence_type: str
    evidence_count: float
    weight: float
    score: float
    confidence: float
    source: str


@dataclass
class ColourLayer:
    """Provenance for one colouring pass."""

    layer_id: str
    source: str
    method: str
    database: str
    database_version: str
    parameters: dict
    taxonomy_version: str
    timestamp: str
    operation: str


def _now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


OPERATIONS = ("replace", "merge", "intersect", "subtract")
AGGREGATIONS = ("count", "weighted_count", "majority", "lca", "probability_sum")


@dataclass
class EvidenceGraph:
    """Node and edge colour tables. Topology is not stored here."""

    node_colours: dict[str, list[Colour]] = field(default_factory=dict)
    edge_colours: dict[str, list[Colour]] = field(default_factory=dict)
    layers: list[ColourLayer] = field(default_factory=list)

    def _table(self, kind: str) -> dict[str, list[Colour]]:
        if kind == "node":
            return self.node_colours
        if kind == "edge":
            return self.edge_colours
        raise ValueError(f"unknown colour target: {kind}")

    def apply(
        self,
        kind: str,
        incoming: dict[str, list[Colour]],
        layer: ColourLayer,
    ) -> None:
        """Combine ``incoming`` with the current table using ``layer.operation``."""
        if layer.operation not in OPERATIONS:
            raise ValueError(f"unknown colour operation: {layer.operation}")
        table = self._table(kind)
        if table and layer.operation == "replace" and not layer.layer_id:
            raise ValueError("replace requires a layer_id")
        keys = set(table) | set(incoming)
        updated: dict[str, list[Colour]] = {}
        for key in keys:
            updated[key] = _combine(table.get(key, []), incoming.get(key, []), layer.operation)
        if kind == "node":
            self.node_colours = updated
        else:
            self.edge_colours = updated
        self.layers.append(layer)

    def distribution(self, kind: str, element_id: str) -> dict[int, float]:
        """Positive weights for one element. Empty when the element has no colours."""
        rows = self._table(kind).get(element_id, [])
        return {row.taxon_id: row.weight for row in rows if row.weight > 0}


def make_layer(
    *,
    layer_id: str,
    source: str,
    method: str,
    database: str,
    database_version: str,
    parameters: dict,
    taxonomy_version: str,
    operation: str,
) -> ColourLayer:
    """Build a layer record with a UTC timestamp."""
    if method not in AGGREGATIONS and method not in {"observed", "resolved"}:
        raise ValueError(f"unknown colour method: {method}")
    return ColourLayer(
        layer_id=layer_id,
        source=source,
        method=method,
        database=database,
        database_version=database_version,
        parameters=dict(parameters),
        taxonomy_version=taxonomy_version,
        timestamp=_now(),
        operation=operation,
    )


def colours_from_weights(
    weights: dict[int, float],
    *,
    evidence_type: str,
    source: str,
    counts: dict[int, float] | None = None,
) -> list[Colour]:
    """Build colour rows. ``weight`` is the stored score for that taxon."""
    rows = []
    total = sum(max(value, 0.0) for value in weights.values())
    for taxon_id, weight in sorted(weights.items()):
        if weight < 0:
            raise ValueError(f"negative weight for taxon {taxon_id}")
        count = float(counts.get(taxon_id, weight)) if counts else float(weight)
        confidence = (weight / total) if total > 0 else 0.0
        rows.append(
            Colour(
                taxon_id=int(taxon_id),
                evidence_type=evidence_type,
                evidence_count=count,
                weight=float(weight),
                score=float(weight),
                confidence=confidence,
                source=source,
            )
        )
    return rows


def _combine(existing: list[Colour], incoming: list[Colour], operation: str) -> list[Colour]:
    if operation == "replace":
        return list(incoming)
    left = {row.taxon_id: row for row in existing}
    right = {row.taxon_id: row for row in incoming}
    if operation == "merge":
        keys = set(left) | set(right)
    elif operation == "intersect":
        keys = set(left) & set(right)
    elif operation == "subtract":
        keys = set(left) - set(right)
    else:
        raise ValueError(f"unknown colour operation: {operation}")
    rows = []
    for taxon_id in sorted(keys):
        if taxon_id in right and operation == "merge":
            base = right[taxon_id]
            if taxon_id in left:
                prior = left[taxon_id]
                rows.append(
                    Colour(
                        taxon_id=taxon_id,
                        evidence_type=base.evidence_type,
                        evidence_count=prior.evidence_count + base.evidence_count,
                        weight=prior.weight + base.weight,
                        score=prior.score + base.score,
                        confidence=base.confidence,
                        source=base.source,
                    )
                )
            else:
                rows.append(base)
        else:
            rows.append(left[taxon_id])
    return rows

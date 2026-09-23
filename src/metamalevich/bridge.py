"""Project colours onto a MetaMetro CFA and compact it to a ToCUMG.

MetaMetro colour sets are integer presence masks. Weights stay in the
metamalevich colour tables. This module does not implement a second
compaction rule.
"""

from __future__ import annotations

import math
import os
import sys
from pathlib import Path

from metamalevich.tables import sanitize_sequence


def metametro_src(root: Path) -> Path:
    """Locate a MetaMetro checkout.

    Search ``METAMETRO_SRC``, then ``external/MetaMetro`` under ``root`` and
    under the current directory.
    """
    candidates = []
    if os.environ.get("METAMETRO_SRC"):
        candidates.append(Path(os.environ["METAMETRO_SRC"]))
    candidates.append(root / "external" / "MetaMetro" / "src")
    candidates.append(Path("external") / "MetaMetro" / "src")
    seen: set[Path] = set()
    for path in candidates:
        try:
            key = path.resolve()
        except OSError:
            key = path
        if key in seen:
            continue
        seen.add(key)
        if (path / "metametro" / "__init__.py").is_file():
            return path
    raise FileNotFoundError(
        "MetaMetro source was not found. Set METAMETRO_SRC or clone it to external/MetaMetro."
    )


def ensure_metametro(root: Path) -> None:
    """Put MetaMetro on ``sys.path`` once."""
    src = str(metametro_src(root))
    if src not in sys.path:
        sys.path.insert(0, src)


def _composition(sequence: str) -> tuple[str, str]:
    length = len(sequence) or 1
    counts = [sequence.count(base) for base in "ACGTN"]
    gc = (counts[2] + counts[1]) / length
    probs = [count / length for count in counts if count]
    entropy = -sum(prob * math.log(prob) for prob in probs)
    return f"{gc:.6f}", f"{entropy:.6f}"


def _display_path(path: Path, relative_to: Path | None) -> str:
    """Path relative to ``relative_to`` when that stays inside the tree."""
    if relative_to is None:
        return str(path)
    try:
        return str(path.resolve().relative_to(relative_to.resolve()))
    except ValueError:
        return str(path)


def export_tocumg(
    *,
    root: Path,
    graph_id: str,
    sequences: dict[str, str],
    edges: list[dict],
    node_taxa: dict[str, list[int]],
    edge_taxa: dict[str, list[int]],
    taxonomy_names: dict[int, str],
    destination: Path,
    relative_to: Path | None = None,
) -> dict:
    """Colour a CFA with ``colour_cfa`` and write CFA plus CDBG directories.

    Returned ``cfa`` and ``cdbg`` paths are relative to ``relative_to`` when
    that root contains the destination. Colouring must leave node and edge
    order unchanged.
    """
    ensure_metametro(root)
    from metametro.contracts.colouring import colour_cfa
    from metametro.converters import cfa_to_cdbg
    from metametro.formats.cdbg import dump_cdbg
    from metametro.formats.cfa import dump_cfa
    from metametro.formats.cfa.model import CfaGraph

    taxon_ids = sorted({taxon for taxa in list(node_taxa.values()) + list(edge_taxa.values()) for taxon in taxa})
    color_of = {taxon_id: index for index, taxon_id in enumerate(taxon_ids)}
    colors = [
        {
            "color_id": str(color_of[taxon_id]),
            "namespace": "taxon",
            "value": taxonomy_names.get(taxon_id, str(taxon_id)),
        }
        for taxon_id in taxon_ids
    ]
    nodes = []
    stored = {}
    for node_id, sequence in sequences.items():
        clean = sanitize_sequence(sequence)
        stored[node_id] = clean
        gc, entropy = _composition(clean)
        nodes.append({"node_id": node_id, "gc": gc, "entropy": entropy})
    edge_rows = []
    for edge in edges:
        edge_rows.append(
            {
                "edge_id": edge["edge_id"],
                "source": edge["source"],
                "target": edge["target"],
                "orientation": edge.get("orientation", "++"),
                "weight": f"{float(edge.get('weight', 1.0)):.6f}",
            }
        )
    before_nodes = [row["node_id"] for row in nodes]
    before_edges = [(row["edge_id"], row["source"], row["target"]) for row in edge_rows]
    graph = CfaGraph(
        metadata={
            "schema_version": "1.0",
            "graph_id": graph_id,
            "graph_type": "knn",
            "contract": "tca_to_tocumg",
            "contract_version": "1.0",
            "features": {
                "node": {"gc": "float", "entropy": "float"},
                "edge": {"orientation": "orientation", "weight": "float"},
            },
            "source": {"format": "metamalevich", "graph_type": "knn"},
        },
        sequences=stored,
        nodes=nodes,
        edges=edge_rows,
        colors=colors or None,
        node_header=["node_id", "gc", "entropy"],
        edge_header=["edge_id", "source", "target", "orientation", "weight"],
    )
    node_colors = {
        node_id: [color_of[taxon_id] for taxon_id in taxa if taxon_id in color_of]
        for node_id, taxa in node_taxa.items()
    }
    edge_colors = {
        edge["edge_id"]: [color_of[taxon_id] for taxon_id in edge_taxa.get(edge["edge_id"], []) if taxon_id in color_of]
        for edge in edges
    }
    coloured = colour_cfa(graph, node_colors, edge_colors, operation="replace", colors=colors or None)
    after_nodes = [row["node_id"] for row in coloured.nodes]
    after_edges = [(row["edge_id"], row["source"], row["target"]) for row in coloured.edges]
    if after_nodes != before_nodes or after_edges != before_edges:
        raise RuntimeError("TCA changed graph topology")
    cdbg = cfa_to_cdbg(coloured)
    cfa_dir = destination / "cfa"
    cdbg_dir = destination / "cdbg"
    dump_cfa(coloured, cfa_dir)
    dump_cdbg(cdbg, cdbg_dir)
    return {
        "cfa": _display_path(cfa_dir, relative_to),
        "cdbg": _display_path(cdbg_dir, relative_to),
        "n_nodes": len(after_nodes),
        "n_edges": len(after_edges),
        "n_colours": len(taxon_ids),
        "topology_unchanged": True,
    }

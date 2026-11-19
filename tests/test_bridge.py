"""Optional ToCUMG export through the installed MetaMetro package."""

from __future__ import annotations

from pathlib import Path

import pytest

from metamalevich.bridge import export_tocumg

pytestmark = pytest.mark.optional


def test_tocumg_colouring_keeps_topology(tmp_path: Path) -> None:
    """colour_cfa then cfa_to_cdbg leaves the node and edge order unchanged."""
    destination = tmp_path / "tocumg"
    result = export_tocumg(
        graph_id="toy",
        sequences={"a": "ACGTACGT", "b": "GGGGCCCC"},
        edges=[{"edge_id": "e1", "source": "a", "target": "b", "orientation": "++", "weight": 0.5}],
        node_taxa={"a": [10], "b": [20]},
        edge_taxa={"e1": [10]},
        taxonomy_names={10: "speciesA", 20: "speciesB"},
        destination=destination,
        relative_to=tmp_path,
    )
    assert result["topology_unchanged"] is True
    assert result["n_nodes"] == 2
    assert result["n_edges"] == 1
    assert (destination / "cfa").is_dir()
    assert (destination / "cdbg").is_dir()
    assert not result["cfa"].startswith("/")

"""Mandatory cache signatures for rebuilt intermediates."""

from __future__ import annotations

from pathlib import Path

import pytest

from metamalevich.bench import cache_is_fresh, cache_signature, write_cache_signature

pytestmark = pytest.mark.mandatory


def test_cache_misses_when_parameters_change(tmp_path: Path) -> None:
    """A new top-k must not reuse edges built for a different top-k."""
    edges = tmp_path / "knn_edges.tsv"
    edges.write_text("edge_id\tsource\ttarget\n", encoding="utf-8")
    first = cache_signature({"top_k": 8, "min_sim": 0.15})
    write_cache_signature(edges, first)
    assert cache_is_fresh(edges, first) is True
    assert cache_is_fresh(edges, cache_signature({"top_k": 4, "min_sim": 0.15})) is False


def test_cache_misses_without_a_sidecar(tmp_path: Path) -> None:
    """An edge list left by an older run is rebuilt until it has a signature."""
    edges = tmp_path / "knn_edges.tsv"
    edges.write_text("edge_id\tsource\ttarget\n", encoding="utf-8")
    assert cache_is_fresh(edges, cache_signature({"top_k": 8})) is False

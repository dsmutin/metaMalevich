"""Mandatory checks for unsupervised dropping of unsupported taxon labels."""

from __future__ import annotations

import pytest

from metamalevich.label_drop import drop_unsupported_labels, false_positive_report, logistic_label_drop

pytestmark = pytest.mark.mandatory


def test_unsupported_taxon_with_a_pure_neighbour_is_dropped() -> None:
    """A minor taxon with no neighbour mass is removed when it is below both cutoffs."""
    distributions = {"a": {1: 0.9, 2: 0.1}, "b": {1: 1.0}}
    original = {"a": {1: 0.9, 2: 0.1}, "b": {1: 1.0}}
    edges = [{"source": "a", "target": "b"}]
    dropped = drop_unsupported_labels(distributions, edges, min_own=0.2, min_neighbour=0.2)
    assert dropped["a"] == {1: 1.0}
    assert 2 not in dropped["a"]
    assert distributions == original


def test_argmax_is_kept_when_neighbours_disagree() -> None:
    """The argmax stays even when own mass and neighbour support are both below the cutoffs."""
    distributions = {"a": {1: 0.55, 2: 0.45}, "b": {2: 1.0}}
    edges = [{"source": "a", "target": "b"}]
    dropped = drop_unsupported_labels(distributions, edges, min_own=0.6, min_neighbour=1.01)
    assert dropped["a"] == {1: 1.0}
    assert distributions["a"] == {1: 0.55, 2: 0.45}


def test_neighbour_support_keeps_a_taxon_below_the_own_cutoff() -> None:
    """Neighbour mass at or above min_neighbour retains a taxon the node itself barely holds."""
    distributions = {"a": {1: 0.8, 2: 0.2}, "b": {2: 1.0}}
    edges = [{"source": "a", "target": "b", "weight": 1.0}]
    dropped = drop_unsupported_labels(distributions, edges, min_own=0.5, min_neighbour=0.3)
    assert dropped["a"][1] == pytest.approx(0.8)
    assert dropped["a"][2] == pytest.approx(0.2)


def test_empty_distribution_stays_empty() -> None:
    """A node with no taxa is returned empty by both filters."""
    assert drop_unsupported_labels({"a": {}}, []) == {"a": {}}
    assert logistic_label_drop({"a": {}}, []) == {"a": {}}


def test_logistic_label_drop_keeps_argmax_on_a_pure_pair() -> None:
    """A two-node pure graph does not delete either node's argmax."""
    distributions = {"a": {1: 1.0}, "b": {1: 1.0}}
    original = {"a": {1: 1.0}, "b": {1: 1.0}}
    edges = [{"source": "a", "target": "b"}, {"source": "b", "target": "a"}]
    dropped = logistic_label_drop(distributions, edges)
    assert dropped["a"] == {1: 1.0}
    assert dropped["b"] == {1: 1.0}
    assert distributions == original


def test_logistic_drops_a_taxon_neighbours_never_choose() -> None:
    """With enough pseudo-labels, a minor taxon that no neighbour argmaxes is removed."""
    distributions = {
        "a": {1: 0.9, 2: 0.1},
        "b": {1: 0.9, 2: 0.1},
        "c": {1: 0.9, 3: 0.1},
        "d": {1: 0.9, 3: 0.1},
    }
    nodes = ("a", "b", "c", "d")
    edges = [{"source": source, "target": target} for source in nodes for target in nodes if source != target]
    dropped = logistic_label_drop(distributions, edges)
    assert set(dropped["a"]) == {1}
    assert dropped["a"][1] == pytest.approx(1.0)
    assert distributions["a"][2] == 0.1


def test_false_positive_report_one_node_mismatch() -> None:
    """One wrong node has a false-positive rate of 1."""
    report = false_positive_report({"a": {1: 1.0}}, {"a": 2})
    assert report["n"] == 1
    assert report["n_false_positive"] == 1
    assert report["false_positive_rate"] == 1
    assert report["n_false_negative"] == 1
    assert report["micro_precision"] == 0


def test_false_positive_report_breaks_argmax_ties_toward_the_smaller_id() -> None:
    """Equal masses select the smaller taxon id before comparing with truth."""
    matched = false_positive_report({"a": {5: 0.5, 2: 0.5}}, {"a": 2})
    missed = false_positive_report({"a": {5: 0.5, 2: 0.5}}, {"a": 5})
    assert matched["false_positive_rate"] == 0
    assert missed["false_positive_rate"] == 1


def test_logistic_threshold_outside_unit_interval() -> None:
    """A threshold outside [0, 1] is rejected. The endpoints are accepted."""
    with pytest.raises(ValueError):
        logistic_label_drop({}, [], threshold=1.5)
    with pytest.raises(ValueError):
        logistic_label_drop({}, [], threshold=-0.1)
    assert logistic_label_drop({}, [], threshold=0.0) == {}
    assert logistic_label_drop({}, [], threshold=1.0) == {}

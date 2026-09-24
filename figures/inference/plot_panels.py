"""Baseline and method panels for L1, Pearson R2, and presence F1.

``read_metrics.tsv`` is written by ``score_metrics.py``. Baselines are Kraken2
and Kaiju. Method panels are the graph transfers tried on top of each baseline.
The high100 assembly-hypothesis panel reads ``examples/high100/work/metrics.tsv``.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "figures" / "inference"
GENUS_ORDER = ["low75", "low75half", "low75_enriched", "low75half_enriched"]
FAMILY_ORDER = ["high100", "high100_enriched", "half100half_enriched"]
BASELINE_ORDER = ["Kraken2", "Kaiju"]
METHOD_ORDER = [
    "Kraken2 + k141",
    "Kraken2 + k21",
    "Kraken2 + 4-mer",
    "Kaiju + k141",
    "Kaiju + k21",
    "Kaiju + 4-mer",
]
BASELINE_COLORS = ["#0072B2", "#E69F00"]
METHOD_COLORS = ["#0072B2", "#56B4E9", "#D55E00", "#009E73", "#F0E442", "#CC79A7"]


def _style() -> None:
    import cnsplots as cns

    cns.settings.font_family = "sans-serif"
    cns.settings.font_sans_serif = ("DejaVu Sans", "Arial", "Helvetica")
    cns.settings.axes_spines_top = False
    cns.settings.axes_spines_right = False
    cns.settings.legend_frameon = False


def _save(stem: Path) -> None:
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig = plt.gcf()
    fig.patch.set_facecolor("white")
    for ax in fig.axes:
        ax.set_facecolor("white")
    for ext, extra in (("pdf", {}), ("svg", {}), ("png", {"dpi": 600})):
        fig.savefig(stem.with_suffix(f".{ext}"), bbox_inches="tight", facecolor="white", transparent=False, **extra)
    plt.close("all")


def _panel_bars(frame: pd.DataFrame, examples: list[str], methods: list[str], colors: list[str], ylabel: str, title: str) -> None:
    import cnsplots as cns

    present = [method for method in methods if method in set(frame["method"])]
    ax = cns.barplot(
        data=frame,
        x="example",
        y="value",
        order=examples,
        hue="method",
        hue_order=present,
        palette=colors[: len(present)],
        legend=True,
    )
    ax.set_xlabel("")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.tick_params(axis="x", labelrotation=20)
    cns.setup_ax(ax)


def _three(frame: pd.DataFrame, examples: list[str], methods: list[str], colors: list[str], stem: Path, title: str) -> None:
    import cnsplots as cns

    _style()
    metrics = [("l1", "L1"), ("r2", "Pearson R2"), ("f1", "Presence F1")]
    mp = cns.multipanel(max_width=980, title=title, title_fontweight="regular")
    for index, (column, label) in enumerate(metrics):
        sub = frame.rename(columns={column: "value"})
        mp.panel(chr(ord("A") + index), width=300, height=220, pad_left=48, pad_top=16, margin_right=8)
        _panel_bars(sub, examples, methods, colors, label, label)
    _save(stem)


def _resolvers() -> None:
    import cnsplots as cns

    path = ROOT / "examples" / "high100" / "work" / "metrics.tsv"
    rows = list(csv.DictReader(path.open(encoding="utf-8"), delimiter="\t"))
    kept = []
    for row in rows:
        if row["weight"] != "node_length" and row["method"] not in {"kraken2", "kaiju"}:
            continue
        if row["method"] in {"kraken2", "kaiju"} and row["weight"] != "read_count":
            continue
        kept.append(row)
    frame = pd.DataFrame(kept)
    frame["l1"] = frame["l1"].astype(float)
    frame["r2"] = frame["r2"].astype(float)
    frame["f1"] = frame["presence_f1"].astype(float)
    frame["group"] = frame["method"].map(lambda name: "baseline" if name in {"kraken2", "kaiju"} else "method")
    order = ["kraken2", "kaiju"] + sorted(
        (name for name in frame["method"] if name not in {"kraken2", "kaiju"}),
        key=lambda name: float(frame.loc[frame["method"] == name, "l1"].iloc[0]),
    )
    _style()
    mp = cns.multipanel(max_width=1100, title="high100 family: baselines and assembly hypotheses", title_fontweight="regular")
    for index, (column, label) in enumerate((("l1", "L1"), ("r2", "Pearson R2"), ("f1", "Presence F1"))):
        sub = frame.rename(columns={column: "value"})
        mp.panel(chr(ord("A") + index), width=340, height=280, pad_left=110, pad_top=14, margin_right=6)
        ax = cns.barplot(
            data=sub,
            x="value",
            y="method",
            order=order,
            hue="group",
            hue_order=["baseline", "method"],
            palette=["#0072B2", "#009E73"],
            legend=True,
        )
        ax.set_xlabel(label)
        ax.set_ylabel("")
        ax.set_title(label)
        cns.setup_ax(ax)
    _save(OUT / "high100_hypotheses")


def main() -> None:
    scores = pd.read_csv(OUT / "read_metrics.tsv", sep="\t")
    genus = scores[scores["rank"] == "genus"]
    family = scores[scores["rank"] == "family"]
    _three(genus[genus["group"] == "baseline"], GENUS_ORDER, BASELINE_ORDER, BASELINE_COLORS, OUT / "genus_baselines", "Genus baselines")
    _three(genus[genus["group"] == "method"], GENUS_ORDER, METHOD_ORDER, METHOD_COLORS, OUT / "genus_methods", "Genus graph methods")
    _three(family[family["group"] == "baseline"], FAMILY_ORDER, BASELINE_ORDER, BASELINE_COLORS, OUT / "family_baselines", "Family baselines")
    _three(family[family["group"] == "method"], FAMILY_ORDER, METHOD_ORDER, METHOD_COLORS, OUT / "family_methods", "Family graph methods")
    _resolvers()
    print(OUT)


if __name__ == "__main__":
    main()

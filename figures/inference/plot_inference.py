"""Publication bars for read-level inference on the key communities.

Numbers are copied from the scored TSV files. Outputs are PDF, SVG, and PNG
under ``figures/inference/``.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "figures" / "inference"

METHOD_ORDER = [
    "Kraken2",
    "Kaiju",
    "Kaiju + k141",
    "Kaiju + k21",
    "Kaiju + 4-mer",
]
METHOD_COLOR = {
    "Kraken2": "#0072B2",
    "Kaiju": "#E69F00",
    "Kaiju + k141": "#009E73",
    "Kaiju + k21": "#56B4E9",
    "Kaiju + 4-mer": "#D55E00",
}
GENUS_ORDER = ["low75", "low75half", "low75_enriched", "low75half_enriched"]
FAMILY_ORDER = ["high100", "high100_enriched", "half100half_enriched"]


def _rows() -> list[dict[str, str]]:
    """Read-level L1 and unclassified fraction already measured for each method."""
    raw = [
        ("low75", "genus", "Kraken2", 1.54014, 0.76599),
        ("low75", "genus", "Kaiju", 1.01624, 0.49412),
        ("low75", "genus", "Kaiju + k21", 0.9688, 0.46864),
        ("low75", "genus", "Kaiju + 4-mer", 0.78304, 0.31074),
        ("low75half", "genus", "Kraken2", 1.56126, 0.77663),
        ("low75half", "genus", "Kaiju", 1.1342, 0.5583),
        ("low75half", "genus", "Kaiju + k21", 1.01604, 0.4983),
        ("low75half", "genus", "Kaiju + 4-mer", 0.75986, 0.24982),
        ("low75_enriched", "genus", "Kraken2", 1.54117, 0.766586),
        ("low75_enriched", "genus", "Kaiju", 1.0121, 0.492109),
        ("low75_enriched", "genus", "Kaiju + k141", 0.785012, 0.374918),
        ("low75_enriched", "genus", "Kaiju + k21", 0.831512, 0.395649),
        ("low75_enriched", "genus", "Kaiju + 4-mer", 0.544312, 0.111869),
        ("low75half_enriched", "genus", "Kraken2", 1.56177, 0.776901),
        ("low75half_enriched", "genus", "Kaiju", 1.13852, 0.560632),
        ("low75half_enriched", "genus", "Kaiju + k141", 0.717654, 0.331246),
        ("low75half_enriched", "genus", "Kaiju + k21", 0.911108, 0.439982),
        ("low75half_enriched", "genus", "Kaiju + 4-mer", 0.42235, 0.044543),
        ("high100", "family", "Kraken2", 1.6323, 0.81317),
        ("high100", "family", "Kaiju", 1.2605, 0.61366),
        ("high100", "family", "Kaiju + k141", 1.24614, 0.60173),
        ("high100", "family", "Kaiju + k21", 1.23956, 0.60137),
        ("high100", "family", "Kaiju + 4-mer", 1.2089, 0.47124),
        ("high100_enriched", "family", "Kraken2", 1.63308, 0.814037),
        ("high100_enriched", "family", "Kaiju", 1.26107, 0.613819),
        ("high100_enriched", "family", "Kaiju + k141", 1.13009, 0.533185),
        ("high100_enriched", "family", "Kaiju + k21", 1.20249, 0.579648),
        ("high100_enriched", "family", "Kaiju + 4-mer", 0.835244, 0.257421),
        ("half100half_enriched", "family", "Kraken2", 1.67935, 0.838908),
        ("half100half_enriched", "family", "Kaiju", 1.3402, 0.663717),
        ("half100half_enriched", "family", "Kaiju + k141", 0.92092, 0.426942),
        ("half100half_enriched", "family", "Kaiju + k21", 1.20726, 0.584535),
        ("half100half_enriched", "family", "Kaiju + 4-mer", 0.635004, 0.118426),
    ]
    return [
        {
            "example": example,
            "rank": rank,
            "method": method,
            "l1": f"{l1:.6g}",
            "unclassified": f"{unclassified:.6g}",
        }
        for example, rank, method, l1, unclassified in raw
    ]


def _density() -> list[dict[str, str]]:
    raw = [
        ("low75", "k141", 2685, 24),
        ("low75", "k21", 122002, 26914),
        ("low75half", "k141", 4119, 296),
        ("low75half", "k21", 108039, 64134),
        ("low75_enriched", "k141", 51535, 1292),
        ("low75_enriched", "k21", 527112, 668818),
        ("low75half_enriched", "k141", 28455, 6578),
        ("low75half_enriched", "k21", 356823, 558598),
        ("high100", "k141", 989, 4),
        ("high100", "k21", 149792, 25674),
        ("high100_enriched", "k141", 58137, 620),
        ("high100_enriched", "k21", 876751, 823034),
        ("half100half_enriched", "k141", 31394, 1238),
        ("half100half_enriched", "k21", 569445, 659448),
    ]
    rows = []
    for example, graph, nodes, edges in raw:
        rows.append(
            {
                "example": example,
                "graph": graph,
                "nodes": str(nodes),
                "edges": str(edges),
                "edges_per_node": f"{edges / nodes:.6g}",
            }
        )
    return rows


def _write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _save(stem: Path) -> None:
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig = plt.gcf()
    fig.patch.set_facecolor("white")
    for ax in fig.axes:
        ax.set_facecolor("white")
    for ext, extra in (("pdf", {}), ("svg", {}), ("png", {"dpi": 600})):
        fig.savefig(
            stem.with_suffix(f".{ext}"),
            bbox_inches="tight",
            facecolor="white",
            transparent=False,
            **extra,
        )
    plt.close("all")


def _bars(frame: pd.DataFrame, examples: list[str], title: str, ylabel: str, stem: Path) -> None:
    import cnsplots as cns

    present = [method for method in METHOD_ORDER if method in set(frame["method"])]
    colors = [METHOD_COLOR[method] for method in present]
    cns.settings.font_family = "sans-serif"
    cns.settings.font_sans_serif = ("DejaVu Sans", "Arial", "Helvetica")
    cns.settings.axes_spines_top = False
    cns.settings.axes_spines_right = False
    cns.settings.legend_frameon = False
    cns.figure(width=720, height=320)
    ax = cns.barplot(
        data=frame,
        x="example",
        y="value",
        order=examples,
        hue="method",
        hue_order=present,
        palette=colors,
        legend=True,
    )
    ax.set_xlabel("")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.tick_params(axis="x", labelrotation=25)
    cns.setup_ax(ax)
    _save(stem)


def main() -> None:
    score_rows = _rows()
    density_rows = _density()
    _write_tsv(OUT / "read_metrics.tsv", score_rows)
    _write_tsv(OUT / "graph_density.tsv", density_rows)
    scores = pd.DataFrame(score_rows)
    scores["l1"] = scores["l1"].astype(float)
    scores["unclassified"] = scores["unclassified"].astype(float)
    genus = scores[scores["rank"] == "genus"]
    family = scores[scores["rank"] == "family"]
    _bars(
        genus.rename(columns={"l1": "value"}),
        GENUS_ORDER,
        "Genus L1 on reads",
        "L1",
        OUT / "genus_l1",
    )
    _bars(
        family.rename(columns={"l1": "value"}),
        FAMILY_ORDER,
        "Family L1 on reads",
        "L1",
        OUT / "family_l1",
    )
    _bars(
        genus.rename(columns={"unclassified": "value"}),
        GENUS_ORDER,
        "Unclassified read fraction at genus",
        "Unclassified fraction",
        OUT / "genus_unclassified",
    )
    _bars(
        family.rename(columns={"unclassified": "value"}),
        FAMILY_ORDER,
        "Unclassified read fraction at family",
        "Unclassified fraction",
        OUT / "family_unclassified",
    )
    density = pd.DataFrame(density_rows)
    density["edges_per_node"] = density["edges_per_node"].astype(float)
    import cnsplots as cns

    cns.figure(width=720, height=320)
    ax = cns.barplot(
        data=density,
        x="example",
        y="edges_per_node",
        hue="graph",
        hue_order=["k141", "k21"],
        palette=["#009E73", "#0072B2"],
        legend=True,
    )
    ax.set_yscale("log")
    ax.set_xlabel("")
    ax.set_ylabel("Edges per node")
    ax.set_title("MEGAHIT graph density")
    ax.tick_params(axis="x", labelrotation=25)
    cns.setup_ax(ax)
    _save(OUT / "graph_density")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()

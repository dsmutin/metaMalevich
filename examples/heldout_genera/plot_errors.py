"""Altair views of where colouring and resolution disagree with truth.

Samovar10 numbers come from the committed benchmark tables. The held-out
community chart is written when ``work/comparison.tsv`` exists.
"""

from __future__ import annotations

import csv
from pathlib import Path


def _frame(path: Path):
    import pandas as pd

    with path.open(encoding="utf-8") as handle:
        return pd.DataFrame(list(csv.DictReader(handle, delimiter="\t")))


def write_samovar10_error_charts(benchmark: Path, destination: Path) -> list[Path]:
    """Plot species abundance error and the Enterococcus split on samovar10."""
    import altair as alt

    table = benchmark / "summary" / "abundance_compare.tsv"
    if not table.is_file():
        raise FileNotFoundError(f"missing samovar10 abundance table: {table}")
    frame = _frame(table)
    frame["truth"] = frame["truth"].astype(float)
    frame["estimated"] = frame["estimated"].astype(float)
    frame["error"] = frame["estimated"] - frame["truth"]
    destination.mkdir(parents=True, exist_ok=True)
    written = []
    error = (
        alt.Chart(frame)
        .mark_bar()
        .encode(
            x=alt.X("name:N", sort="-y", title="Species"),
            y=alt.Y("error:Q", title="Estimated minus truth (contig-base share)"),
            color=alt.Color("hypothesis:N", title="Hypothesis"),
            xOffset="hypothesis:N",
            tooltip=["hypothesis", "name", "truth", "estimated", "error"],
        )
        .properties(title="samovar10 species abundance error", width=520, height=280)
    )
    focus_ids = {"1351", "2774762"}
    focus = frame[frame["taxon_id"].isin(focus_ids)]
    split = (
        alt.Chart(focus)
        .mark_bar()
        .encode(
            x=alt.X("hypothesis:N", title="Hypothesis"),
            y=alt.Y("estimated:Q", title="Estimated share"),
            color=alt.Color("name:N", title="Species"),
            xOffset="name:N",
            tooltip=["hypothesis", "name", "truth", "estimated"],
        )
        .properties(
            title="samovar10 Enterococcus: faecalis absorbs DIV2432",
            width=420,
            height=280,
        )
    )
    for name, chart in (("samovar10_species_error", error), ("samovar10_enterococcus", split)):
        html = destination / f"{name}.html"
        payload = destination / f"{name}.json"
        chart.save(str(html))
        payload.write_text(chart.to_json(), encoding="utf-8")
        written.extend((html, payload))
    return written


def write_heldout_chart(comparison: Path, destination: Path) -> list[Path]:
    """Plot held-out abundance error for every method in the comparison table."""
    import altair as alt

    if not comparison.is_file():
        raise FileNotFoundError(f"missing held-out comparison: {comparison}")
    frame = _frame(comparison)
    for column in ("truth", "estimated"):
        frame[column] = frame[column].astype(float)
    frame["error"] = frame["estimated"] - frame["truth"]
    destination.mkdir(parents=True, exist_ok=True)
    chart = (
        alt.Chart(frame)
        .mark_point(filled=True)
        .encode(
            x=alt.X("truth:Q", title="True read fraction"),
            y=alt.Y("estimated:Q", title="Estimated fraction"),
            color=alt.Color("method:N", title="Method"),
            tooltip=["method", "name", "truth", "estimated", "error"],
        )
        .properties(title="Held-out genera: estimated versus true read fraction", width=420, height=320)
    )
    html = destination / "heldout_abundance.html"
    payload = destination / "heldout_abundance.json"
    chart.save(str(html))
    payload.write_text(chart.to_json(), encoding="utf-8")
    return [html, payload]


def main() -> None:
    """Write the samovar10 charts and the held-out chart when its table exists."""
    root = Path(__file__).resolve().parents[2]
    figures = Path(__file__).resolve().parent / "figures"
    write_samovar10_error_charts(root / "benchmark", figures)
    comparison = Path(__file__).resolve().parent / "work" / "comparison.tsv"
    if comparison.is_file():
        write_heldout_chart(comparison, figures)


if __name__ == "__main__":
    main()

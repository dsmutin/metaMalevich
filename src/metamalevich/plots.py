"""Altair charts for hypothesis benchmarks."""

from __future__ import annotations

from pathlib import Path


def write_summary_charts(rows: list[dict], destination: Path) -> list[str]:
    """Write HTML and Vega-Lite JSON charts. Returns the file names."""
    import altair as alt
    import pandas as pd

    destination.mkdir(parents=True, exist_ok=True)
    frame = pd.DataFrame(rows)
    written = []
    l1 = (
        alt.Chart(frame)
        .mark_bar()
        .encode(
            x=alt.X("hypothesis:N", sort=None, title="Hypothesis", axis=alt.Axis(labelAngle=-40)),
            y=alt.Y("l1:Q", title="L1 abundance error"),
            color=alt.Color("dataset:N", title="Dataset"),
            xOffset="dataset:N",
            tooltip=["dataset", "hypothesis", "l1", "f1", "accuracy"],
        )
        .properties(title="Abundance error by colouring and resolution hypothesis", width=720, height=280)
    )
    f1 = (
        alt.Chart(frame)
        .mark_bar()
        .encode(
            x=alt.X("hypothesis:N", sort=None, title="Hypothesis", axis=alt.Axis(labelAngle=-40)),
            y=alt.Y("f1:Q", title="Macro F1 at species"),
            color=alt.Color("dataset:N", title="Dataset"),
            xOffset="dataset:N",
            tooltip=["dataset", "hypothesis", "f1", "accuracy", "l1"],
        )
        .properties(title="Species assignment by hypothesis", width=720, height=280)
    )
    for name, chart in (("l1", l1), ("f1", f1)):
        html_path = destination / f"{name}.html"
        json_path = destination / f"{name}.json"
        chart.save(str(html_path))
        chart.save(str(json_path))
        written.extend([html_path.name, json_path.name])
    return written


def write_abundance_chart(rows: list[dict], destination: Path) -> None:
    """Compare truth and estimated relative abundance across hypotheses."""
    import altair as alt
    import pandas as pd

    destination.mkdir(parents=True, exist_ok=True)
    frame = pd.DataFrame(rows)
    chart = (
        alt.Chart(frame)
        .mark_circle(size=80)
        .encode(
            x=alt.X("truth:Q", title="True relative abundance"),
            y=alt.Y("estimated:Q", title="Estimated relative abundance"),
            color=alt.Color("hypothesis:N", title="Hypothesis"),
            shape=alt.Shape("dataset:N", title="Dataset"),
            tooltip=["dataset", "hypothesis", "name", "truth", "estimated"],
        )
        .properties(title="Closed-community abundance", width=420, height=320)
    )
    rule = alt.Chart(pd.DataFrame({"truth": [0, 1], "estimated": [0, 1]})).mark_line(strokeDash=[4, 4], color="#888888").encode(
        x="truth:Q", y="estimated:Q"
    )
    (chart + rule).save(str(destination / "abundance.html"))
    (chart + rule).save(str(destination / "abundance.json"))

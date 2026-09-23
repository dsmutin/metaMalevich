"""Pipeline entry kept compatible with the scaffold contract keys."""

from __future__ import annotations

from metamalevich.toy import run_toy


def run_pipeline(input_path: str | None = None) -> dict:
    """Run the toy community, or one dataset directory when ``input_path`` is set.

    The result always contains ``status``, ``ok``, and ``input_path``.
    """
    if input_path:
        from metamalevich.bench import DATASETS, run_benchmark
        from pathlib import Path

        root = Path(input_path)
        names = [name for name, spec in DATASETS.items() if (root / spec["fasta"]).is_file()]
        if not names:
            raise FileNotFoundError(f"no known dataset under {input_path}")
        summary = run_benchmark(root, names, root / "benchmark", intermediate=root / "intermediate")
        beats = [row for row in summary["beats"] if row["beats_initial"]]
        return {
            "status": "ok",
            "ok": True,
            "input_path": input_path,
            "beats_initial": bool(beats),
            "n_rows": summary["n_rows"],
        }
    result = run_toy()
    result["input_path"] = input_path
    return result

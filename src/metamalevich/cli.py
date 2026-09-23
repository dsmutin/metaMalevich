"""Command-line entry for metamalevich."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from metamalevich import __version__
from metamalevich.baseline import run_pipeline


def main(argv: list[str] | None = None) -> int:
    """Parse CLI arguments and run the toy pipeline or a benchmark."""
    parser = argparse.ArgumentParser(
        prog="metamalevich",
        description="Improve metagenomic taxonomic re-profiling with coloured assembly graphs",
    )
    parser.add_argument("--version", action="store_true", help="print version and exit")
    parser.add_argument("-o", "--output", default="-", help="output path or - for stdout")
    parser.add_argument("command", nargs="?", choices=["bench", "pipeline"], help="run the dataset benchmark")
    parser.add_argument("--data-root", default=".", help="directory that contains samovar10 and samovar10_ont1b")
    parser.add_argument("--dataset", action="append", dest="datasets", help="dataset name; repeat to run several")
    parser.add_argument("--benchmark", default="benchmark", help="directory for hypothesis outputs")
    args = parser.parse_args(argv)
    if args.version:
        print(__version__)
        return 0
    if args.command in {"bench", "pipeline"}:
        from metamalevich.bench import DATASETS, run_benchmark

        names = args.datasets or list(DATASETS)
        summary = run_benchmark(Path(args.data_root), names, Path(args.benchmark))
        text = json.dumps({"status": "ok", "ok": True, "input_path": args.data_root, **summary}, indent=2, default=str)
        _emit(text, args.output)
        return 0
    result = run_pipeline(None)
    _emit(json.dumps(result, indent=2), args.output)
    return 0


def _emit(text: str, output: str) -> None:
    if output in {"", "-"}:
        print(text)
    else:
        with open(output, "w", encoding="utf-8") as handle:
            handle.write(text + "\n")


if __name__ == "__main__":
    raise SystemExit(main())

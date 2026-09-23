#!/usr/bin/env python3
"""Run metamalevich baseline on the toy example and check the contract."""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from metamalevich.cli import main  # noqa: E402


def run() -> int:
    """Execute the toy CLI and require a lower L1 than the hard colouring."""
    out = Path(__file__).resolve().parent / "data" / "toy_result.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    code = main(["-o", str(out)])
    if code != 0:
        return code
    payload = json.loads(out.read_text(encoding="utf-8"))
    if payload.get("status") != "ok" or payload.get("ok") is not True or payload.get("beats_initial") is not True:
        print("unexpected payload:", payload, file=sys.stderr)
        return 1
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(run())

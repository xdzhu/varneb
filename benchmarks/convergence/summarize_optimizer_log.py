"""Summarize one ASE optimizer log for the VARNEB acceleration benchmark."""

from __future__ import annotations

import argparse
from datetime import datetime, timedelta
import json
from pathlib import Path
import re


ROW = re.compile(
    r"^\s*(?P<optimizer>\S+):\s+(?P<step>\d+)\s+"
    r"(?P<time>\d\d:\d\d:\d\d)\s+(?P<energy>[-+0-9.eE]+)\s+"
    r"(?P<fmax>[-+0-9.eE]+)\s*$"
)


def parse(path: Path):
    rows = []
    for line in path.read_text(encoding="utf-8").splitlines():
        match = ROW.match(line)
        if match:
            row = match.groupdict()
            rows.append(
                {
                    "optimizer": row["optimizer"],
                    "step": int(row["step"]),
                    "time": row["time"],
                    "energy_eV": float(row["energy"]),
                    "fmax_eV_per_A": float(row["fmax"]),
                }
            )
    return rows


def elapsed_seconds(start: str, end: str) -> float:
    epoch = datetime(2000, 1, 1)
    first = datetime.combine(epoch.date(), datetime.strptime(start, "%H:%M:%S").time())
    last = datetime.combine(epoch.date(), datetime.strptime(end, "%H:%M:%S").time())
    if last < first:
        last += timedelta(days=1)
    return (last - first).total_seconds()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("log", type=Path)
    parser.add_argument("--fmax", type=float, default=0.10)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    rows = parse(args.log)
    if not rows:
        raise SystemExit(f"No optimizer rows found in {args.log}")
    crossing = next((row for row in rows if row["fmax_eV_per_A"] <= args.fmax), None)
    report = {
        "log": str(args.log.resolve()),
        "optimizer": rows[-1]["optimizer"],
        "rows": len(rows),
        "initial_fmax_eV_per_A": rows[0]["fmax_eV_per_A"],
        "minimum_fmax_eV_per_A": min(row["fmax_eV_per_A"] for row in rows),
        "last_fmax_eV_per_A": rows[-1]["fmax_eV_per_A"],
        "target_fmax_eV_per_A": args.fmax,
        "first_crossing_step": None if crossing is None else crossing["step"],
        "seconds_to_first_crossing": None
        if crossing is None
        else elapsed_seconds(rows[0]["time"], crossing["time"]),
        "seconds_observed": elapsed_seconds(rows[0]["time"], rows[-1]["time"]),
        "history": rows,
    }
    text = json.dumps(report, indent=2, ensure_ascii=False) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    print(text, end="")


if __name__ == "__main__":
    main()


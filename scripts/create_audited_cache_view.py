"""Create a non-overwriting symlink view of exactly one audited DFT point set.

An append-only work directory may later acquire new points. This view lets a
historical read-only curvature audit inspect its original point set without
deleting or copying any calculator output. It never launches a calculator.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re


NAME = re.compile(r"eval-\d{6}-[0-9a-f]{12}")


def create_view(workdir: Path, audit: Path, view: Path) -> int:
    source = workdir.resolve(strict=True)
    if not source.is_dir() or view.exists() or view.is_symlink():
        raise ValueError("source must be a directory and view must not exist")
    report = json.loads(audit.read_text(encoding="utf-8"))
    records = report.get("evaluations")
    expected = report.get("n_individually_audited_DFT_points")
    if (not isinstance(records, list) or type(expected) is not int
            or expected <= 0 or len(records) != expected):
        raise ValueError("audit point count is incomplete")
    names = [item.get("directory") for item in records]
    if len(set(names)) != expected or any(
        not isinstance(name, str) or NAME.fullmatch(name) is None for name in names
    ):
        raise ValueError("audit has duplicate or unsafe point names")
    targets = []
    for name in names:
        target = (source / name).resolve(strict=True)
        if (target.parent != source or not target.is_dir()
                or not (target / "result.json").is_file()):
            raise ValueError(f"audited point is missing or outside source: {name}")
        targets.append((name, target))
    view.mkdir(parents=True, exist_ok=False)
    for name, target in targets:
        (view / name).symlink_to(target, target_is_directory=True)
    return len(targets)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument("--view", type=Path, required=True)
    args = parser.parse_args()
    count = create_view(args.workdir, args.audit, args.view)
    print(json.dumps({"status": "audited_symlink_view_created", "points": count,
                      "view": str(args.view)}, sort_keys=True))


if __name__ == "__main__":
    main()

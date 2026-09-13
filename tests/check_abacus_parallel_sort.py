"""Regression check for the ASE-ABACUS global sort-file race workaround."""

from __future__ import annotations

import os
from pathlib import Path
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import vcneb.abacus  # noqa: F401  (installs the adapter patch)


def main() -> None:
    try:
        import ase.io.abacus as abacus_io
    except Exception:
        print("abacus_sort_writer_regression=skipped")
        return
    writer = abacus_io.write_input_stru_sort
    with tempfile.TemporaryDirectory() as name:
        root = Path(name)
        old = Path.cwd()
        os.chdir(root)
        try:
            writer([0, 1, 2, 3])
            if (root / "ase_sort.dat").exists():
                raise SystemExit("identity sort unexpectedly wrote global ase_sort.dat")
            writer([1, 0, 2, 3])
            if (root / "ase_sort.dat").read_text(encoding="utf-8").split() != ["1", "0", "2", "3"]:
                raise SystemExit("non-identity sort writer changed semantics")
        finally:
            os.chdir(old)
    print("abacus_sort_writer_regression=ok")


if __name__ == "__main__":
    main()

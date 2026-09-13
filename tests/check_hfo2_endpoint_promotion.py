"""Regression checks for the HfO2 endpoint promotion gate."""

from __future__ import annotations

import json
from pathlib import Path
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.promote_hfo2_endpoints import promote


def _write_endpoint(path: Path, force: float, stress: float) -> None:
    path.mkdir(parents=True, exist_ok=True)
    (path / "CONTCAR").write_text("valid endpoint placeholder\n", encoding="utf-8")
    (path / "relax_summary.json").write_text(
        json.dumps(
            {
                "returncode": 0,
                "converged": True,
                "natoms": 12,
                "final_natoms": 12,
                "composition": {"Hf": 4, "O": 8},
                "max_generalized_force_eV_per_A": force,
                "max_abs_stress_kbar": stress,
            }
        ),
        encoding="utf-8",
    )


def main() -> None:
    with tempfile.TemporaryDirectory() as root_name:
        root = Path(root_name)
        t = root / "T"
        po = root / "PO"
        target = root / "validation"
        _write_endpoint(t, 0.01, 0.05)
        _write_endpoint(po, 0.015, 0.09)
        report = promote(t, po, target, force_thr=0.02, stress_thr=0.1, natoms=12)
        if report["status"] != "promoted" or not (target / "relaxed_T" / "CONTCAR").exists():
            raise SystemExit("valid endpoints were not promoted")

        _write_endpoint(po, 0.015, 0.2)
        try:
            promote(t, po, root / "rejected", force_thr=0.02, stress_thr=0.1, natoms=12)
        except RuntimeError:
            pass
        else:
            raise SystemExit("over-stressed endpoint passed promotion gate")
    print("hfo2_endpoint_promotion_regression=ok")


if __name__ == "__main__":
    main()

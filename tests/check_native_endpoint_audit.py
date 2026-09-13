"""Regression for the native endpoint threshold auditor."""

from __future__ import annotations

import json
from pathlib import Path
from tempfile import TemporaryDirectory
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.audit_native_endpoint import audit


def main() -> None:
    payload = {
        "returncode": 0,
        "converged": True,
        "final_natoms": 12,
        "final_composition": {"Hf": 4, "O": 8},
        "last_max_atom_force_eV_per_A": 0.012,
        "last_max_stress_kbar": 0.04,
    }
    with TemporaryDirectory(prefix="vcneb-native-audit-") as tmp:
        path = Path(tmp) / "native_relax_summary.json"
        path.write_text(json.dumps(payload), encoding="utf-8")
        report = audit(path, force_thr=0.02, stress_thr=0.1, natoms=12, composition="Hf=4,O=8")
        assert report["status"] == "ok", report
        payload["last_max_stress_kbar"] = 0.2
        path.write_text(json.dumps(payload), encoding="utf-8")
        report = audit(path, force_thr=0.02, stress_thr=0.1, natoms=12, composition="Hf=4,O=8")
        assert report["status"] == "failed"
        assert any("max_stress" in item for item in report["issues"])
    print("native_endpoint_audit_regression=ok")


if __name__ == "__main__":
    main()

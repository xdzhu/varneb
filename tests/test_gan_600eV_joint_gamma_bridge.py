"""Guard the GaN paper panel against cross-cutoff mode comparisons."""

from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path

import pytest

from scripts.plot_gan_joint_mode_600eV import load_report, source_rows


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "benchmarks/numerical_integrity/gan_600eV_joint_gamma_bridge_20260928.json"


def test_same_contract_report_and_source_rows() -> None:
    report, digest = load_report(REPORT)
    assert report["input_contract_sha256"]["INCAR"] == (
        "83ba34d4b14ff4ea641bd74dec26e462ea1cc24b575db79a07138ccdfd552772"
    )
    rows = source_rows(report, digest)
    assert len(rows) == 10
    assert {row["electronic_contract"] for row in rows} == {"VASP 600 eV"}
    assert sum(row["value"] for row in rows if row["panel"] == "b" and row["phase"] == "B4") == pytest.approx(100)


def test_plot_rejects_unverified_status(tmp_path: Path) -> None:
    report = deepcopy(json.loads(REPORT.read_text(encoding="utf-8")))
    report["status"] = "mixed_1000_and_600_eV"
    candidate = tmp_path / "changed.json"
    candidate.write_text(json.dumps(report), encoding="utf-8")
    with pytest.raises(ValueError, match="600-eV"):
        load_report(candidate)

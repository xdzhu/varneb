"""Regression guards for the VASP static-baseline production gate."""

from __future__ import annotations

import hashlib
import runpy
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "validate_vasp_static_baseline_gate.py"


def _inputs(root: Path) -> dict[str, str]:
    values = {"INCAR": "ENCUT = 500\n", "KPOINTS": "Gamma\n", "POTCAR": "licensed bytes\n"}
    for name, content in values.items():
        (root / name).write_text(content, encoding="utf-8")
    return {name: hashlib.sha256((root / name).read_bytes()).hexdigest() for name in values}


def _summary(fingerprints: dict[str, str]) -> dict:
    return {
        "status": "completed",
        "execution_mode": "fixed_initial_endpoint_static_scf",
        "n_images": 7,
        "n_interior_images": 5,
        "evaluated_image_index": 0,
        "git_revision": "8501daa",
        "endpoint_structures": {"initial": {"sha256": "vasp-initial"}},
        "licensed_input_fingerprints": {name: {"sha256": digest} for name, digest in fingerprints.items()},
    }


def _endpoint_gate() -> dict:
    return {"matches": True, "endpoints": {"initial": {"matches": True, "candidate_sha256": "vasp-initial"}}}


def test_vasp_static_gate_accepts_matching_input_basis(tmp_path: Path) -> None:
    validate = runpy.run_path(str(SCRIPT))["validate"]
    fingerprints = _inputs(tmp_path)
    validate(_summary(fingerprints), _endpoint_gate(), initial_dir=tmp_path, git_revision="8501daa")


def test_vasp_static_gate_rejects_changed_potcar_and_unknown_revision(tmp_path: Path) -> None:
    validate = runpy.run_path(str(SCRIPT))["validate"]
    fingerprints = _inputs(tmp_path)
    with pytest.raises(ValueError, match="explicit VCNEB_GIT_REVISION"):
        validate(_summary(fingerprints), _endpoint_gate(), initial_dir=tmp_path, git_revision="remote-sync-unknown")
    (tmp_path / "POTCAR").write_text("different licensed bytes\n", encoding="utf-8")
    with pytest.raises(ValueError, match="POTCAR differs"):
        validate(_summary(fingerprints), _endpoint_gate(), initial_dir=tmp_path, git_revision="8501daa")

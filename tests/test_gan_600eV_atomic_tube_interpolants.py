"""Keep the GaN central tube's post-hoc interpolation diagnosis honest."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.compare_gan_600eV_atomic_tube_interpolants import compare


ROOT = Path(__file__).resolve().parents[1]
FIRST = ROOT / "benchmarks/numerical_integrity/gan_600eV_atomic_tube_central_20260928.json"
REFINED = ROOT / "benchmarks/numerical_integrity/gan_600eV_atomic_tube_refinement_20260928.json"
REPORT = ROOT / "benchmarks/numerical_integrity/gan_600eV_atomic_tube_interpolants_20260928.json"


def test_posthoc_models_do_not_promote_failed_contour_gate() -> None:
    computed = compare(FIRST, REFINED)
    saved = json.loads(REPORT.read_text(encoding="utf-8"))
    assert computed == saved
    assert computed["anchors"] == [5, 6, 7, 8, 11, 14, 17, 18, 19, 20, 21, 22]
    assert computed["predeclared_linear_gate_pass"] is False
    assert computed["contour_certified_by_this_comparison"] is False
    for summary in computed["method_summaries"].values():
        assert summary["n_holdouts"] == 20
        assert summary["maximum_absolute_error_meV_per_GaN"] > 1.0


def test_first_audit_hash_is_required(tmp_path: Path) -> None:
    altered = json.loads(FIRST.read_text(encoding="utf-8"))
    altered["cases"][0]["delta_enthalpy_meV_per_GaN"] += 0.1
    altered_path = tmp_path / "altered_first.json"
    altered_path.write_text(json.dumps(altered), encoding="utf-8")
    with pytest.raises(ValueError, match="comparison inputs"):
        compare(altered_path, REFINED)

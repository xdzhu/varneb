"""Guard the raw-audited local GaN 600-eV joint cut and figure claims."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import pytest


pytest.importorskip("matplotlib")

from scripts.plot_gan_600eV_local_joint_cut import (  # noqa: E402
    DEFAULT_HESSIAN_ROOT,
    DEFAULT_PREFIX,
    DEFAULT_REPORT,
    load,
    model_meV_per_GaN,
    source_rows,
)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_local_model_and_actual_holdouts() -> None:
    report = load(DEFAULT_REPORT, DEFAULT_HESSIAN_ROOT)
    rows = source_rows(report, digest(DEFAULT_REPORT))
    assert len(rows) == 13
    assert sum(row["kind"] == "grid" for row in rows) == 8
    assert sum(row["kind"] == "axial_holdout" for row in rows) == 4
    assert model_meV_per_GaN(report, -0.02, 0) < 0
    assert model_meV_per_GaN(report, 0.02, 0) < 0
    assert model_meV_per_GaN(report, 0, -0.0125) > 0
    assert model_meV_per_GaN(report, 0, 0.0125) > 0
    assert max(abs(row["residual_meV_per_GaN"]) for row in rows
               if row["kind"] == "axial_holdout") == pytest.approx(
                   0.0497377730281912, abs=1e-8)


def test_figure_bundle_does_not_claim_ts_or_full_path_surface() -> None:
    qa_path = DEFAULT_PREFIX.with_name(DEFAULT_PREFIX.name + "_qa.json")
    csv_path = DEFAULT_PREFIX.with_name(DEFAULT_PREFIX.name + "_source_data.csv")
    qa = json.loads(qa_path.read_text(encoding="utf-8"))
    with csv_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 13
    assert qa["predeclared_local_gate_pass"] == {
        "negative_u_positive_v_at_both_scales": True,
        "holdout_model_error_below_0p2_meV_per_GaN": True,
    }
    assert qa["maximum_axial_holdout_error_meV_per_GaN"] < 0.2
    assert qa["TS_certified"] is False
    assert qa["whole_path_2D_surface_certified"] is False
    assert qa["source_sha256"]["source_csv"] == digest(csv_path)
    for suffix in (".pdf", ".svg", ".png"):
        path = DEFAULT_PREFIX.with_suffix(suffix)
        assert qa["exports_sha256"][suffix.lstrip(".")] == digest(path)
    local_tiff = DEFAULT_PREFIX.with_suffix(".tiff")
    if local_tiff.exists():
        assert qa["exports_sha256"]["tiff"] == digest(local_tiff)
    assert b"<text" in DEFAULT_PREFIX.with_suffix(".svg").read_bytes()

"""Keep the GaN sample figure traceable and explicitly non-interpolated."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import pytest


pytest.importorskip("matplotlib")

from scripts.plot_gan_600eV_atomic_tube_samples import (  # noqa: E402
    DEFAULT_COMPONENTS,
    DEFAULT_FIRST,
    DEFAULT_PREFIX,
    DEFAULT_REFINED,
    load,
    source_rows,
)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_source_rows_contain_only_computed_points() -> None:
    first, refined, _ = load(DEFAULT_FIRST, DEFAULT_REFINED, DEFAULT_COMPONENTS)
    rows = source_rows(first, refined, digest(DEFAULT_FIRST), digest(DEFAULT_REFINED))
    assert sum(row["kind"] == "archived_path" for row in rows) == 29
    assert sum(row["kind"] == "archived_path_center" for row in rows) == 18
    assert sum(row["kind"] == "offpath_static" for row in rows) == 28
    assert {abs(float(row["q_atom_A"])) for row in rows
            if row["kind"] == "offpath_static"} == {0.025, 0.05}
    assert {row["source_report_sha256"] for row in rows
            if row["kind"] == "offpath_static"} == {
                digest(DEFAULT_FIRST), digest(DEFAULT_REFINED)}
    assert all(row["case"] and row["raw_outcar_sha256"] for row in rows
               if row["kind"] == "offpath_static")
    assert max(float(row["enthalpy_rel_B4_eV_per_GaN"]) for row in rows
               if row["kind"] == "archived_path") == pytest.approx(0.33849086614368673)


def test_export_bundle_records_failed_contour_gate() -> None:
    qa_path = DEFAULT_PREFIX.with_name(DEFAULT_PREFIX.name + "_qa.json")
    csv_path = DEFAULT_PREFIX.with_name(DEFAULT_PREFIX.name + "_source_data.csv")
    qa = json.loads(qa_path.read_text(encoding="utf-8"))
    with csv_path.open(newline="", encoding="utf-8") as handle:
        exported = list(csv.DictReader(handle))
    assert qa["contour_certified"] is False
    assert qa["max_linear_LOO_error_meV_per_GaN"] > qa[
        "predeclared_LOO_gate_meV_per_GaN"]
    assert qa["panel_map"]["b"].endswith("no surface interpolation")
    assert len(exported) == 75
    assert qa["source_sha256"]["source_csv"] == digest(csv_path)
    for suffix in (".pdf", ".svg", ".png"):
        path = DEFAULT_PREFIX.with_suffix(suffix)
        assert qa["exports_sha256"][suffix.lstrip(".")] == digest(path)
    local_tiff = DEFAULT_PREFIX.with_suffix(".tiff")
    if local_tiff.exists():
        assert qa["exports_sha256"]["tiff"] == digest(local_tiff)
    assert b"<text" in DEFAULT_PREFIX.with_suffix(".svg").read_bytes()

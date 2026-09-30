"""Keep the main-text BTO and GaN mode figures tied to measured DFT evidence.

The test checks the committed figure source tables against independent raw
audit records.  A smooth contour is not counted as additional DFT data.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT / "paper" / "VARNEB_CPC"
FIGURES = PAPER / "figures"
AUDITS = ROOT / "benchmarks" / "numerical_integrity"


def _json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as source:
        return list(csv.DictReader(source))


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _assert_complete_grid(
    rows: list[dict[str, str]], x: str, y: str, nx: int, ny: int
) -> None:
    coordinates = [(float(row[x]), float(row[y])) for row in rows]
    assert len(coordinates) == len(set(coordinates)) == nx * ny
    xs, ys = {point[0] for point in coordinates}, {point[1] for point in coordinates}
    assert len(xs) == nx and len(ys) == ny
    assert set(coordinates) == {(x_value, y_value) for x_value in xs for y_value in ys}


def test_bto_main_text_grid_is_289_measured_statics() -> None:
    stem = "bto_frozen_soft_mode_289_20260929"
    csv_path = FIGURES / f"{stem}_source_data.csv"
    qa = _json(FIGURES / f"{stem}_qa.json")
    audit = _json(AUDITS / "bto_transverse_soft_frozen289_20260929.json")
    rows = _rows(csv_path)
    measured = [row for row in rows if row["record"] == "frozen_DFT"]
    path = [row for row in rows if row["record"] == "VCNEB_image"]

    assert (PAPER / "varneb_CPC.tex").read_text(encoding="utf-8").count(
        f"{{figures/{stem}.pdf}}"
    ) == 1
    assert (FIGURES / f"{stem}.pdf").is_file()
    assert qa["interpolation"]["n_measured_DFT_samples"] == 289
    assert qa["source_sha256"]["analysis289"] == _sha256(
        AUDITS / "bto_transverse_soft_frozen289_20260929.json"
    )
    assert audit["status"] == "BTO_frozen_soft_mode_17x17_raw_audited"
    assert audit["n_reused_measured_points"] == 81
    assert audit["n_new_raw_audited_statics"] == len(audit["new_samples"]) == 208
    assert len(measured) == 289
    assert len(path) == 7  # projections of actual VCNEB images, not grid statics
    _assert_complete_grid(
        measured, "q_parallel_sqrt_amu_A", "q_transverse_sqrt_amu_A", 17, 17
    )

    by_coordinate = {
        (float(row["q_parallel_sqrt_amu_A"]), float(row["q_transverse_sqrt_amu_A"])):
        float(row["E_minus_C_eV_per_BTO"])
        for row in measured
    }
    for old_point in _rows(FIGURES / "bto_frozen_soft_mode_81_20260928_source_data.csv"):
        if old_point["record"] != "frozen_DFT":
            continue
        coordinate = (
            float(old_point["q_parallel_sqrt_amu_A"]),
            float(old_point["q_transverse_sqrt_amu_A"]),
        )
        assert math.isclose(
            by_coordinate[coordinate], float(old_point["E_minus_C_eV_per_BTO"]),
            abs_tol=1e-9,
        )
    for new_point in audit["new_samples"]:
        coordinate = (
            float(new_point["q1"]), float(new_point["q2"]),
        )
        assert math.isclose(
            by_coordinate[coordinate],
            new_point["energy_minus_C_eV_per_BTO"],
            abs_tol=1e-9,
        )
    assert math.isclose(
        qa["interpolation"]["prospective_208_nested_node_max_abs_error_meV_per_BTO"],
        audit["prospective_max_abs_error_meV_per_BTO"],
        abs_tol=1e-9,
    )
    assert audit["prospective_max_abs_error_meV_per_BTO"] < 0.5
    assert qa["interpolation"]["interpolated_downward_overshoot_meV_per_BTO"] > 0


def test_gan_main_local_grid_is_289_measured_statics() -> None:
    stem = "gan_600eV_local_joint_dft289_20260930_v2"
    csv_path = FIGURES / f"{stem}_source_data.csv"
    qa = _json(FIGURES / f"{stem}_qa.json")
    audit_path = AUDITS / "gan_600eV_ts_2d_17x17_refinement_20260930.json"
    audit = _json(audit_path)
    rows = _rows(csv_path)

    assert (PAPER / "varneb_CPC.tex").read_text(encoding="utf-8").count(
        f"{{figures/{stem}.pdf}}"
    ) == 1
    assert (FIGURES / f"{stem}.pdf").is_file()
    assert _sha256(csv_path) == qa["source_sha256"]["source_data"]
    assert _sha256(audit_path) == qa["source_sha256"]["dense17_audit"]
    assert audit["status"] == "GaN_600eV_local_joint_17x17_refinement_raw_audited"
    assert qa["pressure_GPa"] == audit["pressure_GPa"] == 45.7
    assert "not" in qa["claim_limit"].lower()
    assert audit["n_prior_DFT_points"] == 81
    assert audit["n_new_DFT_points"] == len(audit["cases"]) == 208
    assert audit["n_total_DFT_points"] == qa["n_measured_DFT_points"] == len(rows) == 289
    _assert_complete_grid(rows, "q_u_A", "q_v_A", 17, 17)
    assert {row["kind"] for row in rows} == {
        "center", "grid", "axial_holdout", "offaxis_refinement",
        "dense9_refinement", "dense17_midpoint"
    }

    old_rows = _rows(FIGURES / "gan_600eV_local_joint_dft81_v2_20260928_source_data.csv")
    by_case = {row["case"]: row for row in rows}
    for old_row in old_rows:
        assert by_case[old_row["case"]] == old_row

    new_rows = {row["case"]: row for row in rows if row["kind"] == "dense17_midpoint"}
    assert len(new_rows) == 208
    for new_point in audit["cases"]:
        row = new_rows[new_point["case"]]
        assert row["raw_OUTCAR_sha256"] == new_point["outcar_sha256"]
        assert math.isclose(
            float(row["measured_delta_H_meV_per_GaN"]),
            new_point["delta_enthalpy_meV_per_GaN"],
            abs_tol=1e-9,
        )
    assert math.isclose(
        qa["prospective_max_abs_error_meV_per_GaN"],
        audit["prior_9x9_prospective_max_abs_error_meV_per_GaN"],
        abs_tol=1e-9,
    )
    assert audit["prior_9x9_prospective_gate"]["pass"] is True

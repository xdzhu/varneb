"""Keep the two 81-point CPC mode figures tied to measured DFT evidence.

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


def _assert_complete_nine_by_nine(rows: list[dict[str, str]], x: str, y: str) -> None:
    coordinates = [(float(row[x]), float(row[y])) for row in rows]
    assert len(coordinates) == len(set(coordinates)) == 81
    xs, ys = {point[0] for point in coordinates}, {point[1] for point in coordinates}
    assert len(xs) == len(ys) == 9
    assert set(coordinates) == {(x_value, y_value) for x_value in xs for y_value in ys}


def test_bto_main_text_grid_is_81_measured_statics() -> None:
    stem = "bto_frozen_soft_mode_81_20260928"
    csv_path = FIGURES / f"{stem}_source_data.csv"
    qa = _json(FIGURES / f"{stem}_qa.json")
    audit = _json(AUDITS / "bto_transverse_soft_frozen81_20260928.json")
    rows = _rows(csv_path)
    measured = [row for row in rows if row["record"] == "frozen_DFT"]
    path = [row for row in rows if row["record"] == "VCNEB_image"]

    assert (PAPER / "varneb_CPC.tex").read_text(encoding="utf-8").count(
        f"{{figures/{stem}.pdf}}"
    ) == 1
    assert (FIGURES / f"{stem}.pdf").is_file()
    assert _sha256(csv_path) == qa["source_data_csv_sha256"]
    assert qa["interpolation"]["n_measured_DFT_samples"] == 81
    assert audit["n_prior_DFT_points"] == 59
    assert audit["n_new_DFT_points"] == len(audit["new_holdouts"]) == 22
    assert audit["n_total_DFT_points"] == len(measured) == 81
    assert len(path) == 7  # projections of actual VCNEB images, not grid statics
    _assert_complete_nine_by_nine(
        measured, "q_parallel_sqrt_amu_A", "q_transverse_sqrt_amu_A"
    )

    by_coordinate = {
        (float(row["q_parallel_sqrt_amu_A"]), float(row["q_transverse_sqrt_amu_A"])):
        float(row["E_minus_C_eV_per_BTO"])
        for row in measured
    }
    for new_point in audit["new_holdouts"]:
        coordinate = (
            float(new_point["q_parallel_sqrt_amu_A"]),
            float(new_point["q_transverse_sqrt_amu_A"]),
        )
        assert math.isclose(
            by_coordinate[coordinate],
            new_point["energy_eV_per_BTO"] - audit["reference_energy_eV_per_BTO"],
            abs_tol=1e-9,
        )
    assert math.isclose(
        qa["interpolation"]["independent_22_missing_node_max_abs_error_meV_per_BTO"],
        audit["prospective_22_point_max_abs_error_meV_per_BTO"],
        abs_tol=1e-9,
    )
    assert qa["interpolation"]["interpolated_downward_overshoot_meV_per_BTO"] > 0


def test_gan_main_local_grid_is_81_measured_statics() -> None:
    stem = "gan_600eV_local_joint_dft81_v2_20260928"
    csv_path = FIGURES / f"{stem}_source_data.csv"
    qa = _json(FIGURES / f"{stem}_qa.json")
    audit = _json(AUDITS / "gan_600eV_ts_2d_9x9_refinement_20260928.json")
    rows = _rows(csv_path)

    assert (PAPER / "varneb_CPC.tex").read_text(encoding="utf-8").count(
        f"{{figures/{stem}.pdf}}"
    ) == 1
    assert (FIGURES / f"{stem}.pdf").is_file()
    assert _sha256(csv_path) == qa["source_sha256"]["source_csv"]
    assert audit["status"] == "GaN_600eV_local_joint_9x9_refinement_raw_audited"
    assert qa["pressure_GPa"] == audit["pressure_GPa"] == 45.7
    assert qa["TS_certified"] is False
    assert qa["whole_path_2D_surface_certified"] is False
    assert audit["n_prior_DFT_points"] == 25
    assert audit["n_new_DFT_points"] == len(audit["cases"]) == 56
    assert audit["n_total_DFT_points"] == len(rows) == 81
    _assert_complete_nine_by_nine(rows, "q_u_A", "q_v_A")
    assert {row["kind"] for row in rows} == {
        "center", "grid", "axial_holdout", "offaxis_refinement", "dense9_refinement"
    }

    new_rows = {row["case"]: row for row in rows if row["kind"] == "dense9_refinement"}
    assert len(new_rows) == 56
    for new_point in audit["cases"]:
        row = new_rows[new_point["case"]]
        assert row["raw_OUTCAR_sha256"] == new_point["outcar_sha256"]
        assert math.isclose(
            float(row["measured_delta_H_meV_per_GaN"]),
            new_point["delta_enthalpy_meV_per_GaN"],
            abs_tol=1e-9,
        )
    assert math.isclose(
        qa["prospective_5x5_max_abs_error_meV_per_GaN"],
        audit["prior_5x5_prospective_max_abs_error_meV_per_GaN"],
        abs_tol=1e-9,
    )

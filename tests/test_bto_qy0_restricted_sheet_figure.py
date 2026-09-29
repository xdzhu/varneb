"""Source-bound checks for the BTO manuscript sheet; no DFT is launched."""

import csv
import hashlib
import json
from pathlib import Path

import numpy as np

from scripts.plot_bto_qy0_restricted_sheet import evidence, same_structure_up_to_translation


ROOT = Path(__file__).resolve().parents[1]


def test_main_text_27_node_sheet_is_measured_and_stress_gated():
    stem = "bto_qy0_restricted_sheet_27_20260929"
    figures = ROOT / "paper/VARNEB_CPC/figures"
    report = json.loads((ROOT / "benchmarks/numerical_integrity/bto_qy0_27_node_sheet_20260929.json").read_text(encoding="utf-8"))
    qa = json.loads((figures / f"{stem}_qa.json").read_text(encoding="utf-8"))
    source = figures / f"{stem}_source_data.csv"
    with source.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    manuscript = (ROOT / "paper/VARNEB_CPC/varneb_CPC.tex").read_text(encoding="utf-8")

    assert f"{{figures/{stem}.pdf}}" in manuscript
    assert (figures / f"{stem}.pdf").is_file()
    assert qa["source_sha256"]["source_data"] == hashlib.sha256(source.read_bytes()).hexdigest()
    assert qa["n_measured_DFT_nodes"] == report["n_total_measured_nodes"] == len(rows) == 27
    assert report["n_reused_measured_nodes"] == 12
    assert report["n_new_raw_audited_nodes"] == len(report["new"]) == 15
    assert report["n_new_raw_audited_DFT_evaluations"] == 215
    assert report["new_node_model_gate_pass"] is True
    assert report["new_node_model_max_abs_error_meV_per_BTO"] < 2.0
    assert {row["source"] for row in rows} == {"audited_DFT"}
    assert np.asarray(report["energy_minus_C_meV_per_BTO"]).shape == (9, 3)
    coordinates = {(float(row["Qz_sqrt_amu_A"]), float(row["Qx_sqrt_amu_A"])) for row in rows}
    assert coordinates == {
        (qz, qx) for qz in report["grid_qz_sqrt_amu_A"]
        for qx in report["grid_qx_sqrt_amu_A"]
    }
    measured = {
        (float(row["Qz_sqrt_amu_A"]), float(row["Qx_sqrt_amu_A"])):
        float(row["E_minus_C_meV_per_BTO"]) for row in rows
    }
    for node in report["new"]:
        assert np.isclose(measured[tuple(node["q_sqrt_amu_A"])],
                          node["energy_minus_C_meV_per_BTO"], atol=1e-9)
        assert node["orthogonal_gradient_eV_per_sqrt_amu_A"] < 0.003
        assert node["maximum_absolute_stress_kbar"] < 2.0


def test_bto_sheet_has_audited_blind_predictions_and_endpoints():
    qa, points, path, t_qz, t_qy = evidence()
    assert qa["status"] == "restricted_Qy0_local_model_two_prospective_gates_pass_not_global_PES"
    assert qa["model_training_nodes"] == 9
    assert qa["prospective_holdout_nodes"] == 2
    assert qa["maximum_absolute_prospective_error_meV_per_BTO"] < qa["gate_meV_per_BTO"]
    assert np.isclose(t_qz, 1.2043298776588882)
    assert abs(t_qy) < 1e-8
    assert abs(path[-1]["q_parallel_sqrt_amu_A"]) < 1e-8
    assert [p["class"] for p in points].count("fit_node") == 9
    assert [p["class"] for p in points].count("retrospective") == 4
    assert [p["class"] for p in points].count("prospective") == 2
    assert [p["class"] for p in points].count("T_endpoint_retrospective") == 1


def test_endpoint_identity_rejects_nonrigid_shift():
    from ase import Atoms

    atoms = Atoms("H2", positions=[[0, 0, 0], [0.5, 0.5, 0.5]],
                  cell=np.eye(3), pbc=True)
    record = {"n_atoms": 2, "species_order": ["H", "H"], "cell_A": np.eye(3).tolist(),
              "fractional_positions_wrapped": [[0.1, 0.1, 0.1], [0.6, 0.6, 0.6]]}
    assert same_structure_up_to_translation(atoms, record)
    record["fractional_positions_wrapped"][1][0] += 0.01
    assert not same_structure_up_to_translation(atoms, record)

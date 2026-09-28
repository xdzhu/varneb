"""Source-bound checks for the BTO manuscript sheet; no DFT is launched."""

import numpy as np

from scripts.plot_bto_qy0_restricted_sheet import evidence, same_structure_up_to_translation


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

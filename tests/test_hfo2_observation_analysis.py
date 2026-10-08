import json
from pathlib import Path
import shutil

import numpy as np
import pytest
from ase.io import read

from scripts.analyze_hfo2_chain_observations import analyze, analyze_observation
from scripts.analyze_hfo2_parent_patterns import rotated_t_triplet


CASE = Path(__file__).resolve().parents[1] / "benchmarks/hfo2_channels/20261008"


def test_actual_HF_observations_replay_without_DFT(tmp_path):
    expected = json.loads((CASE / "chain_observations/analysis.json").read_text())
    # Freeze the exact seven observations named by the historical analysis.
    # New audit snapshots must neither silently enter this comparison nor
    # force a rewrite of the original numeric/provenance report.
    historical_scope = tmp_path / "historical_observations"
    historical_scope.mkdir()
    for record in expected["observations"]:
        name = record["observation"]
        shutil.copytree(CASE / "chain_observations" / name, historical_scope / name)
    actual = analyze(historical_scope, CASE / "reference_variants",
                     CASE / "gamma_analysis/T_d0.01.npz", tmp_path / "analysis.json")
    assert actual["new_DFT_calls"] == 0 and not actual["physical_parameters_changed"]
    assert len(actual["observations"]) == len(expected["observations"]) >= 7
    names = [r["observation"] for r in actual["observations"]]
    assert names == sorted(names, key=str.casefold)
    for r, e in zip(actual["observations"], expected["observations"]):
        assert r["observation"] == e["observation"]
        assert r["replayed_fmax_eV_A"] == pytest.approx(e["replayed_fmax_eV_A"], abs=1e-10)
        assert r["provisional_discrete_forward_barrier_meV_fu"] - r["provisional_discrete_reverse_barrier_meV_fu"] == pytest.approx(r["reaction_energy_meV_fu"], abs=1e-10)
        assert r["all_images_inside_local_parent_chart"]
        assert all(i["T_Gamma_full_basis_reconstruction_residual_sqrt_amu_A"] < 1e-12 for i in r["images"])
    po_m = next(r for r in actual["observations"] if r["observation"] == "PO_M_lifted_step10")
    assert po_m["status"] == "unconverged_observation"
    assert po_m["residual_dominant_image"] == 5 and po_m["residual_dominant_block"] == "atomic"
    # A complete parent-pattern plane must not be claimed when the actual
    # projected peak has nearly half its displacement outside that subspace.
    peak = po_m["images"][po_m["highest_image_index"]]
    assert .50 < peak["parent_pattern_captured_squared_norm_fraction"] < .55
    assert peak["parent_pattern_residual_A"] > .5


@pytest.mark.parametrize("tamper", ["trajectory", "raw_energy", "snapshot", "pressure"])
def test_frozen_material_evidence_changes_fail_closed(tmp_path, tamper):
    folder = tmp_path / "observation"
    shutil.copytree(CASE / "chain_observations/PO_M_lifted_step10", folder)
    if tamper == "trajectory":
        with (folder / "evaluated_chain.traj").open("ab") as f:
            f.write(b"changed")
    elif tamper == "snapshot":
        with (folder / "POSCAR_03").open("a") as f:
            f.write("changed")
    else:
        source = folder / "observation.json"
        report = json.loads(source.read_text())
        if tamper == "raw_energy":
            report["raw_image_evaluations"][3]["energy_eV_cell"] += .001
        else:
            report["pressure_GPa"] = 1.
        source.write_text(json.dumps(report))
    t = read(CASE / "reference_variants/T.vasp")
    parent, patterns, _ = rotated_t_triplet(t)
    with np.load(CASE / "gamma_analysis/T_d0.01.npz") as gamma, pytest.raises(ValueError):
        analyze_observation(folder, t, parent, patterns, gamma)

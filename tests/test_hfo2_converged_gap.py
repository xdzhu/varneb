"""Replay the first G1 residual pass, without promoting it to a TS certificate."""

import json
from pathlib import Path

import numpy as np
import pytest
from ase.io import read

from examples.hfo2_fixed_input_factory import same_ordered_geometry
from scripts.analyze_hfo2_channel_network import read_evaluated_observation


ROOT = Path(__file__).resolve().parents[1] / "benchmarks/hfo2_channels/20261008"


def test_raw_audited_gap_chain_passes_ordinary_residual_and_preserves_endpoints():
    folder = ROOT / "converged_gap/gap_converged_step06"
    images, record, _ = read_evaluated_observation(folder)
    summary = json.loads((ROOT / "converged_gap/production_summary.json").read_text())
    assert len(images) == 10 and record["n_active_images"] == 8
    assert record["source_job_id"] == "28300425" and record["snapshot_step"] == 6
    assert summary["converged"] and summary["termination"] == "force_threshold"
    assert summary["fmax_target_eV_per_A"] == .10
    assert not summary["climbing_image_requested"]
    assert not summary["climbing_image_active_final"]
    assert record["replayed_fmax_eV_A"] == pytest.approx(.05988161569025288, abs=1e-12)
    assert same_ordered_geometry(images[0], read(ROOT / "reference_variants/T.vasp"))
    assert same_ordered_geometry(images[-1], read(ROOT / "reference_variants/PO.vasp"))
    energies = np.array([image.get_potential_energy() for image in images])
    assert (energies.max() - energies[0]) * 250 == pytest.approx(33.88897141439884)
    assert (energies.max() - energies[-1]) * 250 == pytest.approx(115.21016416145358)


def test_descriptive_analysis_retains_TS_and_sampling_limitations():
    analysis = json.loads((ROOT / "converged_gap/analysis.json").read_text())
    observation = analysis["observations"][0]
    assert observation["status"] == "NEB_residual_passed_TS_and_sampling_gates_pending"
    assert observation["highest_image_index"] == 3
    assert observation["residual_dominant_image"] == 2
    assert observation["residual_dominant_block"] == "atomic"
    assert observation["all_images_inside_local_parent_chart"]
    assert observation["reaction_energy_meV_fu"] == pytest.approx(-81.32119274705474)
    assert not analysis["physical_parameters_changed"] and analysis["new_DFT_calls"] == 0
    assert any("TS audits" in limitation for limitation in analysis["limitations"])

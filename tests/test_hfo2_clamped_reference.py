"""Portable native-data and declared-frame tests; no author mode arrays needed."""
import json
from pathlib import Path
import shutil
import subprocess

import numpy as np
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import read
import pytest

from scripts.analyze_hfo2_clamped_reference import (
    original_chart_images, pinned_json, terminal_images, REGISTRATION_SHA256,
)
from scripts.audit_hfo2_static_replica import sha256
from scripts.prepare_hfo2_clamped_endpoints import orient_long_axis_x
from vcneb.continuous_projection import continuous_reference_coordinates, project_reference_basis

ROOT = Path(__file__).resolve().parents[1]
CASE = ROOT/"benchmarks/hfo2_channels/20261008"
TERMINAL = CASE/"clamped_M_terminal_E068_20261011"
AUDIT_SHA = "56b92e8b15195d6f619a8402c7eee83a774b2e364e012478faa3fdc5d2da431c"


def orientation():
    return json.loads((CASE/"clamped_endpoint_seeds/clamped_seed_manifest.json").read_text())["orientation"]


def test_fixed_frame_roundtrip_preserves_unwrapped_order_and_sources():
    t = read(CASE/"reference_variants/T.vasp")
    original = [t.copy() for _ in range(3)]
    for i, image in enumerate(original):
        q = image.get_scaled_positions(wrap=False)
        q[1, 0] += 1
        q[2, 2] -= 2
        q[:, 1] += i*.01
        image.set_scaled_positions(q)
    oriented = [orient_long_axis_x(image) for image in original]
    for image in oriented:
        image.calc = SinglePointCalculator(image, energy=-1)
    before = [(a.positions.copy(), a.cell.array.copy(), a.calc) for a in oriented]
    restored = original_chart_images(oriented, orientation(), t)
    for result, expected, source, (positions, cell, calc) in zip(restored, original, oriented, before):
        np.testing.assert_allclose(result.cell.array, expected.cell.array, atol=1e-14, rtol=0)
        np.testing.assert_allclose(result.positions, expected.positions, atol=1e-14, rtol=0)
        np.testing.assert_allclose(result.get_scaled_positions(wrap=False),
                                   expected.get_scaled_positions(wrap=False), atol=1e-14, rtol=0)
        assert result.calc is None  # unrotated E/F/stress must not be attached
        assert np.array_equal(source.positions, positions) and np.array_equal(source.cell.array, cell)
        assert source.calc is calc


@pytest.mark.parametrize("field,value", [
    ("new_axes_from_old", [0, 1, 2]), ("cell_row_permutation", [0, 1, 2]),
    ("atom_permutation", list(reversed(range(12)))), ("determinant", -1),
    ("cartesian_rotation_rows", np.eye(3).tolist()),
])
def test_unregistered_frames_or_atom_rematching_are_rejected(field, value):
    t = read(CASE/"reference_variants/T.vasp")
    wrong = orientation()
    wrong[field] = value
    with pytest.raises(ValueError, match="registered fixed"):
        original_chart_images([orient_long_axis_x(t)], wrong, t)


def test_common_chart_projection_matches_declared_roundtrip_not_lab_axis_mix():
    t = read(CASE/"reference_variants/T.vasp")
    images = [t.copy() for _ in range(3)]
    rng = np.random.default_rng(20261011)
    for i, image in enumerate(images):
        image.positions += (i+1)*rng.normal(size=(12, 3))*.01
    basis = np.linalg.qr(rng.normal(size=(36, 4)))[0]
    true_chart = continuous_reference_coordinates(images, t)
    converted = original_chart_images([orient_long_axis_x(a) for a in images], orientation(), t)
    actual_chart = continuous_reference_coordinates(converted, t)
    a = project_reference_basis(true_chart["displacements_A"], basis)
    b = project_reference_basis(actual_chart["displacements_A"], basis)
    np.testing.assert_allclose(a["amplitudes"], b["amplitudes"], atol=1e-12, rtol=0)
    np.testing.assert_allclose(a["residual_norm"], b["residual_norm"], atol=1e-12, rtol=0)
    np.testing.assert_allclose(true_chart["green_strains"], actual_chart["green_strains"], atol=1e-12, rtol=0)


def test_real_clamped_terminal_is_reparsed_without_any_external_launch(monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("a zero-DFT analysis must not invoke a scheduler or calculator")
    monkeypatch.setattr(subprocess, "run", forbidden)
    images, report, residual, receipt = terminal_images(TERMINAL, AUDIT_SHA)
    assert len(images) == 9 and all(isinstance(a.calc, SinglePointCalculator) for a in images)
    assert residual["fmax_eV_A"] == pytest.approx(.0992140881476476, abs=1e-10, rel=0)
    assert report["mechanical_boundary"]["fixed_cell_rows"] == [0, 1]
    assert receipt["fresh_interior_SCFs"] == 91  # historical HF audit, not 91 new local reads


@pytest.mark.parametrize("digest", ["0"*64, "A"*64, "x", True, None])
def test_unpinned_or_invalid_evidence_digest_is_rejected(digest):
    with pytest.raises(ValueError, match="digest"):
        pinned_json(TERMINAL/"audit_receipt.json", digest)


@pytest.mark.parametrize("name", ["INPUT", "KPT", "STRU", "OUT.ABACUS/running_scf.log", "POSCAR"])
def test_portable_native_or_snapshot_tamper_is_not_accepted(tmp_path, name):
    copied = tmp_path/"terminal"
    shutil.copytree(TERMINAL, copied, ignore=shutil.ignore_patterns("__pycache__"))
    observation = copied/"observations/step_0013"
    path = observation/"POSCAR_03" if name == "POSCAR" else observation/"raw/image_0003"/name
    path.write_bytes(path.read_bytes()+b"\n# changed\n")
    with pytest.raises(ValueError, match="differs"):
        terminal_images(copied, AUDIT_SHA)


def test_registered_reference_report_keeps_all_frames_and_claim_boundaries():
    path = CASE/"clamped_reference_E071_20261011/analysis.json"
    result = json.loads(path.read_text())
    prior = pinned_json(CASE/"cmma_path_mapping/analysis.json", REGISTRATION_SHA256)
    assert [r["frame_id"] for r in result["Cmma_reference_frames"]] == [r["id"] for r in prior["selected_frame_registrations"]]
    assert result["fixed_G2_orientation"] == orientation()
    assert result["strain_training_condition"] == 0
    assert result["source_job_id"] == "28722320" and result["highest_image_index"] == 3
    assert result["new_DFT_calls"] == 0 and result["physical_parameters_changed"] is False
    assert result["complete_six_physical_bytes_rechecked_locally"] is False
    assert result["portable_native_frames_reparsed"] == 9
    assert result["T_triplet_rank3_fraction"][3] == pytest.approx(.5792373469688469, abs=1e-12, rel=0)
    assert result["T_triplet_rank3_fraction"][-1] == pytest.approx(.2443997899474907, abs=1e-12, rel=0)
    assert result["T_lowest_two_complete_doublets_rank4_indices"] == [3, 4, 5, 6]
    assert all(r["basis_audit"]["real_rank"] == 4 for r in result["Cmma_reference_frames"])
    assert all(r["basis_audit"]["ASR_or_translation_removal_from_basis"] is False for r in result["Cmma_reference_frames"])
    assert "not_prediction" in result["status"]
    assert sha256(TERMINAL/"audit_receipt.json") == result["terminal_audit_sha256"]


def test_executed_source_capsules_retain_their_actual_byte_pins():
    folder = CASE/"clamped_reference_E071_20261011"
    result = json.loads((folder/"analysis.json").read_text())
    rejected = json.loads((folder/"rejected_untransported_trial.json").read_text())
    assert sha256(folder/"executed_analyzer.py") == result["script_sha256"]
    assert sha256(folder/"executed_untransported_trial.py") == rejected["script_sha256"]


def test_rejected_frame_trial_is_not_an_energetic_or_production_change():
    folder = CASE/"clamped_reference_E071_20261011"
    result = json.loads((folder/"analysis.json").read_text())
    rejected = json.loads((folder/"rejected_untransported_trial.json").read_text())
    assert rejected["T_triplet_rank3_fraction"][3] < .02
    assert result["T_triplet_rank3_fraction"][3] > .57
    for key in ("terminal_audit_sha256", "trajectory_sha256", "relative_energy_meV_fu", "replayed_fmax_eV_A"):
        assert result[key] == rejected[key]
    assert "fixed_G2_orientation" not in rejected
    assert result["new_DFT_calls"] == rejected["new_DFT_calls"] == 0

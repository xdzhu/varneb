"""The published GaN joint-curvature evidence must remain reproducible."""

import hashlib
import json
from pathlib import Path

import numpy as np

from scripts.prepare_gan_joint_curvature import prepare
from scripts.prepare_gan_ts_bracketed_trial import (
    choose_guarded_trial, cross_axis_skew_angstrom,
)
from scripts.prepare_gan_equivalent_basis_probe import TRANSFORMS, physically_equivalent
from vcneb.joint_curvature import JointCurvatureCoordinates
from ase.io import read


ROOT = Path(__file__).resolve().parents[1]
EVIDENCE = ROOT / "outputs/gan_b4_b1_joint_curvature_20260926"
TRAJECTORY = ROOT / "outputs/gan_b4_b1_gamma_1x1x1_20260926/gan_vasp_tetragonal_final_29_images.traj"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_published_two_step_hessian_provenance_and_qualifier(tmp_path) -> None:
    summary = EVIDENCE / "source_summary.json"
    comparison_path = EVIDENCE / "comparison_0p01_0p02_v1.json"
    comparison = json.loads(comparison_path.read_text(encoding="utf-8"))
    assert "not_TS_certificate" in comparison["status"]
    assert comparison["negative_count_below_minus_0p1_eV_per_A2"] == [1, 1]
    assert comparison["lowest_mode_squared_overlap_between_steps"] > 0.99999
    assert comparison["relative_hessian_frobenius_difference"] < 0.002
    assert min(comparison["energy_gradient_max_abs_difference_eV_per_A"]) > 0.02
    for index, (label, distance) in enumerate((("0p01", 0.01), ("0p02", 0.02))):
        manifest_path = EVIDENCE / f"step_{label}_manifest.json"
        audit_path = EVIDENCE / f"step_{label}_audit.json"
        archive_path = EVIDENCE / f"step_{label}_joint_hessian.npz"
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        audit = json.loads(audit_path.read_text(encoding="utf-8"))
        assert manifest["status"] == "inputs_finalized_no_DFT"
        assert manifest["step_A"] == distance
        assert manifest["supercell_matrix"] == np.eye(3, dtype=int).tolist()
        assert len(manifest["cases"]) == audit["n_static_displacements"] == 36
        assert manifest["source_sha256"]["trajectory"] == _sha256(TRAJECTORY)
        assert manifest["source_sha256"]["summary"] == _sha256(summary)
        assert manifest["source_sha256"]["preparer"] == _sha256(
            ROOT / "scripts/prepare_gan_joint_curvature.py"
        )
        assert manifest["source_sha256"]["joint_curvature"] == _sha256(
            ROOT / "vcneb/joint_curvature.py"
        )
        assert audit["source_sha256"]["manifest"] == _sha256(manifest_path)
        assert audit["source_sha256"]["auditor"] == _sha256(
            ROOT / "scripts/audit_gan_joint_curvature.py"
        )
        assert audit["center_gradient_translation_free_eV_per_A"] > 0.05
        assert audit["translation_free_eigenvalues_eV_per_A2"][0] < -4.0
        assert audit["translation_free_eigenvalues_eV_per_A2"][1] > 2.0
        assert len(audit["cases"]) == 36
        assert all(len(item["outcar_sha256"]) == 64 for item in audit["cases"])
        with np.load(archive_path, allow_pickle=False) as arrays:
            assert arrays["hessian"].shape == (18, 18)
            assert arrays["eigenvectors"].shape == (18, 15)
            np.testing.assert_allclose(
                arrays["eigenvalues"], audit["translation_free_eigenvalues_eV_per_A2"],
            )
        assert comparison["source_sha256"]["audits"][index] == _sha256(audit_path)
        assert comparison["source_sha256"]["hessian_archives"][index] == _sha256(archive_path)

        # Recreate every four-atom displaced POSCAR from the public source
        # chain and compare the exact geometry hashes, without any DFT call.
        regenerated = prepare(TRAJECTORY, summary, tmp_path / label, distance)
        assert [case["POSCAR_sha256"] for case in regenerated["cases"]] == [
            case["POSCAR_sha256"] for case in manifest["cases"]
        ]


def test_published_800ev_probe_is_isolated_and_not_misread_as_convergence() -> None:
    manifest_path = EVIDENCE / "strain_encut800_manifest.json"
    report_path = EVIDENCE / "strain_encut800_diagnostic_v1.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    report = json.loads(report_path.read_text(encoding="utf-8"))
    assert manifest["purpose"] == "GaN_image15_VASP_strain_energy_stress_cutoff_sensitivity_only"
    assert manifest["encut_eV"] == 800 and manifest["baseline_encut_eV"] == 600
    assert len(manifest["cases"]) == len(report["cases"]) == 7
    assert manifest["source_600_manifest_sha256"] == _sha256(EVIDENCE / "step_0p02_manifest.json")
    assert manifest["preparer_sha256"] == _sha256(
        ROOT / "scripts/prepare_gan_strain_cutoff_diagnostic.py"
    )
    assert report["source_sha256"]["manifest_800"] == _sha256(manifest_path)
    assert report["source_sha256"]["audit_600"] == _sha256(EVIDENCE / "step_0p02_audit.json")
    assert report["source_sha256"]["archive_600"] == _sha256(
        EVIDENCE / "step_0p02_joint_hessian.npz"
    )
    assert report["source_sha256"]["auditor"] == _sha256(
        ROOT / "scripts/audit_gan_strain_cutoff_diagnostic.py"
    )
    assert report["mean_absolute_gradient_mismatch_eV_per_A"]["800"] < (
        report["mean_absolute_gradient_mismatch_eV_per_A"]["600"]
    )
    assert report["mean_absolute_gradient_mismatch_eV_per_A"]["800"] > 0.01
    assert all(len(item["outcar_sha256"]) == 64 for item in report["cases"])
    assert "not_production" in report["status"]


def test_published_three_cutoff_stress_series_remains_local_diagnostic() -> None:
    manifest_path = EVIDENCE / "strain_encut1000_manifest.json"
    series_path = EVIDENCE / "strain_encut1000_series_v2.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    series = json.loads(series_path.read_text(encoding="utf-8"))
    assert manifest["encut_eV"] == 1000 and manifest["baseline_encut_eV"] == 600
    assert len(manifest["cases"]) == len(series["cases"]) == 7
    assert manifest["source_600_manifest_sha256"] == _sha256(EVIDENCE / "step_0p02_manifest.json")
    assert manifest["preparer_sha256"] == _sha256(
        ROOT / "scripts/prepare_gan_strain_cutoff_diagnostic.py"
    )
    assert series["source_sha256"]["manifest_1000"] == _sha256(manifest_path)
    assert series["source_sha256"]["report_800"] == _sha256(
        EVIDENCE / "strain_encut800_diagnostic_v1.json"
    )
    assert series["source_sha256"]["auditor"] == _sha256(
        ROOT / "scripts/audit_gan_strain_cutoff_series.py"
    )
    mean = series["mean_absolute_gradient_mismatch_eV_per_A"]
    assert mean["600"] > mean["800"] > mean["1000"] > 0
    assert max(abs(value) for value in series["equivalent_stress_difference_kbar"]["1000"]) < 1
    gradients = series["center_gradient_translation_free_eV_per_A"]
    assert gradients["1000"] > gradients["600"] > 0.05
    assert "not_production_or_TS_certificate" in series["status"]


def test_published_newton_probe_preserves_failure_as_missing_result() -> None:
    manifest_path = EVIDENCE / "newton_probe1000_manifest.json"
    report_path = EVIDENCE / "newton_probe1000_partial_audit_v1.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    report = json.loads(report_path.read_text(encoding="utf-8"))
    assert manifest["evaluation_encut_eV"] == 1000
    assert manifest["source_hessian_encut_eV"] == 600
    assert manifest["status"] == "inputs_finalized_no_DFT"
    assert [case["name"] for case in manifest["cases"]] == ["newton_half", "newton_full"]
    assert manifest["source_sha256"]["trajectory"] == _sha256(TRAJECTORY)
    assert manifest["source_sha256"]["hessian_archive"] == _sha256(
        EVIDENCE / "step_0p02_joint_hessian.npz"
    )
    assert manifest["source_sha256"]["preparer"] == _sha256(
        ROOT / "scripts/prepare_gan_ts_newton_probe.py"
    )
    assert report["source_sha256"]["manifest"] == _sha256(manifest_path)
    assert report["source_sha256"]["auditor"] == _sha256(
        ROOT / "scripts/audit_gan_ts_newton_probe.py"
    )
    assert "not_TS_certificate" in report["status"]
    assert report["best_trial"] == "newton_half"
    half, full = report["cases"]
    assert half["status"] == "complete_converged_static"
    assert half["translation_free_gradient_eV_per_A"] < 0.05
    assert 0.4 < half["gradient_ratio_to_center"] < 0.6
    assert len(half["outcar_sha256"]) == 64
    assert "failure_before_SCF" in full["status"]
    assert full["energy_eV_per_cell"] is None
    assert full["translation_free_gradient_eV_per_A"] is None
    assert full["gradient_ratio_to_center"] is None
    assert len(full["outcar_sha256"]) == len(full["stdout_sha256"]) == 64
    for record, staged in zip(report["cases"], manifest["cases"]):
        assert record["input_sha256"] == staged["input_sha256"]


def test_bracketed_trial_preflight_rejects_risky_three_quarter_step() -> None:
    manifest = json.loads((EVIDENCE / "newton_probe1000_manifest.json").read_text())
    reference = read(TRAJECTORY, index=15)
    coordinates = JointCurvatureCoordinates(reference, manifest["cell_scale_A"])
    correction = np.asarray(manifest["correction_A"], dtype=float)
    skews = {
        "center": cross_axis_skew_angstrom(reference.cell.array),
        "completed_half": cross_axis_skew_angstrom(
            coordinates.displaced(0.5 * correction).cell.array
        ),
        "failed_full": cross_axis_skew_angstrom(
            coordinates.displaced(correction).cell.array
        ),
    }
    fraction, trial, skew, limit, options = choose_guarded_trial(
        coordinates, correction, skews,
    )
    assert fraction == 0.625
    assert options[0]["fraction"] == 0.75
    assert not options[0]["passes_empirical_margin"]
    assert options[1]["passes_empirical_margin"]
    assert skew < limit
    assert trial.get_volume() > 0
    assert trial.get_chemical_symbols() == reference.get_chemical_symbols()


def test_equivalent_basis_probes_preserve_full_step_physics_and_kmesh() -> None:
    manifest = json.loads((EVIDENCE / "newton_probe1000_manifest.json").read_text())
    reference = read(TRAJECTORY, index=15)
    coordinates = JointCurvatureCoordinates(reference, manifest["cell_scale_A"])
    original = coordinates.displaced(np.asarray(manifest["correction_A"], dtype=float))
    for matrix in TRANSFORMS.values():
        candidate = original.copy()
        candidate.set_cell(np.asarray(matrix) @ original.cell.array, scale_atoms=False)
        assert physically_equivalent(original, candidate, np.asarray(matrix))
    unsuitable = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]])
    candidate = original.copy()
    candidate.set_cell(unsuitable @ original.cell.array, scale_atoms=False)
    assert not physically_equivalent(original, candidate, unsuitable)


def test_published_bracketed_trial_has_only_audited_gradient_descent() -> None:
    manifest_path = EVIDENCE / "bracketed_0p625_manifest.json"
    audit_path = EVIDENCE / "bracketed_0p625_audit_v1.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    assert manifest["evaluation_encut_eV"] == 1000
    assert manifest["status"] == "inputs_finalized_no_DFT"
    assert manifest["cases"][0]["fraction_from_original_center"] == 0.625
    assert not manifest["empirical_bravais_risk"]["candidate_options"][0][
        "passes_empirical_margin"
    ]
    assert manifest["cases"][0]["cross_axis_skew_A"] < manifest[
        "empirical_bravais_risk"
    ]["allowed_skew_A"]
    assert manifest["source_sha256"]["preparer"] == _sha256(
        ROOT / "scripts/prepare_gan_ts_bracketed_trial.py"
    )
    assert manifest["source_sha256"]["previous_manifest"] == _sha256(
        EVIDENCE / "newton_probe1000_manifest.json"
    )
    assert audit["source_sha256"]["manifest"] == _sha256(manifest_path)
    assert audit["source_sha256"]["auditor"] == _sha256(
        ROOT / "scripts/audit_gan_ts_bracketed_trial.py"
    )
    assert audit["result"]["status"] == "complete_converged_static_not_TS_certificate"
    assert 0.03 < audit["result"]["translation_free_gradient_eV_per_A"] < 0.04
    assert 0.7 < audit["result"]["gradient_ratio_to_previous_half"] < 0.8
    assert audit["result"]["input_sha256"] == manifest["cases"][0]["input_sha256"]
    assert len(audit["result"]["outcar_sha256"]) == 64


def test_published_equivalent_basis_probe_does_not_count_unrun_case_as_failure() -> None:
    manifest_path = EVIDENCE / "equivalent_basis_full_manifest.json"
    audit_path = EVIDENCE / "equivalent_basis_full_two_trials_audit_v1.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    assert manifest["evaluation_encut_eV"] == 1000
    assert manifest["kmesh"] == [8, 8, 6]
    assert manifest["source_sha256"]["preparer"] == _sha256(
        ROOT / "scripts/prepare_gan_equivalent_basis_probe.py"
    )
    assert audit["source_sha256"]["manifest"] == _sha256(manifest_path)
    assert audit["source_sha256"]["auditor"] == _sha256(
        ROOT / "scripts/audit_gan_equivalent_basis_probe.py"
    )
    assert [result["status"] for result in audit["results"]] == [
        "VASP_6p3p2_Bravais_classification_failure_before_SCF",
        "not_evaluated",
        "VASP_6p3p2_Bravais_classification_failure_before_SCF",
    ]
    for staged, result in zip(manifest["cases"], audit["results"]):
        assert staged["input_sha256"] == result["input_sha256"]
    assert all(result["energy_eV_per_cell"] is None for result in
               (audit["results"][0], audit["results"][2]))

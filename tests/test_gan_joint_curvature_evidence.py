"""The published GaN joint-curvature evidence must remain reproducible."""

import hashlib
import json
from pathlib import Path

import numpy as np

from scripts.prepare_gan_joint_curvature import prepare


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

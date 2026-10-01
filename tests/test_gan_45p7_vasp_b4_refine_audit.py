"""A native endpoint stop is not automatically a 2-kbar stress pass."""

import pytest

from scripts import audit_gan_45p7_vasp_b4_refine as refinement


SOURCE = refinement.Path(
    "paper/VARNEB_CPC/evidence/gan_vasp_b4_endpoint_refine_20261001"
)
STRICT_SOURCE = refinement.Path(
    "paper/VARNEB_CPC/evidence/gan_vasp_b4_endpoint_strict_20261001"
)


def test_first_native_b4_refinement_stopped_without_moving():
    report = refinement.audit(SOURCE, slurm_job=27812043)
    assert report["n_evaluated_ionic_frames"] == 1
    assert report["vasp_reports_ionic_convergence"]
    assert report["final_max_atomic_force_eV_per_A"] == pytest.approx(0.005532, abs=1e-6)
    assert report["final_max_raw_stress_residual_kbar"] == pytest.approx(2.91154)
    assert not report["stress_pass_2kbar"]
    assert abs(report["final_minus_original_H_meV_per_GaN"]) < 1e-6
    assert report["start_to_final_cell_max_abs_A"] < 1e-8
    assert report["start_to_final_position_max_abs_A"] < 1e-5


def test_strict_same_electronic_contract_refinement_passes_pressure_gate():
    report = refinement.audit(STRICT_SOURCE, slurm_job=27812112)
    assert report["n_evaluated_ionic_frames"] == 5
    assert report["vasp_reports_ionic_convergence"]
    assert report["force_pass_0p02eV_per_A"]
    assert report["stress_pass_2kbar"]
    assert report["final_max_raw_stress_residual_kbar"] == pytest.approx(0.05888)
    assert report["final_minus_original_H_meV_per_GaN"] == pytest.approx(-0.0359775)
    assert report["final_vs_original_B4"]["relative_displacement_max_A"] < 0.003
    assert report["final_vs_original_B4"]["relative_volume_difference"] < 0.0001


def test_refinement_audit_rejects_input_manifest_drift(tmp_path):
    manifest = (SOURCE / "manifest.json").read_text(encoding="utf-8")
    for name in ("INCAR", "KPOINTS", "POSCAR", "OUTCAR", "CONTCAR",
                 "sha256.inputs.json"):
        (tmp_path / name).write_bytes((SOURCE / name).read_bytes())
    (tmp_path / "manifest.json").write_text(
        manifest.replace('"ENCUT_eV": 600', '"ENCUT_eV": 700'),
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="not a pinned 600-eV"):
        refinement.audit(tmp_path)

"""Published GaN Gamma evidence must remain tied to its 1x1x1 DFT audit."""

from __future__ import annotations

import hashlib
import csv
import json
from pathlib import Path

import numpy as np
from ase.io import read
from ase.units import GPa

from examples.analyze_path_gamma_modes import make_report
from vcneb.phonons import load_gamma_force_constants, load_phonopy_gamma_eigenpairs


ROOT = Path(__file__).resolve().parents[1]
EVIDENCE = ROOT / "outputs" / "gan_b4_b1_gamma_1x1x1_20260926"


def test_endpoint_gamma_evidence_contract_and_numerical_stability() -> None:
    manifest = json.loads((EVIDENCE / "manifest.json").read_text(encoding="utf-8"))
    audit = json.loads((EVIDENCE / "audit.json").read_text(encoding="utf-8"))
    subspaces = json.loads((EVIDENCE / "subspace_audit_v2.json").read_text(encoding="utf-8"))
    assert manifest["supercell_matrix"] == np.eye(3, dtype=int).tolist()
    assert len(manifest["cases"]) == 32
    assert len({item["name"] for item in manifest["cases"]}) == 32
    assert audit["source_manifest_sha256"] == hashlib.sha256(
        (EVIDENCE / "manifest.json").read_bytes()
    ).hexdigest()
    assert audit["status"] == "computed_GaN_endpoint_atomic_Gamma_1x1x1_requires_numerical_interpretation"
    assert subspaces["status"] == "audited_GaN_endpoint_Gamma_atomic_subspaces_not_TS_modes"
    assert len(audit["source_code_sha256"]) == 3
    assert all(len(value) == 64 for value in audit["source_code_sha256"].values())
    for role, source in {
        "preparer": ROOT / "scripts/prepare_gan_gamma_phonopy_1x1x1.py",
        "launcher": ROOT / "cluster/hf_gan_gamma_phonopy_1x1x1.slurm",
        "auditor": ROOT / "scripts/audit_gan_gamma_phonopy_1x1x1.py",
    }.items():
        assert hashlib.sha256(source.read_bytes()).hexdigest() == audit["source_code_sha256"][role]
    for phase in ("B4", "B1"):
        endpoint = read(EVIDENCE / f"{phase}_CONTCAR", format="vasp")
        assert endpoint.get_chemical_symbols() == ["Ga", "Ga", "N", "N"]
        assert endpoint.get_volume() > 0
        comparison = audit["step_size_comparison"][phase]
        assert comparison["relative_fc_frobenius_difference"] < 0.002
        assert comparison["max_optical_frequency_difference_thz"] < 0.02
        for distance in ("0.01", "0.02"):
            label = f"{phase}_d{distance}"
            family = audit["families"][label]
            assert family["n_displacements"] == 8
            assert len(family["cases"]) == 8
            assert family["max_abs_acoustic_frequency_thz"] < 0.05
            assert family["max_force_constant_asr_drift_eV_per_A2"] < 0.001
            assert min(family["frequencies_thz"][3:]) > 0
            assert all(len(case["outcar_sha256"]) == 64 for case in family["cases"])
            with np.load(EVIDENCE / f"{label}.npz", allow_pickle=False) as data:
                assert data["force_constants"].shape == (4, 3, 4, 3)
                assert data["eigenvectors_mass_weighted"].shape == (12, 12)
                np.testing.assert_allclose(data["frequencies_thz"], family["frequencies_thz"])
        phase_analysis = subspaces["phases"][phase]
        assert phase_analysis["three_group_captured_atomic_path_squared_norm_fraction"] > 0.99999999
        assert len(phase_analysis["groups_ranked_by_path_amplitude"]) == 3
        for group in phase_analysis["groups_ranked_by_path_amplitude"]:
            assert group["minimum_principal_subspace_overlap"] > 0.999
            assert group["max_Q_norm_step_size_difference_sqrt_amu_A"] < 0.001
    for filename, expected in subspaces["source_sha256"].items():
        path = EVIDENCE / filename
        if not path.exists():
            path = ROOT / "benchmarks" / "numerical_integrity" / filename
        assert path.exists()
        assert hashlib.sha256(path.read_bytes()).hexdigest() == expected


def test_public_final_chain_reproduces_endpoint_gamma_projections() -> None:
    chain = read(EVIDENCE / "gan_vasp_tetragonal_final_29_images.traj", index=":")
    assert len(chain) == 29
    assert all(image.get_chemical_symbols() == ["Ga", "Ga", "N", "N"] for image in chain)
    with (ROOT / "benchmarks/numerical_integrity/gan_b4_b1_tetragonal_empirical_atomic_strain_20260926.csv").open(
        newline="", encoding="utf-8"
    ) as handle:
        rows = list(csv.DictReader(handle))
    enthalpy = np.array([image.get_potential_energy() + 45.7 * GPa * image.get_volume()
                         for image in chain])
    np.testing.assert_allclose(
        (enthalpy - enthalpy[0]) / 2,
        [float(row["relative_enthalpy_eV_per_GaN"]) for row in rows],
        rtol=0, atol=1e-7,
    )
    for phase in ("B4", "B1"):
        reference = read(EVIDENCE / f"{phase}_CONTCAR", format="vasp")
        for distance in ("0.01", "0.02"):
            label = f"{phase}_d{distance}"
            archive = EVIDENCE / f"{label}.npz"
            fc, masses = load_gamma_force_constants(archive)
            modes = load_phonopy_gamma_eigenpairs(archive)
            reproduced, _ = make_report(
                reference=reference, images=chain, force_constants=fc,
                masses_amu=masses, phonopy_eigenpairs=modes,
                include_translations=False,
            )
            stored = json.loads((EVIDENCE / f"{label}_path_projection.json").read_text(encoding="utf-8"))
            np.testing.assert_allclose(
                reproduced["normal_coordinates_sqrt_amu_A"],
                stored["normal_coordinates_sqrt_amu_A"], rtol=0, atol=1e-10,
            )

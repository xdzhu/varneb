"""Published GaN Gamma evidence must remain tied to its 1x1x1 DFT audit."""

from __future__ import annotations

import hashlib
import csv
import json
from pathlib import Path

import numpy as np
import pytest
from ase.io import read
from ase.units import GPa

from examples.analyze_path_gamma_modes import make_report
from scripts.analyze_gan_gamma_path_subspaces import _dominant_groups, audit as audit_subspaces
from vcneb.phonons import (
    identify_acoustic_modes,
    load_gamma_force_constants,
    load_phonopy_gamma_eigenpairs,
)


ROOT = Path(__file__).resolve().parents[1]
EVIDENCE = ROOT / "outputs" / "gan_b4_b1_gamma_1x1x1_20260926"


def test_soft_optical_group_is_not_discarded_by_acoustic_frequency_cutoff() -> None:
    coordinates = np.zeros((29, 12))
    coordinates[:, 0] = 3.0  # genuine soft optical mode ordered before acoustics
    coordinates[:, 4] = 2.0
    coordinates[:, 5] = 1.0
    groups = _dominant_groups({"degenerate_mode_subspaces": [
        {"mode_indices": [0], "frequency_cm1": -2.0},
        {"mode_indices": [1, 2, 3], "frequency_cm1": 0.2},
        {"mode_indices": [4], "frequency_cm1": 100.0},
        {"mode_indices": [5], "frequency_cm1": 200.0},
    ]}, coordinates, {1, 2, 3})
    assert [group["mode_indices"] for group in groups] == [[0], [4], [5]]
    with pytest.raises(ValueError, match="mixed acoustic/optical"):
        _dominant_groups({"degenerate_mode_subspaces": [
            {"mode_indices": [0, 1], "frequency_cm1": -2.0},
        ]}, coordinates, {1, 2, 3})


def test_geometric_acoustic_identification_preserves_archived_gan_projections() -> None:
    source = ROOT / "benchmarks/numerical_integrity/gan_b4_b1_tetragonal_empirical_atomic_strain_20260926.csv"
    reproduced, rows = audit_subspaces(EVIDENCE, source)
    archived = json.loads((EVIDENCE / "subspace_audit_v2.json").read_text(encoding="utf-8"))
    assert reproduced == archived
    assert len(rows) == 58


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
            assert reproduced["interpretation"]["acoustic_mode_indices"] == [0, 1, 2]
            assert reproduced["interpretation"]["minimum_rigid_translation_subspace_overlap"] > 0.999


def test_gan_three_group_residual_separates_optical_and_acoustic_components() -> None:
    subspaces = json.loads((EVIDENCE / "subspace_audit_v2.json").read_text(encoding="utf-8"))
    figure_csv = ROOT / "paper/VARNEB_CPC/figures/gan_gamma_path_600eV_source_data.csv"
    with figure_csv.open(newline="", encoding="utf-8") as handle:
        plotted = list(csv.DictReader(handle))
    assert len(plotted) == 29
    expected = {"B4": (0.0002362187675, 0.0003597366448),
                "B1": (0.0001826949300, 0.0013915895328)}
    for phase in ("B4", "B1"):
        with np.load(EVIDENCE / f"{phase}_d0.01.npz", allow_pickle=False) as data:
            acoustic, minimum = identify_acoustic_modes(
                data["eigenvectors_mass_weighted"], data["masses_amu"]
            )
        assert minimum > 0.999
        coordinates = np.asarray(json.loads(
            (EVIDENCE / f"{phase}_d0.01_path_projection.json").read_text(encoding="utf-8")
        )["normal_coordinates_sqrt_amu_A"])
        selected = {index for group in subspaces["phases"][phase]["groups_ranked_by_path_amplitude"]
                    for index in group["mode_indices"]}
        acoustic_set = set(acoustic)
        optical_remainder = [i for i in range(12) if i not in selected | acoustic_set]
        total_remainder = [i for i in range(12) if i not in selected]
        optical = np.linalg.norm(coordinates[:, optical_remainder], axis=1)
        leakage = np.linalg.norm(coordinates[:, acoustic], axis=1)
        total = np.linalg.norm(coordinates[:, total_remainder], axis=1)
        np.testing.assert_allclose(total**2, optical**2 + leakage**2, atol=1e-12)
        np.testing.assert_allclose(total,
            [float(row[f"{phase}_three_group_residual_sqrt_amu_A"]) for row in plotted],
            rtol=0, atol=1e-12)
        assert np.max(optical) == pytest.approx(expected[phase][0], abs=1e-9)
        assert np.max(leakage) == pytest.approx(expected[phase][1], abs=1e-9)
    manuscript = (ROOT / "paper/VARNEB_CPC/varneb_CPC.tex").read_text(encoding="utf-8")
    assert "0.000236/0.000183" in manuscript
    assert "0.000360/0.001392" in manuscript

"""No-DFT regression test for the Gamma path-projection CLI."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.io.trajectory import Trajectory


ROOT = Path(__file__).resolve().parents[1]
DRIVER = ROOT / "examples" / "analyze_path_gamma_modes.py"


def test_gamma_mode_driver_projects_a_complete_path_without_dft(tmp_path: Path) -> None:
    reference = Atoms("H", cell=[4, 4, 4], pbc=True)
    reference_path = tmp_path / "reference.traj"
    Trajectory(reference_path, "w").write(reference)
    path = tmp_path / "vcneb.traj"
    with Trajectory(path, "w") as trajectory:
        for fraction in (0.0, 0.05, 0.10):
            image = reference.copy()
            image.set_scaled_positions([[fraction, 0.0, 0.0]])
            trajectory.write(image)
    archive = tmp_path / "force_constants.npz"
    np.savez_compressed(archive, force_constants=np.diag([1.0, 4.0, 9.0]), masses_amu=reference.get_masses())
    output = tmp_path / "modes.json"
    result = subprocess.run(
        [
            sys.executable,
            str(DRIVER),
            "--reference",
            str(reference_path),
            "--trajectory",
            str(path),
            "--force-constants",
            str(archive),
            "--n-images",
            "3",
            "--output",
            str(output),
        ],
        text=True,
        capture_output=True,
        check=True,
    )
    assert "projected 3 images" in result.stdout
    report = json.loads(output.read_text(encoding="utf-8"))
    assert report["n_images"] == 3
    assert report["interpretation"]["cell_degrees_of_freedom"].startswith("excluded")
    assert len(report["normal_coordinates_sqrt_amu_A"]) == 3
    assert output.with_suffix(".npz").exists()


def test_reference_permutation_reorders_both_force_constant_axes() -> None:
    module = __import__("runpy").run_path(str(DRIVER))
    reference = Atoms("HHe", cell=[4, 4, 4], pbc=True)
    force_constants = np.arange(36.0).reshape(6, 6)
    reordered, reordered_fc, masses, permutation = module["reorder_reference_force_constants"](
        reference,
        force_constants,
        reference.get_masses(),
        "1,0",
    )
    assert permutation == [1, 0]
    assert reordered.get_chemical_symbols() == ["He", "H"]
    assert np.array_equal(masses, reference.get_masses()[[1, 0]])
    assert np.array_equal(reordered_fc, force_constants[np.ix_([3, 4, 5, 0, 1, 2], [3, 4, 5, 0, 1, 2])])


def test_gamma_mode_report_accepts_standard_atomic_weight_rounding() -> None:
    module = __import__("runpy").run_path(str(DRIVER))
    reference = Atoms("O", cell=[4, 4, 4], pbc=True)
    images = [reference.copy(), reference.copy()]
    report, _ = module["make_report"](
        reference=reference,
        images=images,
        force_constants=np.eye(3),
        masses_amu=np.array([15.9994]),
        include_translations=False,
    )
    assert report["n_atoms"] == 1


def test_gamma_mode_report_records_explicit_force_constant_unit() -> None:
    module = __import__("runpy").run_path(str(DRIVER))
    reference = Atoms("H", cell=[4, 4, 4], pbc=True)
    report, _ = module["make_report"](
        reference=reference,
        images=[reference.copy(), reference.copy()],
        force_constants=np.eye(3),
        masses_amu=reference.get_masses(),
        include_translations=False,
        force_constant_unit="eV/angstrom.au",
    )
    assert report["interpretation"]["force_constant_unit_input"] == "eV/angstrom.au"
    assert report["interpretation"]["force_constant_unit_handling"].startswith("converted")


def test_direct_phonopy_basis_is_used_after_path_atom_permutation() -> None:
    module = __import__("runpy").run_path(str(DRIVER))
    from vcneb.phonons import PhonopyGammaEigenpairs

    reference = Atoms("HHe", cell=[4, 4, 4], pbc=True)
    raw = PhonopyGammaEigenpairs(np.array([1.0] * 6), np.eye(6))
    reordered = module["reorder_phonopy_eigenpairs"](raw, [1, 0])
    assert np.array_equal(reordered.eigenvectors_mass_weighted, np.eye(6)[[3, 4, 5, 0, 1, 2]])
    report, _ = module["make_report"](
        reference=reference[[1, 0]],
        images=[reference[[1, 0]], reference[[1, 0]]],
        force_constants=np.eye(6),
        masses_amu=reference.get_masses()[[1, 0]],
        include_translations=True,
        phonopy_eigenpairs=reordered,
    )
    assert report["interpretation"]["mode_basis"] == "direct Phonopy Gamma eigenpairs"
    assert not report["interpretation"]["translations_projected_before_diagonalization"]


def test_reference_translation_uses_fractional_gauge_without_changing_cell() -> None:
    module = __import__("runpy").run_path(str(DRIVER))
    reference = Atoms("H", scaled_positions=[[0.1, 0.2, 0.3]], cell=[4, 5, 6], pbc=True)
    shifted, translation = module["translate_reference_fractional"](reference, "0.5,-0.25,0.125")
    assert translation == [0.5, -0.25, 0.125]
    assert np.allclose(shifted.cell.array, reference.cell.array)
    assert np.allclose(shifted.get_scaled_positions(wrap=False), [[0.6, -0.05, 0.425]])


def test_path_projection_removes_rigid_translation_and_groups_degenerate_modes() -> None:
    module = __import__("runpy").run_path(str(DRIVER))
    reference = Atoms("HHe", cell=[4, 4, 4], pbc=True)
    images = []
    for shift in (0.0, 0.2):
        image = reference.copy()
        image.set_positions(image.get_positions() + [shift, 0.0, 0.0])
        images.append(image)
    report, _ = module["make_report"](
        reference=reference,
        images=images,
        force_constants=np.eye(6),
        masses_amu=reference.get_masses(),
        include_translations=False,
    )
    assert np.allclose(report["normal_coordinates_sqrt_amu_A"], 0.0)
    assert report["interpretation"]["rigid_translations_removed_from_path"]
    assert sum(group["degeneracy"] for group in report["degenerate_mode_subspaces"]) == 6

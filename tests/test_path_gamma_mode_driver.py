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

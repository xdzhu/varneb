"""Test VASP Bravais initialization using exactly equivalent GaN cell bases.

Only unimodular signed permutations of the two equal-sampled basal vectors
are allowed. The Cartesian structure and Gamma-centered 8x8x6 k-point set
are unchanged; this is an input-representation diagnostic, not relaxation.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil

import numpy as np
from ase.io import read, write

from scripts.prepare_gan_ts_newton_probe import sha256


TRANSFORMS = {
    "basal_rotate_plus": [[0, 1, 0], [-1, 0, 0], [0, 0, 1]],
    "basal_rotate_minus": [[0, -1, 0], [1, 0, 0], [0, 0, 1]],
    "basal_invert": [[-1, 0, 0], [0, -1, 0], [0, 0, 1]],
}


def physically_equivalent(original, candidate, matrix: np.ndarray) -> bool:
    """Check lattice/atom equivalence and the unchanged physical k mesh."""

    integer = np.asarray(matrix, dtype=int)
    divisions = np.diag([8.0, 8.0, 6.0])
    if (integer.shape != (3, 3) or round(np.linalg.det(integer)) != 1
            or not np.allclose(integer @ integer.T, np.eye(3), atol=0, rtol=0)
            or not np.allclose(integer @ divisions, divisions @ integer, atol=0, rtol=0)
            or original.get_chemical_symbols() != candidate.get_chemical_symbols()
            or not np.allclose(candidate.cell.array, integer @ original.cell.array,
                               rtol=0, atol=1e-10)
            or not np.isclose(candidate.get_volume(), original.get_volume(),
                              rtol=0, atol=1e-9)):
        return False
    delta = (candidate.positions - original.positions) @ np.linalg.inv(original.cell.array)
    delta -= np.rint(delta)
    return bool(np.allclose(delta, 0, atol=1e-10, rtol=0))


def prepare(previous_work: Path, previous_manifest: Path, previous_audit: Path,
            work_root: Path) -> dict:
    if work_root.exists():
        raise FileExistsError(work_root)
    manifest = json.loads(previous_manifest.read_text(encoding="utf-8"))
    report = json.loads(previous_audit.read_text(encoding="utf-8"))
    previous_case = previous_work / "cases/newton_full"
    if (sha256(previous_work / "manifest.json") != sha256(previous_manifest)
            or report["source_sha256"]["manifest"] != sha256(previous_manifest)
            or report["cases"][1]["outcar_sha256"] != sha256(previous_case / "OUTCAR")
            or "failure_before_SCF" not in report["cases"][1]["status"]
            or manifest["evaluation_encut_eV"] != 1000):
        raise ValueError("failed full-step provenance is missing or changed")
    previous_input = manifest["cases"][1]["input_sha256"]
    for filename, expected in previous_input.items():
        if sha256(previous_case / filename) != expected:
            raise ValueError(f"source full-step input changed: {filename}")
    incar = (previous_case / "INCAR").read_text(encoding="utf-8")
    kpoints = (previous_case / "KPOINTS").read_text(encoding="utf-8")
    if ("ENCUT = 1000.000000" not in incar or "ISYM = -1" not in incar
            or "SYMPREC = 1.00e-04" not in incar or "Gamma\n8 8 6" not in kpoints):
        raise ValueError("the fixed VASP input contract changed")
    original = read(previous_case / "POSCAR", format="vasp")
    records = []
    for name, data in TRANSFORMS.items():
        transform = np.asarray(data, dtype=int)
        candidate = original.copy()
        candidate.set_cell(transform @ original.cell.array, scale_atoms=False)
        if not physically_equivalent(original, candidate, transform):
            raise ValueError(f"basis transformation is not physically equivalent: {name}")
        directory = work_root / "cases" / name
        directory.mkdir(parents=True)
        write(directory / "POSCAR", candidate, format="vasp", direct=True, sort=False)
        if not physically_equivalent(original, read(directory / "POSCAR", format="vasp"),
                                     transform):
            raise ValueError(f"POSCAR serialization changed the physical geometry: {name}")
        for filename in ("INCAR", "KPOINTS", "POTCAR"):
            shutil.copy2(previous_case / filename, directory / filename)
        hashes = {filename: sha256(directory / filename)
                  for filename in ("POSCAR", "INCAR", "KPOINTS", "POTCAR")}
        (directory / "sha256.inputs.json").write_text(
            json.dumps(hashes, indent=2) + "\n", encoding="utf-8",
        )
        records.append({"name": name, "unimodular_matrix": data,
                        "POSCAR_sha256": hashes["POSCAR"], "input_sha256": hashes,
                        "volume_A3": float(candidate.get_volume())})
    result = {
        "purpose": "GaN_full_step_equivalent_basis_VASP_Bravais_diagnostic_only",
        "status": "inputs_finalized_no_DFT",
        "evaluation_encut_eV": 1000,
        "pressure_GPa": manifest["pressure_GPa"],
        "kmesh": [8, 8, 6],
        "original_full_step_POSCAR_sha256": previous_input["POSCAR"],
        "cases": records,
        "source_sha256": {
            "previous_manifest": sha256(previous_manifest),
            "previous_audit": sha256(previous_audit),
            "full_failed_OUTCAR": sha256(previous_case / "OUTCAR"),
            "preparer": sha256(Path(__file__)),
        },
        "boundary": "No lattice/atom/k-mesh physics changes; VASP classification may still fail.",
    }
    (work_root / "manifest.json").write_text(json.dumps(result, indent=2) + "\n",
                                              encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--previous-work", type=Path, required=True)
    parser.add_argument("--previous-manifest", type=Path, required=True)
    parser.add_argument("--previous-audit", type=Path, required=True)
    parser.add_argument("--work-root", type=Path, required=True)
    args = parser.parse_args()
    result = prepare(args.previous_work, args.previous_manifest, args.previous_audit,
                     args.work_root)
    print(json.dumps({"status": result["status"],
                      "cases": [case["name"] for case in result["cases"]]}))


if __name__ == "__main__":
    main()

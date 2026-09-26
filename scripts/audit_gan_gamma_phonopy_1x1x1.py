"""Audit completed GaN 1x1x1 VASP displacements and build Gamma eigenpairs.

No DFT is launched. Both endpoint/step-size families must be complete before
the output directory is created. Fixed-cell Gamma modes are not a variable-
cell transition-state Hessian and omit non-analytic polar LO-TO corrections.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from ase.io import read
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def phonopy_for_reference(reference, distance_A: float) -> Phonopy:
    unit = PhonopyAtoms(
        symbols=reference.get_chemical_symbols(),
        cell=reference.cell.array,
        scaled_positions=reference.get_scaled_positions(wrap=False),
    )
    phonon = Phonopy(unit, supercell_matrix=np.eye(3, dtype=int),
                     primitive_matrix=np.eye(3), symprec=1e-4, is_symmetry=True)
    phonon.generate_displacements(
        distance=distance_A, is_plusminus=True, is_diagonal=False,
    )
    return phonon


def real_gamma_eigenvectors(vectors: np.ndarray, *, tolerance: float = 1e-8) -> np.ndarray:
    """Fix each arbitrary complex phase at Gamma and reject material imaginary parts."""

    values = np.asarray(vectors, dtype=complex).copy()
    if values.ndim != 2 or values.shape[0] != values.shape[1]:
        raise ValueError("Phonopy Gamma eigenvectors must form a square matrix")
    for column in range(values.shape[1]):
        pivot = np.argmax(np.abs(values[:, column]))
        values[:, column] *= np.exp(-1j * np.angle(values[pivot, column]))
    if np.max(np.abs(values.imag)) > tolerance:
        raise ValueError("Gamma eigenvectors remain materially complex after phase fixing")
    real = values.real
    if not np.allclose(real.T @ real, np.eye(len(real)), rtol=1e-8, atol=1e-8):
        raise ValueError("Gamma eigenvectors are not orthonormal")
    return real


def completed_static(case_dir: Path, expected_hashes: dict) -> dict:
    for name, expected in expected_hashes.items():
        if not (case_dir / name).is_file() or sha256(case_dir / name) != expected:
            raise ValueError(f"VASP input changed from staged manifest: {case_dir / name}")
    outcar = case_dir / "OUTCAR"
    if not outcar.is_file():
        raise FileNotFoundError(f"VASP displacement not completed: {outcar}")
    text = outcar.read_text(encoding="utf-8", errors="replace")
    if ("General timing and accounting informations" not in text
            or "aborting loop because EDIFF is reached" not in text):
        raise ValueError(f"VASP static lacks completed SCF/OUTCAR markers: {outcar}")
    atoms = read(outcar)
    supplied = read(case_dir / "POSCAR", format="vasp")
    delta = atoms.get_scaled_positions(wrap=False) - supplied.get_scaled_positions(wrap=False)
    delta -= np.rint(delta)
    if (atoms.get_chemical_symbols() != supplied.get_chemical_symbols()
            or not np.allclose(atoms.cell.array, supplied.cell.array, atol=2e-5)
            or np.max(np.abs(delta)) > 2e-5):
        raise ValueError(f"VASP output geometry differs from staged static POSCAR: {case_dir}")
    forces = np.asarray(atoms.get_forces(), dtype=float)
    energy = float(atoms.get_potential_energy())
    if forces.shape != (4, 3) or not np.isfinite(forces).all() or not np.isfinite(energy):
        raise ValueError(f"VASP force/energy is invalid: {case_dir}")
    return {
        "forces": forces, "energy_eV_per_cell": energy,
        "maximum_force_eV_per_A": float(np.max(np.linalg.norm(forces, axis=1))),
        "net_force_eV_per_A": np.sum(forces, axis=0).tolist(),
        "outcar_sha256": sha256(outcar),
    }


def audit(work_root: Path, output_dir: Path) -> dict:
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite Gamma audit: {output_dir}")
    manifest_path = work_root / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if (manifest["purpose"] != "GaN_1x1x1_endpoint_Gamma_force_constants_two_step_sizes"
            or manifest["supercell_matrix"] != np.eye(3, dtype=int).tolist()
            or len(manifest["cases"]) != 32):
        raise ValueError("wrong or incomplete GaN 1x1x1 manifest")
    evaluated = {}
    family_data = {}
    for phase in ("B4", "B1"):
        source = manifest["source"][phase]
        source_dir = Path(source["directory"])
        if any(sha256(source_dir / name) != expected
               for name, expected in source["sha256"].items()):
            raise ValueError(f"archived endpoint source changed: {phase}")
        reference = read(source_dir / "CONTCAR", format="vasp")
        baseline = read(source_dir / "OUTCAR")
        for distance in (0.01, 0.02):
            records = [case for case in manifest["cases"]
                       if case["phase"] == phase and case["distance_A"] == distance]
            records.sort(key=lambda item: item["local_index"])
            phonon = phonopy_for_reference(reference, distance)
            if len(records) != len(phonon.supercells_with_displacements) or len(records) != 8:
                raise ValueError(f"Phonopy displacement count changed: {phase} {distance}")
            forces = []
            details = []
            for index, case in enumerate(records):
                if case["local_index"] != index:
                    raise ValueError("Phonopy displacement order changed")
                expected_disp = phonon.dataset["first_atoms"][index]
                if (case["atom_index"] != int(expected_disp["number"])
                        or not np.allclose(case["displacement_A"], expected_disp["displacement"], atol=1e-10)):
                    raise ValueError("manifest displacement differs from regenerated Phonopy set")
                item = completed_static(work_root / "cases" / case["name"], case["input_sha256"])
                forces.append(item["forces"])
                details.append({key: val for key, val in item.items() if key != "forces"}
                               | {"case": case["name"], "input_sha256": case["input_sha256"]})
            phonon.produce_force_constants(forces=np.asarray(forces), show_drift=False)
            fc = np.asarray(phonon.force_constants, dtype=float)
            if fc.shape != (4, 4, 3, 3) or not np.isfinite(fc).all():
                raise ValueError("Phonopy did not produce a finite four-atom force-constant tensor")
            phonon.run_qpoints([[0.0, 0.0, 0.0]], with_eigenvectors=True)
            qpoint = phonon.get_qpoints_dict()
            frequency = np.asarray(qpoint["frequencies"][0], dtype=float)
            eigenvectors = real_gamma_eigenvectors(qpoint["eigenvectors"][0])
            if frequency.shape != (12,) or not np.isfinite(frequency).all():
                raise ValueError("Phonopy Gamma spectrum is invalid")
            label = f"{phase}_d{distance:.2f}"
            family_data[label] = {
                "force_constants": fc.transpose(0, 2, 1, 3).copy(),
                "phonopy_force_constants": fc.copy(),
                "masses_amu": np.asarray(phonon.primitive.masses, dtype=float),
                "frequencies_thz": frequency,
                "eigenvectors_mass_weighted": eigenvectors,
            }
            first_forces = np.asarray(forces)
            evaluated[label] = {
                "n_displacements": len(records),
                "max_reference_force_eV_per_A": float(np.max(np.linalg.norm(baseline.get_forces(), axis=1))),
                "max_displaced_net_force_eV_per_A": float(np.max(np.linalg.norm(np.sum(first_forces, axis=1), axis=1))),
                "max_force_constant_asr_drift_eV_per_A2": float(np.max(np.abs(np.sum(fc, axis=1)))),
                "max_abs_acoustic_frequency_thz": float(np.max(np.abs(frequency[:3]))),
                "max_force_constant_permutation_asymmetry_eV_per_A2": float(
                    np.max(np.abs(fc.transpose(0, 2, 1, 3).reshape(12, 12)
                                  - fc.transpose(0, 2, 1, 3).reshape(12, 12).T))
                ),
                "frequencies_thz": frequency.tolist(),
                "cases": details,
            }
    comparisons = {}
    for phase in ("B4", "B1"):
        small = family_data[f"{phase}_d0.01"]
        large = family_data[f"{phase}_d0.02"]
        difference = small["phonopy_force_constants"] - large["phonopy_force_constants"]
        comparisons[phase] = {
            "max_fc_step_size_difference_eV_per_A2": float(np.max(np.abs(difference))),
            "relative_fc_frobenius_difference": float(
                np.linalg.norm(difference) / max(np.linalg.norm(small["phonopy_force_constants"]), 1e-15)
            ),
            "max_sorted_frequency_difference_thz": float(
                np.max(np.abs(small["frequencies_thz"] - large["frequencies_thz"]))
            ),
            "max_optical_frequency_difference_thz": float(
                np.max(np.abs(small["frequencies_thz"][3:] - large["frequencies_thz"][3:]))
            ),
        }
    output_dir.mkdir(parents=True)
    for label, arrays in family_data.items():
        np.savez_compressed(output_dir / f"{label}.npz", **{
            key: value for key, value in arrays.items() if key != "phonopy_force_constants"
        })
    report = {
        "status": "computed_GaN_endpoint_atomic_Gamma_1x1x1_requires_numerical_interpretation",
        "method": "Phonopy 1x1x1 finite displacements, VASP independent statics, two step sizes",
        "phonopy_version": __import__("phonopy").__version__,
        "source_manifest_sha256": sha256(manifest_path),
        "source_code_sha256": {
            "preparer": sha256(work_root.parent / "prepare_gan_gamma_phonopy_1x1x1.py"),
            "launcher": sha256(work_root.parent / "hf_gan_gamma_phonopy_1x1x1.slurm"),
            "auditor": sha256(Path(__file__)),
        },
        "families": evaluated,
        "step_size_comparison": comparisons,
        "limitations": manifest["limitations"],
    }
    (output_dir / "audit.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    report = audit(args.work_root, args.output_dir)
    print(json.dumps({"status": report["status"],
                      "step_size_comparison": report["step_size_comparison"]}))


if __name__ == "__main__":
    main()

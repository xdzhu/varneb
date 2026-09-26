"""Stage one conservative GaN TS trial inside an observed VASP lattice-failure bracket.

This is a case-local empirical guard, not a general replacement for VASP's
Bravais classifier. It keeps the physical VASP input contract fixed and never
modifies the failed full-step POSCAR.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil

import numpy as np
from ase.io import read, write
from ase.units import GPa

from scripts.audit_gan_joint_curvature import completed_case
from scripts.audit_gan_ts_newton_probe import failed_bravais_case
from scripts.prepare_gan_ts_newton_probe import same_geometry, sha256
from vcneb.joint_curvature import JointCurvatureCoordinates


def cross_axis_skew_angstrom(cell: np.ndarray) -> float:
    """Lateral projection of the third vector onto either basal vector."""

    vectors = np.asarray(cell, dtype=float)
    if vectors.shape != (3, 3) or not np.isfinite(vectors).all():
        raise ValueError("invalid cell")
    lengths = np.linalg.norm(vectors[:2], axis=1)
    if np.any(lengths <= 0) or np.linalg.det(vectors) <= 0:
        raise ValueError("non-positive or degenerate cell")
    return float(np.max(np.abs(vectors[:2] @ vectors[2]) / lengths))


def choose_guarded_trial(
    coordinates: JointCurvatureCoordinates, correction: np.ndarray,
    skews: dict[str, float],
) -> tuple[float, object, float, float, list[dict]]:
    """Backtrack within an observed complete/failure bracket with a margin."""

    if not skews["center"] < skews["completed_half"] < skews["failed_full"]:
        raise ValueError("the observed metric risk bracket is not ordered")
    # One-case empirical backtracking guard: require a 10% margin under the
    # mid-skew of the previous complete/failed cells. It cannot guarantee
    # VASP will initialize, but rejects the already dangerous 0.75 step.
    limit = 0.9 * (skews["completed_half"] + skews["failed_full"]) / 2.0
    options = []
    chosen = None
    for fraction in (0.75, 0.625, 0.5625):
        atoms = coordinates.displaced(fraction * correction)
        skew = cross_axis_skew_angstrom(atoms.cell.array)
        options.append({"fraction": fraction, "cross_axis_skew_A": skew,
                        "passes_empirical_margin": bool(skew < limit)})
        if chosen is None and skew < limit:
            chosen = (fraction, atoms, skew)
    if chosen is None:
        raise ValueError("no trial lies inside the empirical Bravais-risk margin")
    return *chosen, limit, options


def prepare(trajectory: Path, previous_work: Path, previous_manifest: Path,
            previous_audit: Path, center_1000: Path, work_root: Path) -> dict:
    if work_root.exists():
        raise FileExistsError(work_root)
    manifest = json.loads(previous_manifest.read_text(encoding="utf-8"))
    audit = json.loads(previous_audit.read_text(encoding="utf-8"))
    raw_manifest = previous_work / "manifest.json"
    if (sha256(raw_manifest) != sha256(previous_manifest)
            or audit["source_sha256"]["manifest"] != sha256(previous_manifest)
            or audit["source_sha256"]["trajectory"] != sha256(trajectory)
            or manifest["source_sha256"]["center_1000_OUTCAR"] != sha256(center_1000 / "OUTCAR")
            or manifest["evaluation_encut_eV"] != 1000
            or manifest["pressure_GPa"] != 45.7
            or manifest["status"] != "inputs_finalized_no_DFT"
            or audit["best_trial"] != "newton_half"
            or [case["name"] for case in manifest["cases"]] != ["newton_half", "newton_full"]):
        raise ValueError("source Newton probe or raw DFT provenance failed")
    reference = read(trajectory, index=15)
    coordinates = JointCurvatureCoordinates(reference, float(manifest["cell_scale_A"]))
    correction = np.asarray(manifest["correction_A"], dtype=float)
    if correction.shape != (coordinates.size,) or not np.isfinite(correction).all():
        raise ValueError("invalid stored Newton correction")
    source = read(center_1000 / "OUTCAR")
    if not same_geometry(reference, source, tolerance=2e-5):
        raise ValueError("1000 eV source center is not trajectory image 15")
    prior = {case["name"]: case for case in manifest["cases"]}
    successful, _ = completed_case(previous_work / "cases/newton_half", prior["newton_half"])
    failed = failed_bravais_case(previous_work / "cases/newton_full", prior["newton_full"])
    if (failed["outcar_sha256"] != audit["cases"][1]["outcar_sha256"]
            or sha256(previous_work / "cases/newton_half/OUTCAR") != audit["cases"][0]["outcar_sha256"]):
        raise ValueError("previous success/failure evidence changed")
    pressure = 45.7 * GPa
    gradient = coordinates.enthalpy_gradient(
        successful, successful.get_forces(), successful.get_stress(voigt=False), pressure,
    )
    norm = float(np.linalg.norm(coordinates.translation_free_basis().T @ gradient))
    if not np.isclose(norm, audit["cases"][0]["translation_free_gradient_eV_per_A"],
                      rtol=0, atol=1e-7):
        raise ValueError("successful half-step gradient differs from audit")

    cells = {"center": reference.cell.array,
             "completed_half": read(previous_work / "cases/newton_half/POSCAR",
                                    format="vasp").cell.array,
             "failed_full": read(previous_work / "cases/newton_full/POSCAR",
                                 format="vasp").cell.array}
    skews = {name: cross_axis_skew_angstrom(cell) for name, cell in cells.items()}
    fraction, atoms, skew, limit, options = choose_guarded_trial(
        coordinates, correction, skews,
    )
    distances = atoms.get_all_distances(mic=True)
    np.fill_diagonal(distances, np.inf)
    if len(atoms) != 4 or float(np.min(distances)) < 1.4:
        raise ValueError("unsafe candidate atomic geometry")
    incar = (center_1000 / "INCAR").read_text(encoding="utf-8")
    kpoints = (center_1000 / "KPOINTS").read_text(encoding="utf-8")
    if ("ENCUT = 1000.000000" not in incar or "ISYM = -1" not in incar
            or "SYMPREC = 1.00e-04" not in incar or "8 8 6" not in kpoints):
        raise ValueError("1000 eV static input contract changed")
    for filename in ("INCAR", "KPOINTS", "POTCAR"):
        if sha256(center_1000 / filename) != manifest["source_static_sha256"][filename]:
            raise ValueError(f"1000 eV source input changed: {filename}")

    name = f"newton_{str(fraction).replace('.', 'p')}"
    directory = work_root / "cases" / name
    directory.mkdir(parents=True)
    write(directory / "POSCAR", atoms, format="vasp", direct=True, sort=False)
    if not same_geometry(atoms, read(directory / "POSCAR", format="vasp")):
        raise ValueError("POSCAR roundtrip changed trial geometry")
    for filename in ("INCAR", "KPOINTS", "POTCAR"):
        shutil.copy2(center_1000 / filename, directory / filename)
    hashes = {filename: sha256(directory / filename)
              for filename in ("POSCAR", "INCAR", "KPOINTS", "POTCAR")}
    (directory / "sha256.inputs.json").write_text(
        json.dumps(hashes, indent=2) + "\n", encoding="utf-8",
    )
    result = {
        "purpose": "GaN_image15_bracketed_local_Newton_trial_not_TS_certificate",
        "status": "inputs_finalized_no_DFT",
        "n_atoms": 4,
        "pressure_GPa": 45.7,
        "evaluation_encut_eV": 1000,
        "source_hessian_encut_eV": 600,
        "cell_scale_A": manifest["cell_scale_A"],
        "empirical_bravais_risk": {
            "observable": "max_abs_basal_third_vector_projection_A",
            "observed_skew_A": skews,
            "allowed_skew_A": limit,
            "candidate_options": options,
            "warning": "Case-local heuristic; VASP itself may still reject the trial.",
        },
        "cases": [{"name": name, "fraction_from_original_center": fraction,
                   "cross_axis_skew_A": skew,
                   "minimum_distance_A": float(np.min(distances)),
                   "volume_A3": float(atoms.get_volume()),
                   "POSCAR_sha256": hashes["POSCAR"],
                   "input_sha256": hashes}],
        "source_sha256": {
            "trajectory": sha256(trajectory),
            "previous_manifest": sha256(previous_manifest),
            "previous_audit": sha256(previous_audit),
            "half_OUTCAR": sha256(previous_work / "cases/newton_half/OUTCAR"),
            "full_failed_OUTCAR": sha256(previous_work / "cases/newton_full/OUTCAR"),
            "center_1000_OUTCAR": sha256(center_1000 / "OUTCAR"),
            "preparer": sha256(Path(__file__)),
        },
        "limitations": "One guarded local trial, no automatic physical-input changes or TS claim.",
    }
    (work_root / "manifest.json").write_text(json.dumps(result, indent=2) + "\n",
                                              encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trajectory", type=Path, required=True)
    parser.add_argument("--previous-work", type=Path, required=True)
    parser.add_argument("--previous-manifest", type=Path, required=True)
    parser.add_argument("--previous-audit", type=Path, required=True)
    parser.add_argument("--center-1000", type=Path, required=True)
    parser.add_argument("--work-root", type=Path, required=True)
    args = parser.parse_args()
    result = prepare(args.trajectory, args.previous_work, args.previous_manifest,
                     args.previous_audit, args.center_1000, args.work_root)
    print(json.dumps({"status": result["status"], "case": result["cases"][0]["name"],
                      "risk": result["empirical_bravais_risk"]}))


if __name__ == "__main__":
    main()

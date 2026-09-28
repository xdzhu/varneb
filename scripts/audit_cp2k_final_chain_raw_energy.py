"""Read-only CP2K final-chain energy/stress reconciliation against raw text.

The final 29 frames are matched to exact-geometry worker cache entries.
The last FORCE_EVAL and stress blocks in each image's cumulative output are
then compared with those cache entries under one empirically fitted unit
factor. This does not certify unprinted raw atomic forces or missing run-end
markers. No calculator is launched and no file is written.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import numpy as np
from ase.io import read


ENERGY_RE = re.compile(
    r"^\s*ENERGY\| Total FORCE_EVAL \( QS \) energy \[a\.u\.\]:\s*"
    r"([-+]?\d+(?:\.\d*)?(?:[Ee][-+]?\d+)?)\s*$",
    re.MULTILINE,
)
STRESS_ROW_RE = re.compile(
    r"^\s*STRESS\|\s*([xyz])\s+"
    r"([-+]?\d+(?:\.\d*)?[Ee][-+]?\d+)\s+"
    r"([-+]?\d+(?:\.\d*)?[Ee][-+]?\d+)\s+"
    r"([-+]?\d+(?:\.\d*)?[Ee][-+]?\d+)\s*$",
    re.MULTILINE,
)


def cache_digest(index: int, image) -> str:
    """Reproduce the exact key used by ThreadedCalculatorExecutor."""
    digest = hashlib.sha256()
    digest.update(str(index).encode("ascii"))
    digest.update(b"\0")
    digest.update("\0".join(image.get_chemical_symbols()).encode("utf-8"))
    for value in (
        np.asarray(image.get_positions(), dtype="<f8", order="C"),
        np.asarray(image.cell.array, dtype="<f8", order="C"),
        np.asarray(image.pbc, dtype="<?", order="C"),
    ):
        digest.update(np.asarray(value).tobytes(order="C"))
        digest.update(b"\0")
    return digest.hexdigest()


def last_stress_gpa(text: str) -> np.ndarray:
    rows = [(axis, float(x), float(y), float(z))
            for axis, x, y, z in STRESS_ROW_RE.findall(text)]
    if len(rows) < 3 or [row[0] for row in rows[-3:]] != ["x", "y", "z"]:
        raise ValueError("last CP2K stress block is missing or incomplete")
    return np.asarray([row[1:] for row in rows[-3:]], dtype=float)


def audit(root: Path, trajectory: Path) -> dict:
    images = read(trajectory, index="-29:")
    if len(images) != 29:
        raise ValueError(f"expected 29 final images, got {len(images)}")
    records = []
    for index in range(1, 28):
        image = images[index]
        digest = cache_digest(index, image)
        cache_path = root / "image_cache" / f"image_{index:04d}_{digest}.npz"
        if not cache_path.is_file():
            raise FileNotFoundError(f"exact final-image cache missing: {cache_path}")
        with np.load(cache_path, allow_pickle=False) as data:
            energy = float(np.asarray(data["energy"]).reshape(()))
            forces = np.asarray(data["forces"])
            stress = np.asarray(data["stress"])
        chain_energy = float(image.get_potential_energy())
        chain_forces = np.asarray(image.get_forces())
        chain_stress = np.asarray(image.get_stress(voigt=False))
        cache_differences = {
            "energy_ev": abs(energy - chain_energy),
            "force_ev_per_angstrom": float(np.max(np.abs(forces - chain_forces))),
            "stress_ev_per_angstrom3": float(np.max(np.abs(stress - chain_stress))),
        }
        if max(cache_differences.values()) > 1e-8:
            raise ValueError(f"image {index} cache differs from final chain")

        output_path = root / f"image_{index:04d}" / "cp2k.out"
        raw = output_path.read_bytes()
        text = raw.decode("utf-8", errors="replace")
        energies = ENERGY_RE.findall(text)
        if not energies:
            raise ValueError(f"image {index} lacks FORCE_EVAL energy")
        last_energy_offset = text.rfind("ENERGY| Total FORCE_EVAL")
        last_start_offset = text.rfind("PROGRAM STARTED AT", 0, last_energy_offset)
        records.append({
            "image": index,
            "raw_output_sha256": hashlib.sha256(raw).hexdigest(),
            "exact_cache_name": cache_path.name,
            "cache_sha256": hashlib.sha256(cache_path.read_bytes()).hexdigest(),
            "cache_chain_max_differences": cache_differences,
            "n_force_eval_records": len(energies),
            "last_raw_energy_hartree": float(energies[-1]),
            "cache_energy_ev": energy,
            "last_raw_stress_gpa": last_stress_gpa(text).tolist(),
            "cache_stress_ev_per_angstrom3": stress.tolist(),
            "last_run_scf_converged": (
                last_start_offset < text.rfind("SCF run converged", 0, last_energy_offset)
            ),
            "last_run_has_program_end": (
                text.rfind("PROGRAM ENDED AT") > last_energy_offset
            ),
        })

    raw_energy = np.asarray([row["last_raw_energy_hartree"] for row in records])
    cache_energy = np.asarray([row["cache_energy_ev"] for row in records])
    # Lock the conversion on image 15, then test the other 26 images without
    # refitting. This prevents a whole-chain fit from concealing outliers.
    reference = 14
    energy_factor = float(cache_energy[reference] / raw_energy[reference])
    energy_residuals = cache_energy - energy_factor * raw_energy
    raw_stress = np.asarray([row["last_raw_stress_gpa"] for row in records])
    cache_stress = np.asarray([row["cache_stress_ev_per_angstrom3"]
                               for row in records])
    stress_factor = float(
        np.sum(raw_stress[reference] * cache_stress[reference])
        / np.sum(raw_stress[reference] * raw_stress[reference])
    )
    stress_residuals = cache_stress - stress_factor * raw_stress
    max_energy_residual = float(np.max(np.abs(energy_residuals)))
    max_stress_residual = float(np.max(np.abs(stress_residuals)))
    n_converged = sum(row["last_run_scf_converged"] for row in records)
    if max_energy_residual > 1e-8 or max_stress_residual > 1e-9 or n_converged != 27:
        raise ValueError("CP2K raw energy/stress or SCF gate failed")
    return {
        "status": "27_interior_raw_energy_stress_reconciled_no_raw_forces_or_end_markers",
        "trajectory_sha256": hashlib.sha256(trajectory.read_bytes()).hexdigest(),
        "n_total_images": len(images),
        "n_interior_images": len(records),
        "n_raw_last_scf_converged": n_converged,
        "n_raw_last_program_end": sum(row["last_run_has_program_end"] for row in records),
        "conversion_reference_image": 15,
        "energy_hartree_to_ev_factor": energy_factor,
        "max_absolute_energy_residual_after_factor_ev": max_energy_residual,
        "stress_gpa_to_ase_factor": stress_factor,
        "max_absolute_stress_residual_after_factor_ev_per_angstrom3": max_stress_residual,
        "max_absolute_cache_chain_difference": {
            key: max(row["cache_chain_max_differences"][key] for row in records)
            for key in records[0]["cache_chain_max_differences"]
        },
        "records": records,
        "limitations": (
            "Raw atomic forces are not printed in the cumulative CP2K text. "
            "A converged SCF marker and final stress/energy block do not replace "
            "a missing PROGRAM ENDED marker. Image-15 conversion factors are "
            "empirical reconciliations, not independently sourced physical constants."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("trajectory", type=Path)
    parser.add_argument("--summary", action="store_true")
    args = parser.parse_args()
    report = audit(args.root, args.trajectory)
    if args.summary:
        report.pop("records")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()

"""Audit isolated B4 endpoint cell relaxations without replacing VCNEB data."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from ase.io import read
from ase.units import GPa

from scripts.audit_gan_45p7_vasp_endpoints import (
    DEFAULT_CHAIN, EXPECTED_CHAIN_SHA256, PRESSURE_GPA,
)
from scripts.audit_gan_ts_basin_followups import residual_stress_kbar
from scripts.compare_gan_basin_endpoint import relative_displacement_metrics


ALLOWED_PURPOSES = {
    "GaN_45p7_600eV_original_B4_endpoint_cell_relax_not_path_replacement",
    "GaN_45p7_600eV_original_B4_endpoint_strict_ionic_continuation",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit(source_dir: Path, chain_path: Path = DEFAULT_CHAIN,
          slurm_job: int | None = None) -> dict:
    source_dir = Path(source_dir)
    chain_path = Path(chain_path)
    if sha256(chain_path) != EXPECTED_CHAIN_SHA256:
        raise ValueError("reference VASP final chain changed")
    manifest = json.loads((source_dir / "manifest.json").read_text(encoding="utf-8"))
    if (manifest.get("status") != "inputs_finalized_no_DFT"
            or manifest.get("purpose") not in ALLOWED_PURPOSES
            or manifest.get("pressure_GPa") != PRESSURE_GPA
            or manifest.get("electronic_contract", {}).get("ENCUT_eV") != 600):
        raise ValueError("not a pinned 600-eV original-B4 endpoint continuation")
    hashes = json.loads((source_dir / "sha256.inputs.json").read_text(encoding="utf-8"))
    for filename in ("INCAR", "KPOINTS", "POSCAR"):
        if (hashes.get(filename) != manifest["input_sha256"].get(filename)
                or sha256(source_dir / filename) != hashes[filename]):
            raise ValueError(f"audited continuation input changed: {filename}")
    if hashes.get("POTCAR") != manifest["input_sha256"].get("POTCAR"):
        raise ValueError("POTCAR identity differs between manifest and input hashes")
    raw_path = source_dir / "OUTCAR"
    raw_text = raw_path.read_text(encoding="utf-8", errors="replace")
    if ("aborting loop because EDIFF is reached" not in raw_text
            or "General timing and accounting informations for this job"
               not in raw_text):
        raise ValueError("continuation raw output lacks SCF or completed-run marker")
    frames = read(raw_path, index=":")
    if not frames:
        raise ValueError("no evaluated ionic geometry in continuation")
    final = frames[-1]
    contcar = read(source_dir / "CONTCAR", format="vasp")
    start = read(source_dir / "POSCAR", format="vasp")
    original = read(chain_path, index=0)
    if (final.get_chemical_symbols() != ["Ga", "Ga", "N", "N"]
            or contcar.get_chemical_symbols() != final.get_chemical_symbols()
            or start.get_chemical_symbols() != final.get_chemical_symbols()):
        raise ValueError("GaN endpoint atom ordering changed")
    contcar_cell_diff = float(np.max(np.abs(contcar.cell.array - final.cell.array)))
    contcar_position_diff = float(np.max(np.abs(contcar.positions - final.positions)))
    if contcar_cell_diff > 1e-8 or contcar_position_diff > 1e-5:
        raise ValueError("CONTCAR differs from the final evaluated geometry")
    original_h = original.get_potential_energy() + PRESSURE_GPA * GPa * original.get_volume()
    final_h = final.get_potential_energy() + PRESSURE_GPA * GPa * final.get_volume()
    fmax = float(np.linalg.norm(final.get_forces(), axis=1).max())
    stress = residual_stress_kbar(final.get_stress(voigt=False), PRESSURE_GPA)
    return {
        "status": "GaN_45p7_VASP_original_B4_endpoint_relax_raw_audited_not_path_replacement",
        "purpose": manifest["purpose"],
        "slurm_job": slurm_job,
        "n_evaluated_ionic_frames": len(frames),
        "vasp_reports_ionic_convergence":
            "reached required accuracy - stopping structural energy minimisation" in raw_text,
        "pressure_GPa": PRESSURE_GPA,
        "ENCUT_eV": 600,
        "final_energy_eV_per_cell": float(final.get_potential_energy()),
        "final_volume_A3": float(final.get_volume()),
        "final_enthalpy_eV_per_cell": float(final_h),
        "final_minus_original_H_meV_per_GaN": float((final_h - original_h) * 500),
        "final_max_atomic_force_eV_per_A": fmax,
        "final_max_raw_stress_residual_kbar": stress,
        "force_pass_0p02eV_per_A": fmax <= 0.02,
        "stress_pass_2kbar": stress <= 2.0,
        "start_to_final_cell_max_abs_A": float(np.max(np.abs(
            start.cell.array - final.cell.array))),
        "start_to_final_position_max_abs_A": float(np.max(np.abs(
            start.positions - final.positions))),
        "final_vs_original_B4": relative_displacement_metrics(final, original),
        "contcar_vs_final_cell_max_abs_A": contcar_cell_diff,
        "contcar_vs_final_position_max_abs_A": contcar_position_diff,
        "source_sha256": {
            "manifest": sha256(source_dir / "manifest.json"),
            "input_hashes": sha256(source_dir / "sha256.inputs.json"),
            "OUTCAR": sha256(raw_path), "CONTCAR": sha256(source_dir / "CONTCAR"),
            "reference_chain": sha256(chain_path),
            "auditor": sha256(Path(__file__)),
        },
        "limitations": [
            "This is not a recalculated VCNEB barrier or a replacement production endpoint.",
            "POTCAR is not redistributed; its hash is checked in the remote input manifest, not locally against POTCAR bytes.",
            "A native stop and a 2-kbar pressure gate are distinct conditions.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--chain", type=Path, default=DEFAULT_CHAIN)
    parser.add_argument("--slurm-job", type=int)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = audit(args.source_dir, args.chain, args.slurm_job)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in (
        "status", "slurm_job", "n_evaluated_ionic_frames",
        "final_max_atomic_force_eV_per_A", "final_max_raw_stress_residual_kbar",
        "final_minus_original_H_meV_per_GaN",
    )}))


if __name__ == "__main__":
    main()

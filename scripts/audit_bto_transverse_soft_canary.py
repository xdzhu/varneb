"""Independently audit the three BTO soft/soft static branch seeds.

This reads raw ABACUS files but never runs a calculator. A completed canary
only establishes static branch energies; it is not a conditional minimum, PES,
or T-to-C barrier.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.audit_bto_q1q2_conditional_pilot import _raw_force_stress


FINAL_ENERGY = re.compile(r"!FINAL_ETOT_IS\s*=?\s*([-+\d.eE]+)\s*eV")
DSIZE = re.compile(r"\bDSIZE\s*=\s*(\d+)")
SOURCE_NAMES = {
    "report": "bto_t_to_c_q1q2_strain_gradient_preflight_2026-09-25.json",
    "reference": "cubic_CONTCAR",
    "force_constants": "bto_cubic_gamma_force_constants.npz",
    "phonopy_eigenpairs": "bto_cubic_phonopy_gamma_eigenpairs.npz",
    "gamma_provenance": "bto_cubic_gamma_phonon_provenance.json",
    "force_sets": "bto_cubic_gamma_FORCE_SETS",
    "eigenpairs_provenance": "bto_cubic_phonopy_gamma_eigenpairs_provenance.json",
    "grid_result_manifest": "result_manifest_hf_27777454.json",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def branch_energy_diagnostics(rows: list[dict]) -> dict:
    """Separate raw-point validity from the sign of branch energy changes."""
    if [row.get("label") for row in rows] != ["frozen", "+Q_y", "-Q_y"]:
        raise ValueError("expected exactly the three reviewed branch labels")
    values = [float(row["relative_energy_eV_per_BTO"]) for row in rows]
    if not np.all(np.isfinite(values)):
        raise ValueError("branch energies must be finite")
    lowered = values[0] - values[1]
    return {
        "frozen_minus_plus_branch_meV_per_BTO": 1000.0 * lowered,
        "plus_branch_lowers_frozen": bool(lowered > 0),
        "signed_branch_energy_difference_meV_per_BTO": 1000.0 * (values[1] - values[2]),
    }


def audit(inputs: Path, workdir: Path, *, expected_job_id: str,
          preflight_path: Path | None = None) -> dict:
    if preflight_path is None:
        preflight_path = inputs / "bto_transverse_soft_conditional_q060_preflight_2026-09-26.json"
    if preflight_path.resolve().parent != inputs.resolve():
        raise ValueError("preflight must be a direct child of the audited inputs directory")
    canary_path = workdir / "canary_starts_result.json"
    preflight = json.loads(preflight_path.read_text(encoding="utf-8"))
    canary = json.loads(canary_path.read_text(encoding="utf-8"))
    grid = json.loads((inputs / SOURCE_NAMES["grid_result_manifest"]).read_text(encoding="utf-8"))
    if (preflight.get("kind") != "bto_transverse_soft_conditional_preflight_no_dft"
            or preflight.get("phonon_supercell") != [1, 1, 1]
            or preflight.get("electronic_kpoints") != [4, 4, 4]
            or preflight.get("n_relaxed_orthogonal_coordinates") != 16
            or canary.get("kind") != "bto_transverse_soft_three_start_canary_not_PES_or_barrier"
            or canary.get("status") != "three_start_statics_complete_optimizer_not_run"
            or canary.get("preflight_sha256") != sha256(preflight_path)
            or str(canary.get("slurm_job_id")) != expected_job_id
            or canary.get("n_new_evaluations") != 3):
        raise ValueError("canary, preflight, or Slurm identity is inconsistent")
    source_hashes = {}
    for key, name in SOURCE_NAMES.items():
        path = inputs / name
        digest = sha256(path)
        if preflight["source_sha256"].get(key) != digest:
            raise ValueError(f"preflight source hash changed: {key}")
        source_hashes[key] = digest
    gamma = json.loads((inputs / SOURCE_NAMES["gamma_provenance"]).read_text(encoding="utf-8"))
    eigen = json.loads((inputs / SOURCE_NAMES["eigenpairs_provenance"]).read_text(encoding="utf-8"))
    if (gamma.get("supercell") != [1, 1, 1]
            or gamma.get("calculator", {}).get("ecutwfc_Ry") != 100
            or gamma.get("calculator", {}).get("kpoints") != [4, 4, 4]
            or eigen.get("force_sets_sha256") != source_hashes["force_sets"]):
        raise ValueError("Gamma source does not independently establish the BTO contract")
    parameters = grid.get("calculator_parameters", {})
    if (grid.get("status") != "converged" or grid.get("mpi_ranks") != 32
            or parameters.get("ecutwfc") != 100.0
            or parameters.get("kpts") != [4, 4, 4]
            or "Orb-DZP-10au" not in parameters.get("basis_dir", "")):
        raise ValueError("archived frozen-grid calculator contract changed")
    cubic = [point for point in grid["points"]
             if abs(point["q1"]) < 1e-12 and abs(point["q2"]) < 1e-12]
    if len(cubic) != 1:
        raise ValueError("frozen grid must have a unique cubic reference")
    cubic_energy = float(cubic[0]["energy_eV"])
    expected_input = grid["points"][0]["input_sha256"]
    labels = ["frozen", "+Q_y", "-Q_y"]
    if ([row.get("label") for row in preflight["branch_starts"]] != labels
            or [row.get("label") for row in canary["points"]] != labels):
        raise ValueError("canary branch labels or order changed")
    rows = []
    for seed, item in zip(preflight["branch_starts"], canary["points"]):
        directory = workdir / item["evaluation_directory"]
        cached_path = directory / "result.json"
        cached = json.loads(cached_path.read_text(encoding="utf-8"))
        coordinates = np.asarray(seed["coordinates_u_A_eta_voigt"], dtype=float)
        if (not np.array_equal(np.asarray(item["coordinates_u_A_eta_voigt"], dtype=float), coordinates)
                or not np.array_equal(np.asarray(cached["coordinates"], dtype=float), coordinates)
                or item["cached_result_sha256"] != sha256(cached_path)
                or cached.get("status") != "complete"
                or cached.get("contract_sha256") != canary["evaluator_contract_sha256"]
                or cached.get("validation", {}).get("mpi_dsize") != 32
                or str(cached.get("validation", {}).get("slurm_job_id")) != expected_job_id):
            raise ValueError(f"cached DFT point does not match signed seed: {item['label']}")
        raw_hashes = {}
        for name in ("INPUT", "KPT", "STRU"):
            digest = sha256(directory / name)
            if cached["validation"]["input_sha256"][name] != digest:
                raise ValueError(f"raw {name} hash changed: {item['label']}")
            raw_hashes[name] = digest
        if raw_hashes["INPUT"] != expected_input["INPUT"] or raw_hashes["KPT"] != expected_input["KPT"]:
            raise ValueError(f"ABACUS INPUT/KPT differs from the 100Ry frozen grid: {item['label']}")
        log_path = directory / "OUT.ABACUS" / "running_scf.log"
        log = log_path.read_text(encoding="utf-8", errors="replace")
        energies = FINAL_ENERGY.findall(log)
        sizes = DSIZE.findall(log)
        if (cached["validation"]["log_sha256"] != sha256(log_path)
                or "charge density convergence is achieved" not in log
                or len(energies) != 1 or sizes != ["32"]):
            raise ValueError(f"SCF, final energy, or MPI markers missing: {item['label']}")
        energy = float(cached["enthalpy_eV"])
        gradient = np.asarray(cached["gradient"], dtype=float)
        stress = np.asarray(cached["stress_eV_per_A3"], dtype=float)
        raw_forces, raw_stress_kbar = _raw_force_stress(
            log, symbols=["Ba", "Ti", "O", "O", "O"],
        )
        raw_force_max = float(np.max(np.linalg.norm(raw_forces, axis=1)))
        raw_stress_eV_per_A3 = -raw_stress_kbar / 1602.176634
        force_gradient_disagreement = float(np.max(np.abs(
            -raw_forces.ravel() - gradient[:15],
        )))
        stress_disagreement_kbar = float(np.max(np.abs(
            raw_stress_eV_per_A3 - stress,
        )) * 1602.176634)
        if (not np.isfinite(energy) or gradient.shape != (21,) or not np.all(np.isfinite(gradient))
                or stress.shape != (3, 3) or not np.all(np.isfinite(stress))
                or abs(float(energies[0]) - energy) > 1e-6
                or abs(energy - cubic_energy - item["energy_minus_c_eV_per_BTO"]) > 1e-8
                or abs(raw_force_max - cached["maximum_atomic_force_eV_per_A"]) > 1e-5
                or force_gradient_disagreement > 1e-5
                or stress_disagreement_kbar > 1e-4
                or cached["minimum_distance_A"] < preflight["minimum_allowed_atomic_distance_A"]
                or cached["volume_A3"] <= 0.0):
            raise ValueError(f"nonphysical or inconsistent raw energy/gradient/stress: {item['label']}")
        rows.append({
            "label": item["label"],
            "evaluation_directory": directory.name,
            "cached_result_sha256": sha256(cached_path),
            "raw_log_sha256": sha256(log_path),
            "raw_input_sha256": raw_hashes,
            "relative_energy_eV_per_BTO": energy - cubic_energy,
            "max_atomic_force_eV_per_A": raw_force_max,
            "max_absolute_stress_kbar": float(np.max(np.abs(raw_stress_kbar))),
            "max_force_gradient_disagreement_eV_per_A": force_gradient_disagreement,
            "max_stress_disagreement_kbar": stress_disagreement_kbar,
            "minimum_distance_A": float(cached["minimum_distance_A"]),
        })
    branch_diagnostics = branch_energy_diagnostics(rows)
    return {
        "kind": "independently_audited_BTO_transverse_soft_three_start_statics_not_conditional_PES",
        "status": "all_three_raw_DFT_points_validated",
        "slurm_job_id": expected_job_id,
        "source_sha256": {**source_hashes, "preflight": sha256(preflight_path),
                          "canary_result": sha256(canary_path)},
        "q_parallel_q_transverse_sqrt_amu_A": preflight["q_parallel_q_transverse_sqrt_amu_A"],
        "phonon_supercell": [1, 1, 1],
        "electronic_kpoints": [4, 4, 4],
        **branch_diagnostics,
        "points": rows,
        "limitations": "three unrelaxed statics only; no orthogonal stationarity, Hessian, PES or barrier claim",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--preflight", type=Path,
                        help="reviewed preflight JSON directly inside --inputs")
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--job-id", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(f"refusing to overwrite an existing audit: {args.output}")
    report = audit(args.inputs, args.workdir, expected_job_id=args.job_id,
                   preflight_path=args.preflight)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({key: report[key] for key in (
        "status", "frozen_minus_plus_branch_meV_per_BTO",
        "signed_branch_energy_difference_meV_per_BTO",
    )}, indent=2))


if __name__ == "__main__":
    main()

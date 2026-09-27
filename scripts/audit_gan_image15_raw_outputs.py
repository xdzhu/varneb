"""Re-read the archived GaN peak raw outputs against the evaluated chain.

This is a narrow image-15 audit, not a full-chain electronic-log audit. QE
and ABINIT use ASE's output readers; CP2K's cumulative text output is
screened conservatively because ASE's shell calculator does not expose an
equivalent offline text reader.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re

import numpy as np
from ase.io import read
from ase.units import Hartree


ROOT = Path(__file__).resolve().parents[1]
EVIDENCE = ROOT / "paper/VARNEB_CPC/evidence"
CHAIN = EVIDENCE / "gan_45p7_final_chains_20260927"
RAW = EVIDENCE / "gan_45p7_peak_raw_outputs_20260927"
EXPECTED_HASHES = {
    "qe": "8c151175c2e85e087da3181a2c43bd403eec05e6855a5a48f5bd4003a7659cb8",
    "abinit": "eb34bbdb7435490057eb51cfe67959449de96ccd06caa299bae1b4ba5be30100",
    "cp2k": "439ceac169eb5186e09afa3f51123f15a5147e29070e337c9fc2d12d4044607e",
}


def audit_peak_outputs() -> dict:
    findings = {}
    specs = {
        "qe": ("qe.pwo", "espresso-out", "JOB DONE.", "convergence has been achieved"),
        "abinit": ("abinit.abo", "abinit-out", "Calculation completed.", "etot is converged"),
    }
    for backend, (filename, fmt, done, scf) in specs.items():
        path = RAW / filename
        raw = path.read_bytes()
        if hashlib.sha256(raw).hexdigest() != EXPECTED_HASHES[backend]:
            raise ValueError(f"{backend} raw output checksum changed")
        text = raw.decode("utf-8", errors="replace")
        if done not in text or scf not in text:
            raise ValueError(f"{backend} raw output lacks SCF/completion marker")
        output = read(path, format=fmt, index=-1)
        chain = read(CHAIN / f"gan_{backend}_45p7_final_chain.traj", index=15)
        energy_difference = float(output.get_potential_energy() - chain.get_potential_energy())
        force_difference = float(np.max(np.abs(output.get_forces() - chain.get_forces())))
        stress_difference = float(np.max(np.abs(output.get_stress() - chain.get_stress())))
        if (output.get_chemical_symbols() != chain.get_chemical_symbols()
                or abs(energy_difference) > 1e-8
                or force_difference > 1e-8 or stress_difference > 1e-8):
            raise ValueError(f"{backend} parsed raw output does not reproduce final image 15")
        findings[backend] = {
            "status": "raw_output_reproduces_final_peak_energy_force_stress",
            "raw_output_sha256": EXPECTED_HASHES[backend],
            "energy_difference_eV": energy_difference,
            "max_absolute_force_difference_eV_per_A": force_difference,
            "max_absolute_stress_difference_eV_per_A3": stress_difference,
            "scf_marker_present": True,
            "completion_marker_present": True,
            "reader": fmt,
        }
    cp2k_path = RAW / "cp2k.out"
    cp2k_raw = cp2k_path.read_bytes()
    if hashlib.sha256(cp2k_raw).hexdigest() != EXPECTED_HASHES["cp2k"]:
        raise ValueError("CP2K raw output checksum changed")
    cp2k_text = cp2k_raw.decode("utf-8", errors="replace")
    energies = re.findall(
        r"ENERGY\| Total FORCE_EVAL \( QS \) energy \[a\.u\.\]:\s*([-+0-9.Ee]+)",
        cp2k_text,
    )
    if not energies:
        raise ValueError("CP2K output lacks FORCE_EVAL energies")
    energies_ev = [float(energy) * Hartree for energy in energies]
    last_energy = energies_ev[-1]
    chain = read(CHAIN / "gan_cp2k_45p7_final_chain.traj", index=15)
    chain_energy = float(chain.get_potential_energy())
    difference = float(last_energy - chain_energy)
    nearest_index = min(
        range(len(energies_ev)), key=lambda index: abs(energies_ev[index] - chain_energy)
    )
    nearest_difference = float(energies_ev[nearest_index] - chain_energy)
    last_energy_offset = cp2k_text.rfind("ENERGY| Total FORCE_EVAL")
    last_completed_offset = cp2k_text.rfind("PROGRAM ENDED AT")
    last_started_offset = cp2k_text.rfind("PROGRAM STARTED AT", 0, last_energy_offset)
    findings["cp2k"] = {
        "status": "last_visible_text_output_does_not_certify_final_peak",
        "raw_output_sha256": EXPECTED_HASHES["cp2k"],
        "n_force_eval_energy_records": len(energies),
        "last_visible_energy_eV": last_energy,
        "final_chain_peak_energy_eV": chain_energy,
        "last_visible_energy_minus_chain_eV": difference,
        "nearest_record_one_based": nearest_index + 1,
        "nearest_record_energy_eV": energies_ev[nearest_index],
        "nearest_record_minus_chain_eV": nearest_difference,
        "n_energy_records_matching_chain_within_1e-8_eV": sum(
            abs(energy - chain_energy) < 1e-8 for energy in energies_ev
        ),
        "last_program_end_precedes_last_energy": last_completed_offset < last_energy_offset,
        "last_run_scf_converged_before_last_energy": (
            last_started_offset < cp2k_text.rfind("SCF run converged", 0, last_energy_offset)
        ),
    }
    if (findings["cp2k"]["n_energy_records_matching_chain_within_1e-8_eV"] != 0
            or not findings["cp2k"]["last_program_end_precedes_last_energy"]):
        raise ValueError("CP2K visible tail changed; this audit needs review")
    cache_path = RAW / "cp2k_image15_worker_cache.npz"
    cache_digest = hashlib.sha256(cache_path.read_bytes()).hexdigest()
    if cache_digest != "f9999926bb269008fe16815e032f1d628f1c337ae3ca3b40b4606b44cb20f1fd":
        raise ValueError("CP2K image-15 worker cache checksum changed")
    with np.load(cache_path, allow_pickle=False) as cached:
        cache_energy_difference = float(cached["energy"] - chain.get_potential_energy())
        cache_force_difference = float(np.max(np.abs(cached["forces"] - chain.get_forces())))
        cache_stress_difference = float(np.max(np.abs(
            cached["stress"] - chain.get_stress(voigt=False)
        )))
    if max(abs(cache_energy_difference), cache_force_difference, cache_stress_difference) > 1e-8:
        raise ValueError("CP2K exact-geometry worker cache differs from final peak")
    findings["cp2k"].update({
        "exact_geometry_cache_key": (
            "image_0015_cc482fbb3276f4d82daccec61e40f26c4c6d87a17e4566c10f7016afadef530d.npz"
        ),
        "worker_cache_sha256": cache_digest,
        "worker_cache_energy_difference_eV": cache_energy_difference,
        "worker_cache_max_absolute_force_difference_eV_per_A": cache_force_difference,
        "worker_cache_max_absolute_stress_difference_eV_per_A3": cache_stress_difference,
    })
    return {
        "kind": "gan_45p7_image15_original_calculator_output_audit_not_full_chain",
        "status": "qe_abinit_peak_output_matched_cp2k_last_visible_text_inconclusive",
        "findings": findings,
        "limitations": "QE/ABINIT are re-read with ASE output readers, the same parser family used during calculation; this checks the frozen raw bytes against the final chain but is not parser-independent. CP2K's exact-geometry worker cache reproduces the chain energy/forces/stress, but none of the cumulative text output's FORCE_EVAL energies matches that final-chain energy within 1e-8 eV; the last visible run lacks a following PROGRAM ENDED marker. The cache establishes the serialized worker return, not independent raw-electronic-output agreement. Other 28 images and cached endpoints are not covered.",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(f"refusing to overwrite {args.output}")
    report = audit_peak_outputs()
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": report["status"], "cp2k": report["findings"]["cp2k"]}, indent=2))


if __name__ == "__main__":
    main()

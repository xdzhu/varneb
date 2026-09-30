"""Raw-audit 208 GaN VASP statics and prospectively score the 9x9 surface.

The 17x17 sample is a frozen local joint-coordinate enthalpy cut at 45.7 GPa,
not a relaxed minimum-energy surface or a certified transition state.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re

import numpy as np
from ase.units import GPa
from scipy.interpolate import CloughTocher2DInterpolator

from scripts.audit_gan_joint_curvature import completed_case
from scripts.plot_gan_600eV_local_joint_cut import source_rows
from scripts.prepare_gan_600eV_ts_2d_refinement import coordinate_key, sha256
from scripts.prepare_gan_600eV_ts_newton_canary import PRODUCTION_INPUT_SHA256


MAX_PROSPECTIVE_ERROR_MEV_PER_GAN = 0.02
RMS_PROSPECTIVE_ERROR_MEV_PER_GAN = 0.01


def audit(work: Path, pilot_path: Path, refined5_path: Path, dense9_path: Path) -> dict:
    manifest_path = work / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    pilot = json.loads(pilot_path.read_text(encoding="utf-8"))
    refined5 = json.loads(refined5_path.read_text(encoding="utf-8"))
    dense9 = json.loads(dense9_path.read_text(encoding="utf-8"))
    if (manifest.get("status") != "inputs_finalized_no_DFT"
            or manifest.get("purpose") != "GaN_45p7_600eV_local_joint_mode_17x17_frozen_refinement"
            or manifest.get("pressure_GPa") != 45.7
            or manifest.get("formula_units_per_cell") != 2
            or manifest.get("input_contract_sha256") != PRODUCTION_INPUT_SHA256
            or manifest.get("n_prior_coordinates") != 81
            or manifest.get("n_new_coordinates") != 208
            or len(manifest.get("cases", [])) != 208
            or manifest.get("source_sha256", {}).get("pilot_audit") != sha256(pilot_path)
            or manifest.get("source_sha256", {}).get("refinement5_audit") != sha256(refined5_path)
            or manifest.get("source_sha256", {}).get("dense9_audit") != sha256(dense9_path)
            or dense9.get("status") != "GaN_600eV_local_joint_9x9_refinement_raw_audited"):
        raise ValueError("17x17 manifest or preceding 81-point audit is inconsistent")
    prior = source_rows(pilot, sha256(pilot_path), refined5, sha256(refined5_path),
                        dense9, sha256(dense9_path))
    if len(prior) != 81:
        raise ValueError("incomplete predeclared 9x9 DFT prediction grid")
    xy = np.asarray([(row["q_u_A"], row["q_v_A"]) for row in prior], dtype=float)
    old_energy = np.asarray([row["measured_delta_H_meV_per_GaN"] for row in prior], dtype=float)
    predictor = CloughTocher2DInterpolator(xy, old_energy)
    h0 = float(pilot["center_enthalpy_eV_per_cell"])
    new = []
    fresh = set()
    prior_keys = {coordinate_key(*point) for point in xy}
    for item in manifest["cases"]:
        q_u, q_v = float(item["q_u_A"]), float(item["q_v_A"])
        key = coordinate_key(q_u, q_v)
        if (key in fresh or key in prior_keys
                or any(item["input_sha256"].get(name) != digest
                       for name, digest in PRODUCTION_INPUT_SHA256.items())):
            raise ValueError("duplicate 17x17 node or changed VASP electronic contract")
        fresh.add(key)
        case_dir = work / "cases" / item["name"]
        header = (case_dir / "OUTCAR").read_text(encoding="utf-8", errors="replace")[:1200]
        if "vasp.6.3.2" not in header or not re.search(r"running on\s+32 total cores", header):
            raise ValueError(f"unexpected VASP version or MPI width at {item['name']}")
        atoms, detail = completed_case(
            case_dir,
            {**item, "POSCAR_sha256": item["input_sha256"]["POSCAR"]},
        )
        h = float(atoms.get_potential_energy() + 45.7 * GPa * atoms.get_volume())
        measured = float((h - h0) * 500.0)
        predicted = float(predictor(q_u, q_v))
        if not np.isfinite([h, measured, predicted]).all():
            raise ValueError(f"nonfinite 9x9 prediction at {item['name']}")
        new.append({**detail,
                    "grid_index": item["grid_index"], "q_u_A": q_u, "q_v_A": q_v,
                    "enthalpy_eV_per_cell": h,
                    "delta_enthalpy_meV_per_GaN": measured,
                    "prior_9x9_interpolation_meV_per_GaN": predicted,
                    "prospective_error_meV_per_GaN": predicted - measured})
    full = {coordinate_key(float(u), float(v))
            for u in manifest["q_u_grid_A"] for v in manifest["q_v_grid_A"]}
    if (len(full) != 289 or len(fresh) != 208 or fresh | prior_keys != full):
        raise ValueError("17x17 grid is incomplete or overlaps existing DFT coordinates")
    errors = np.asarray([item["prospective_error_meV_per_GaN"] for item in new])
    max_error = float(np.max(np.abs(errors)))
    rms_error = float(np.sqrt(np.mean(errors**2)))
    worst = new[int(np.argmax(np.abs(errors)))]
    return {
        "status": "GaN_600eV_local_joint_17x17_refinement_raw_audited",
        "pressure_GPa": 45.7, "formula_units_per_cell": 2,
        "n_prior_DFT_points": 81, "n_new_DFT_points": 208, "n_total_DFT_points": 289,
        "prior_9x9_prospective_max_abs_error_meV_per_GaN": max_error,
        "prior_9x9_prospective_rms_error_meV_per_GaN": rms_error,
        "prior_9x9_prospective_gate": {
            "max_abs_limit_meV_per_GaN": MAX_PROSPECTIVE_ERROR_MEV_PER_GAN,
            "rms_limit_meV_per_GaN": RMS_PROSPECTIVE_ERROR_MEV_PER_GAN,
            "pass": max_error <= MAX_PROSPECTIVE_ERROR_MEV_PER_GAN
                    and rms_error <= RMS_PROSPECTIVE_ERROR_MEV_PER_GAN,
        },
        "worst_predicted_point": {
            "case": worst["case"], "q_u_A": worst["q_u_A"], "q_v_A": worst["q_v_A"],
            "signed_error_meV_per_GaN": worst["prospective_error_meV_per_GaN"],
        },
        "cases": new,
        "source_sha256": {
            "submitted_manifest": sha256(manifest_path),
            "pilot_audit": sha256(pilot_path),
            "refinement5_audit": sha256(refined5_path),
            "dense9_audit": sha256(dense9_path),
            "auditor": sha256(Path(__file__)),
        },
        "claim_limit": "Frozen local atom-strain enthalpy cut only; not a whole-path surface or certified transition state.",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("work", "pilot-audit", "refinement5-audit", "dense9-audit", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = audit(args.work, args.pilot_audit, args.refinement5_audit,
                   args.dense9_audit)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in (
        "status", "n_total_DFT_points", "prior_9x9_prospective_max_abs_error_meV_per_GaN",
        "prior_9x9_prospective_rms_error_meV_per_GaN", "prior_9x9_prospective_gate",
    )}))


if __name__ == "__main__":
    main()

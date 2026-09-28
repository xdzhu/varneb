"""Raw-audit the 56 new GaN statics and prospectively score the old 5x5 cut."""

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


def audit(work: Path, pilot_path: Path, refined5_path: Path) -> dict:
    manifest_path = work / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    pilot = json.loads(pilot_path.read_text(encoding="utf-8"))
    refined5 = json.loads(refined5_path.read_text(encoding="utf-8"))
    if (manifest.get("status") != "inputs_finalized_no_DFT"
            or manifest.get("purpose") != "GaN_45p7_600eV_local_joint_mode_9x9_frozen_refinement"
            or manifest.get("pressure_GPa") != 45.7
            or manifest.get("formula_units_per_cell") != 2
            or manifest.get("input_contract_sha256") != PRODUCTION_INPUT_SHA256
            or manifest.get("n_prior_coordinates") != 25
            or manifest.get("n_new_coordinates") != 56
            or len(manifest.get("cases", [])) != 56
            or manifest.get("source_sha256", {}).get("pilot_audit") != sha256(pilot_path)
            or manifest.get("source_sha256", {}).get("refinement5_audit") != sha256(refined5_path)
            or pilot.get("status") != "GaN_600eV_local_2D_frozen_enthalpy_pilot_raw_audited"
            or refined5.get("status") != "GaN_600eV_local_joint_5x5_refinement_raw_audited"):
        raise ValueError("9x9 manifest or preceding 25-point audit is inconsistent")
    prior = source_rows(pilot, sha256(pilot_path), refined5, sha256(refined5_path))
    xy = np.asarray([(row["q_u_A"], row["q_v_A"]) for row in prior], dtype=float)
    old_energy = np.asarray([row["measured_delta_H_meV_per_GaN"] for row in prior], dtype=float)
    predictor = CloughTocher2DInterpolator(xy, old_energy)
    h0 = float(pilot["center_enthalpy_eV_per_cell"])
    new = []
    fresh = set()
    for item in manifest["cases"]:
        q_u, q_v = float(item["q_u_A"]), float(item["q_v_A"])
        key = coordinate_key(q_u, q_v)
        if (key in fresh or key in {coordinate_key(*point) for point in xy}
                or any(item["input_sha256"].get(name) != digest
                       for name, digest in PRODUCTION_INPUT_SHA256.items())):
            raise ValueError("duplicate 9x9 node or changed VASP electronic contract")
        fresh.add(key)
        header = (work / "cases" / item["name"] / "OUTCAR").read_text(
            encoding="utf-8", errors="replace")[:1200]
        if "vasp.6.3.2" not in header or not re.search(r"running on\s+32 total cores", header):
            raise ValueError(f"unexpected VASP version or MPI width at {item['name']}")
        atoms, detail = completed_case(
            work / "cases" / item["name"],
            {**item, "POSCAR_sha256": item["input_sha256"]["POSCAR"]},
        )
        h = float(atoms.get_potential_energy() + 45.7 * GPa * atoms.get_volume())
        measured = float((h - h0) * 500)
        predicted = float(predictor(q_u, q_v))
        if not np.isfinite([h, measured, predicted]).all():
            raise ValueError(f"nonfinite 5x5 prediction at {item['name']}")
        new.append({**detail,
                    "grid_index": item["grid_index"], "q_u_A": q_u, "q_v_A": q_v,
                    "enthalpy_eV_per_cell": h,
                    "delta_enthalpy_meV_per_GaN": measured,
                    "prior_5x5_interpolation_meV_per_GaN": predicted,
                    "prospective_error_meV_per_GaN": predicted - measured})
    full = {coordinate_key(float(x), float(y))
            for x in manifest["q_u_grid_A"] for y in manifest["q_v_grid_A"]}
    if (len(full) != 81 or len(fresh) != 56
            or fresh | {coordinate_key(*point) for point in xy} != full):
        raise ValueError("9x9 grid not complete or not disjoint from old points")
    errors = np.asarray([item["prospective_error_meV_per_GaN"] for item in new])
    return {
        "status": "GaN_600eV_local_joint_9x9_refinement_raw_audited",
        "pressure_GPa": 45.7, "formula_units_per_cell": 2,
        "n_prior_DFT_points": 25, "n_new_DFT_points": 56, "n_total_DFT_points": 81,
        "prior_5x5_prospective_max_abs_error_meV_per_GaN": float(np.max(np.abs(errors))),
        "prior_5x5_prospective_rms_error_meV_per_GaN": float(np.sqrt(np.mean(errors**2))),
        "cases": new,
        "source_sha256": {
            "submitted_manifest": sha256(manifest_path),
            "pilot_audit": sha256(pilot_path),
            "refinement5_audit": sha256(refined5_path),
            "auditor": sha256(Path(__file__)),
        },
        "claim_limit": "Frozen local atom-strain enthalpy cut only; not a whole-path surface or certified transition state.",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("work", "pilot-audit", "refinement5-audit", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = audit(args.work, args.pilot_audit, args.refinement5_audit)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in (
        "status", "n_total_DFT_points", "prior_5x5_prospective_max_abs_error_meV_per_GaN",
    )}))


if __name__ == "__main__":
    main()

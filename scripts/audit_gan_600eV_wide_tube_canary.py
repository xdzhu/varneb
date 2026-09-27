"""Audit wider frozen GaN transverse probes at the unchanged 600-eV contract."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from ase.io import read
from ase.units import GPa

from scripts.audit_gan_joint_curvature import completed_case
from scripts.prepare_gan_600eV_ts_newton_canary import PRODUCTION_INPUT_SHA256
from scripts.prepare_gan_600eV_wide_tube_canary import ANCHORS, AMPLITUDE_A


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit(work: Path, trajectory: Path, feasibility_path: Path,
          first_audit_path: Path, refined_audit_path: Path,
          hessian_npz: Path) -> dict:
    manifest_path = work / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    feasibility = json.loads(feasibility_path.read_text(encoding="utf-8"))
    first = json.loads(first_audit_path.read_text(encoding="utf-8"))
    refined = json.loads(refined_audit_path.read_text(encoding="utf-8"))
    if (manifest.get("status") != "inputs_finalized_no_DFT"
            or manifest.get("purpose") != "GaN_45p7_600eV_path_tube_wide_amplitude_canary"
            or manifest.get("anchors") != list(ANCHORS)
            or manifest.get("q_perp_A") != [-AMPLITUDE_A, AMPLITUDE_A]
            or manifest.get("input_contract_sha256") != PRODUCTION_INPUT_SHA256
            or len(manifest.get("cases", [])) != 2 * len(ANCHORS)
            or manifest["source_sha256"].get("trajectory") != sha256(trajectory)
            or manifest["source_sha256"].get("feasibility") != sha256(feasibility_path)
            or manifest["source_sha256"].get("hessian_npz") != sha256(hessian_npz)
            or manifest["source_sha256"].get("preparer")
            != sha256(Path(__file__).with_name("prepare_gan_600eV_wide_tube_canary.py"))
            or first.get("status")
            != "GaN_600eV_path_adapted_frozen_transverse_tube_raw_audited"
            or refined.get("status") != "GaN_600eV_path_tube_refined_raw_audited"
            or refined["source_sha256"].get("first_audit") != sha256(first_audit_path)
            or feasibility["source_sha256"].get("trajectory") != sha256(trajectory)):
        raise ValueError("wide-canary source or unchanged electronic contract failed")
    frames = read(trajectory, index=":")
    pressure = 45.7 * GPa
    h0 = feasibility["path_enthalpy_eV_per_cell"]
    if (len(frames) != 29 or any(
            not np.isclose(float(frame.get_potential_energy() + pressure * frame.get_volume()),
                           h0[i], atol=1e-8, rtol=0)
            for i, frame in enumerate(frames))):
        raise ValueError("archived q=0 GaN path changed")
    narrow: dict[int, dict[str, float]] = {index: {} for index in ANCHORS}
    for record in [*first["cases"], *refined["new_cases"]]:
        index = int(record["image_index"])
        if index not in narrow:
            continue
        side = "plus" if record["q_perp_A"] > 0 else "minus"
        narrow[index][side] = float(record["delta_enthalpy_meV_per_GaN"])
    if any(set(pair) != {"minus", "plus"} for pair in narrow.values()):
        raise ValueError("narrow-tube baseline pair is missing")
    by_anchor: dict[int, dict[str, dict]] = {index: {} for index in ANCHORS}
    cases = []
    failures = []
    for record in manifest["cases"]:
        if any(record.get("input_sha256", {}).get(name) != digest
               for name, digest in PRODUCTION_INPUT_SHA256.items()):
            raise ValueError(f"{record['name']} changed electronic settings")
        directory = work / "cases" / record["name"]
        if any(sha256(directory / name) != digest
               for name, digest in record["input_sha256"].items()):
            raise ValueError(f"{record['name']} changed after input staging")
        stdout_path = directory / "vasp.stdout"
        if not stdout_path.is_file():
            raise FileNotFoundError(stdout_path)
        stdout = stdout_path.read_text(encoding="utf-8", errors="replace")
        if "Inconsistent Bravais lattice types found" in stdout:
            if ("I REFUSE TO CONTINUE WITH THIS SICK JOB" not in stdout
                    or "General timing and accounting informations" in stdout):
                raise ValueError(f"ambiguous VASP Bravais failure at {record['name']}")
            failures.append({
                "case": record["name"],
                "image_index": record["image_index"],
                "q_perp_A": record["q_perp_A"],
                "failure": "VASP_direct_reciprocal_Bravais_classification_conflict",
                "vasp_stdout_sha256": sha256(stdout_path),
                "input_sha256": record["input_sha256"],
                "direct_class": (
                    "base-centered monoclinic" if "Crystalline: base-centered monoclinic" in stdout
                    else "simple monoclinic" if "Crystalline: simple monoclinic" in stdout
                    else "unparsed"
                ),
                "reciprocal_class": (
                    "base-centered monoclinic" if "Reciprocal : base-centered monoclinic" in stdout
                    else "simple monoclinic" if "Reciprocal : simple monoclinic" in stdout
                    else "unparsed"
                ),
            })
            continue
        atoms, detail = completed_case(
            directory,
            {**record, "POSCAR_sha256": record["input_sha256"]["POSCAR"]},
        )
        index = int(record["image_index"])
        q = float(record["q_perp_A"])
        side = "plus" if q > 0 else "minus"
        h = float(atoms.get_potential_energy() + pressure * atoms.get_volume())
        detail.update({
            "image_index": index,
            "arc_fraction_s": feasibility["arc_fraction_s"][index],
            "q_perp_A": q,
            "enthalpy_eV_per_cell": h,
            "delta_enthalpy_meV_per_GaN": float((h - h0[index]) * 500),
        })
        if side in by_anchor[index]:
            raise ValueError("wide-canary signed pair duplicated")
        by_anchor[index][side] = detail
        detail["vasp_stdout_sha256"] = sha256(stdout_path)
        cases.append(detail)
    if len(cases) + len(failures) != 2 * len(ANCHORS):
        raise ValueError("wide-canary outcomes do not cover every submitted case")
    diagnostics = []
    for index in ANCHORS:
        if not by_anchor[index]:
            continue
        if set(by_anchor[index]) != {"minus", "plus"}:
            raise ValueError("wide-canary signed pair partially completed")
        small = narrow[index]
        wide = by_anchor[index]
        q_small = 0.015
        narrow_slope = (small["plus"] - small["minus"]) / (2 * q_small)
        narrow_curvature = (small["plus"] + small["minus"]) / q_small**2
        predictions = {
            side: sign * AMPLITUDE_A * narrow_slope
                  + 0.5 * AMPLITUDE_A**2 * narrow_curvature
            for side, sign in (("minus", -1), ("plus", 1))
        }
        discrepancies = {
            side: float(wide[side]["delta_enthalpy_meV_per_GaN"] - predictions[side])
            for side in ("minus", "plus")
        }
        diagnostics.append({
            "image_index": index,
            "narrow_q_A": q_small,
            "wide_q_A": AMPLITUDE_A,
            "narrow_local_quadratic_prediction_at_wide_q_meV_per_GaN": predictions,
            "observed_wide_delta_H_meV_per_GaN": {
                side: wide[side]["delta_enthalpy_meV_per_GaN"]
                for side in ("minus", "plus")
            },
            "wide_minus_narrow_quadratic_prediction_meV_per_GaN": discrepancies,
        })
    return {
        "status": (
            "GaN_600eV_wide_frozen_tube_canary_partial_Bravais_failure_raw_audited"
            if failures else "GaN_600eV_wide_frozen_tube_canary_raw_audited"
        ),
        "claim_limit": "Four 600-eV path anchors at q=+/-0.2 A; no global surface or conditional minimum claim",
        "pressure_GPa": 45.7,
        "formula_units_per_cell": 2,
        "anchors": list(ANCHORS),
        "q_perp_A": [-AMPLITUDE_A, AMPLITUDE_A],
        "n_completed_static_cases": len(cases),
        "n_failed_Bravais_cases": len(failures),
        "cases": cases,
        "failures": failures,
        "nonlinearity_diagnostics": diagnostics,
        "source_sha256": {
            "manifest": sha256(manifest_path),
            "trajectory": sha256(trajectory),
            "feasibility": sha256(feasibility_path),
            "first_audit": sha256(first_audit_path),
            "refined_audit": sha256(refined_audit_path),
            "hessian_npz": sha256(hessian_npz),
            "auditor": sha256(Path(__file__)),
        },
        "limitations": [
            "A safe empirical geometry screen is not a guarantee of smooth branch identity.",
            "The wider points are frozen single-point statics, not constrained minima.",
            "A narrow quadratic extrapolation to 0.2 A is diagnostic only, not an interpolation model.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("work", "trajectory", "feasibility", "first-audit",
                 "refined-audit", "hessian-npz", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    report = audit(args.work, args.trajectory, args.feasibility,
                   args.first_audit, args.refined_audit, args.hessian_npz)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": report["status"],
        "cases": report["n_completed_static_cases"],
        "Bravais_failures": report["n_failed_Bravais_cases"],
        "observed_meV_per_GaN": {
            str(row["image_index"]): row["observed_wide_delta_H_meV_per_GaN"]
            for row in report["nonlinearity_diagnostics"]
        },
    }))


if __name__ == "__main__":
    main()

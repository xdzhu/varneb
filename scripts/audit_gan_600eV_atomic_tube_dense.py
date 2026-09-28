"""Raw-audit and prospectively validate the dense GaN central (s, q) grid."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from ase.io import read
from ase.units import GPa

from scripts.audit_gan_joint_curvature import completed_case
from scripts.prepare_gan_600eV_atomic_tube_dense import (
    CENTRAL_IMAGES,
    PROSPECTIVE_HOLDOUT_IMAGES,
    TRANSVERSE_Q_A,
    missing_points,
)
from scripts.prepare_gan_600eV_ts_newton_canary import PRODUCTION_INPUT_SHA256


Q_HOLDOUT_GATE_MEV_PER_GAN = 1.0


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def qkey(q: float) -> float:
    return round(float(q), 6)


def validation_errors(arc: list[float], delta: dict[int, dict[float, float]]) -> tuple[list[dict], list[dict]]:
    """Test withheld s anchors and interior q values before fitting a contour."""
    width = 0.05

    def quadratic_from_outer(index: int, q: float) -> float:
        plus = delta[index][width]
        minus = delta[index][-width]
        slope = (plus - minus) / (2 * width)
        curvature_half = (plus + minus) / (2 * width**2)
        return slope * q + curvature_half * q**2

    s_errors = []
    training = [index for index in CENTRAL_IMAGES if index not in PROSPECTIVE_HOLDOUT_IMAGES]
    for index in PROSPECTIVE_HOLDOUT_IMAGES:
        lower = max(candidate for candidate in training if candidate < index)
        upper = min(candidate for candidate in training if candidate > index)
        fraction = (arc[index] - arc[lower]) / (arc[upper] - arc[lower])
        for q in TRANSVERSE_Q_A:
            if q == 0.0:
                continue
            predicted = ((1 - fraction) * quadratic_from_outer(lower, q)
                         + fraction * quadratic_from_outer(upper, q))
            observed = delta[index][q]
            s_errors.append({
                "image_index": index,
                "q_atom_A": q,
                "bracketing_images": [lower, upper],
                "observed_excess_meV_per_GaN": observed,
                "predicted_excess_meV_per_GaN": predicted,
                "absolute_error_meV_per_GaN": abs(observed - predicted),
            })
    q_errors = []
    for index in CENTRAL_IMAGES:
        for q in (-0.025, 0.025):
            predicted = quadratic_from_outer(index, q)
            observed = delta[index][q]
            q_errors.append({
                "image_index": index,
                "q_atom_A": q,
                "observed_excess_meV_per_GaN": observed,
                "quadratic_from_outer_pair_meV_per_GaN": predicted,
                "absolute_error_meV_per_GaN": abs(observed - predicted),
            })
    return s_errors, q_errors


def audit(work: Path, trajectory: Path, first_path: Path, refinement_path: Path) -> dict:
    manifest_path = work / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    first = json.loads(first_path.read_text(encoding="utf-8"))
    refinement = json.loads(refinement_path.read_text(encoding="utf-8"))
    expected = missing_points(first, refinement)
    if (manifest.get("status") != "inputs_finalized_no_DFT"
            or manifest.get("purpose") != "GaN_45p7_600eV_atomic_only_central_dense_2D_grid"
            or manifest.get("central_image_indices") != list(CENTRAL_IMAGES)
            or manifest.get("q_atom_A") != list(TRANSVERSE_Q_A)
            or manifest.get("prospective_s_holdout_images") != list(PROSPECTIVE_HOLDOUT_IMAGES)
            or manifest.get("predeclared_s_holdout_gate_meV_per_GaN") != 1.0
            or manifest.get("input_contract_sha256") != PRODUCTION_INPUT_SHA256
            or manifest.get("n_new_statics") != 44
            or [(int(case["image_index"]), qkey(case["q_atom_A"]))
                for case in manifest.get("cases", [])] != expected
            or manifest["source_sha256"].get("trajectory") != sha256(trajectory)
            or manifest["source_sha256"].get("first_audit") != sha256(first_path)
            or manifest["source_sha256"].get("refinement_audit") != sha256(refinement_path)
            or manifest["source_sha256"].get("preparer")
            != sha256(Path(__file__).with_name("prepare_gan_600eV_atomic_tube_dense.py"))):
        raise ValueError("dense-grid provenance or original 600-eV contract failed")
    frames = read(trajectory, index=":")
    arc = refinement["arc_fraction_s"]
    h0 = refinement["path_enthalpy_eV_per_cell"]
    pressure = 45.7 * GPa
    if (len(frames) != 29 or len(arc) != 29 or len(h0) != 29
            or any(not np.isclose(
                float(frame.get_potential_energy() + pressure * frame.get_volume()),
                h0[index], atol=1e-8, rtol=0,
            ) for index, frame in enumerate(frames))):
        raise ValueError("audited 600-eV path centers changed")
    delta: dict[int, dict[float, float]] = {index: {0.0: 0.0} for index in CENTRAL_IMAGES}
    source_records = first["cases"] + refinement["new_cases"]
    for record in source_records:
        index = int(record["image_index"])
        q = qkey(record["q_atom_A"])
        if q in delta[index]:
            raise ValueError(f"duplicate existing image {index}, q={q}")
        delta[index][q] = float(record["delta_enthalpy_meV_per_GaN"])
    new_records = []
    for record in manifest["cases"]:
        if any(record.get("input_sha256", {}).get(name) != digest
               for name, digest in PRODUCTION_INPUT_SHA256.items()):
            raise ValueError(f"electronic setting changed in {record['name']}")
        atoms, detail = completed_case(
            work / "cases" / record["name"],
            {**record, "POSCAR_sha256": record["input_sha256"]["POSCAR"]},
        )
        index = int(record["image_index"])
        q = qkey(record["q_atom_A"])
        if (q in delta[index]
                or not np.allclose(atoms.cell.array, frames[index].cell.array,
                                   atol=2e-5, rtol=0)):
            raise ValueError(f"duplicate point or altered image cell: {record['name']}")
        h = float(atoms.get_potential_energy() + pressure * atoms.get_volume())
        excess = float((h - h0[index]) * 500)
        delta[index][q] = excess
        detail.update({
            "image_index": index,
            "q_atom_A": q,
            "role": record["role"],
            "arc_fraction_s": arc[index],
            "enthalpy_eV_per_cell": h,
            "excess_enthalpy_meV_per_GaN": excess,
        })
        new_records.append(detail)
    if any(set(delta[index]) != set(TRANSVERSE_Q_A) for index in CENTRAL_IMAGES):
        raise ValueError("dense grid is incomplete")
    s_errors, q_errors = validation_errors(arc, delta)
    max_s = max(item["absolute_error_meV_per_GaN"] for item in s_errors)
    max_q = max(item["absolute_error_meV_per_GaN"] for item in q_errors)
    s_gate = max_s <= manifest["predeclared_s_holdout_gate_meV_per_GaN"]
    q_gate = max_q <= Q_HOLDOUT_GATE_MEV_PER_GAN
    return {
        "status": "GaN_600eV_central_dense_atomic_tube_raw_audited",
        "claim_limit": "Frozen central (s,q_atom) cut, not whole-path two-phonon PES or certified TS",
        "pressure_GPa": 45.7,
        "formula_units_per_cell": 2,
        "central_image_indices": list(CENTRAL_IMAGES),
        "arc_fraction_s": [float(arc[index]) for index in CENTRAL_IMAGES],
        "path_enthalpy_eV_per_cell": [float(h0[index]) for index in CENTRAL_IMAGES],
        "q_atom_A": list(TRANSVERSE_Q_A),
        "excess_enthalpy_meV_per_GaN": [
            [delta[index][q] for q in TRANSVERSE_Q_A] for index in CENTRAL_IMAGES
        ],
        "n_reused_offpath_statics": len(source_records),
        "n_reused_path_centers": len(CENTRAL_IMAGES),
        "n_new_raw_audited_statics": len(new_records),
        "prospective_s_holdout_errors": s_errors,
        "maximum_prospective_s_holdout_error_meV_per_GaN": max_s,
        "predeclared_s_holdout_gate_meV_per_GaN": 1.0,
        "s_holdout_gate_pass": s_gate,
        "q_halfstep_errors": q_errors,
        "maximum_q_halfstep_error_meV_per_GaN": max_q,
        "q_halfstep_gate_fixed_before_output_audit_meV_per_GaN": Q_HOLDOUT_GATE_MEV_PER_GAN,
        "q_halfstep_gate_pass": q_gate,
        "central_contour_gate_pass": bool(s_gate and q_gate),
        "new_cases": new_records,
        "source_sha256": {
            "manifest": sha256(manifest_path),
            "trajectory": sha256(trajectory),
            "first_audit": sha256(first_path),
            "refinement_audit": sha256(refinement_path),
            "auditor": sha256(Path(__file__)),
        },
        "limitations": [
            "The s holdouts test linear interpolation of off-path excess at two preselected images only.",
            "The q half-step gate tests a quadratic from q=0 and ±0.05 Å, not all anharmonicity.",
            "The chart is atomic-only transverse and excludes endpoint neighborhoods with mode-identity seams.",
            "The same-600-eV strain energy/stress-gradient inconsistency remains a TS-certification limit.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("work", "trajectory", "first-audit", "refinement-audit", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    report = audit(args.work, args.trajectory, args.first_audit, args.refinement_audit)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": report["status"],
        "max_s_holdout_meV_per_GaN": report["maximum_prospective_s_holdout_error_meV_per_GaN"],
        "max_q_halfstep_meV_per_GaN": report["maximum_q_halfstep_error_meV_per_GaN"],
        "contour_gate_pass": report["central_contour_gate_pass"],
    }))


if __name__ == "__main__":
    main()

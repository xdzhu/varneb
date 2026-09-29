"""Raw-audit 208 new ABACUS statics against the old 81-point BTO surface."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import numpy as np
from ase.units import GPa
from scipy.interpolate import RectBivariateSpline

from scripts.audit_bto_transverse_soft_missing22 import _number_rows


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def audit(points_root: Path, work_root: Path, old_manifest_path: Path,
          old_audit_path: Path) -> dict:
    plan_path = points_root / "plan.json"
    plan, old_manifest, old_audit = map(load, (
        plan_path, old_manifest_path, old_audit_path,
    ))
    if (plan.get("status") != "208_new_points_prepared_no_DFT"
            or plan.get("full_grid_shape") != [17, 17]
            or plan.get("old_measured_points_reused") != 81
            or plan.get("new_points") != 208
            or plan.get("source_sha256", {}).get("old_manifest") != sha256(old_manifest_path)
            or plan.get("source_sha256", {}).get("old_audit") != sha256(old_audit_path)
            or old_manifest.get("n_points") != 81
            or old_audit.get("n_total_DFT_points") != 81
            or old_audit.get("source_sha256", {}).get("dense_manifest")
               != sha256(old_manifest_path)):
        raise ValueError("BTO 17x17 source provenance changed")
    q1, q2 = np.asarray(plan["q1"], dtype=float), np.asarray(plan["q2"], dtype=float)
    if (q1.shape != (17,) or q2.shape != (17,)
            or not np.allclose(np.diff(q1), 0.15, atol=1e-12, rtol=0)
            or not np.allclose(np.diff(q2), 0.075, atol=1e-12, rtol=0)):
        raise ValueError("nested coordinate spacing changed")
    old = {(int(row["dense_i_q1"]), int(row["dense_j_q2"])): row
           for row in old_audit["samples"]}
    if len(old) != 81:
        raise ValueError("old 9x9 source has duplicate nodes")
    old_energy = np.empty((9, 9))
    combined = {}
    for (i, j), row in old.items():
        old_energy[i, j] = float(row["energy_minus_C_eV_per_BTO"])
        combined[(2 * i, 2 * j)] = {
            "grid_index_q1_q2": [2 * i, 2 * j], "q1": q1[2 * i], "q2": q2[2 * j],
            "energy_minus_C_eV_per_BTO": old_energy[i, j], "source": "reused_old_raw_audit",
        }
    predictor = RectBivariateSpline(q1[::2], q2[::2], old_energy, kx=3, ky=3, s=0)
    seen = set(combined)
    samples = []
    contract = None
    for shard in range(8):
        point_dir = points_root / f"shard-{shard:02d}"
        manifest_path = point_dir / "manifest.json"
        work = work_root / f"run-shard-{shard:02d}-001"
        result_path = work / "result_manifest.json"
        manifest, result = load(manifest_path), load(result_path)
        if (sha256(manifest_path) != plan["shard_manifest_sha256"][point_dir.name]
                or manifest.get("n_points") != 26
                or manifest.get("full_grid_shape") != [17, 17]
                or manifest.get("source_sha256") != plan["source_sha256"]
                or result.get("status") != "converged"
                or result.get("dft_requested") is not True
                or result.get("mpi_ranks") != 32
                or result.get("source_manifest_sha256") != sha256(manifest_path)
                or len(result.get("points", [])) != 26):
            raise ValueError(f"BTO shard {shard} incomplete or changed")
        parameters = result["calculator_parameters"]
        if (parameters.get("calculation") != "scf"
                or parameters.get("basis_type") != "lcao"
                or parameters.get("dft_functional") != "pbe"
                or parameters.get("ecutwfc") != 100.0
                or parameters.get("kpts") != [4, 4, 4]
                or parameters.get("cal_force") != 1
                or parameters.get("cal_stress") != 1
                or "Orb-DZP-10au" not in parameters.get("basis_dir", "")):
            raise ValueError(f"BTO shard {shard} changed the original electronic contract")
        identity = (parameters, result["asset_sha256"], result["abacus_binary_sha256"])
        if contract is None:
            contract = identity
        elif identity != contract:
            raise ValueError("BTO shards used different ABACUS contracts")
        records = {row["name"]: row for row in result["points"]}
        if len(records) != 26:
            raise ValueError(f"duplicate result in shard {shard}")
        for point in manifest["points"]:
            name = point["name"]
            row = records[name]
            i, j = point["grid_index_q1_q2"]
            key = (i, j)
            folder = work / name
            structure = point_dir / point["structure"]
            if (key in seen or (i % 2 == 0 and j % 2 == 0)
                    or row.get("status") != "converged"
                    or row.get("mpi_dsize") != 32
                    or sha256(structure) != point["structure_sha256"]
                    or row.get("source_structure_sha256") != point["structure_sha256"]
                    or sha256(folder / "source_POSCAR") != point["structure_sha256"]
                    or not np.isclose(point["q1"], q1[i], atol=1e-12, rtol=0)
                    or not np.isclose(point["q2"], q2[j], atol=1e-12, rtol=0)):
                raise ValueError(f"invalid new BTO node {name}")
            seen.add(key)
            if any(sha256(folder / filename) != digest
                   for filename, digest in row["input_sha256"].items()):
                raise ValueError(f"ABACUS input changed after preflight: {name}")
            log_path = folder / "OUT.ABACUS" / "running_scf.log"
            if sha256(log_path) != row["log_sha256"]:
                raise ValueError(f"ABACUS log hash mismatch: {name}")
            log = log_path.read_text(encoding="utf-8", errors="replace")
            if log.count("charge density convergence is achieved") != 1:
                raise ValueError(f"SCF did not uniquely converge: {name}")
            dsize = re.findall(r"\bDSIZE\s*=\s*(\d+)", log)
            totals = re.findall(r"!FINAL_ETOT_IS\s+([-+\d.eE]+)\s+eV", log)
            if dsize != ["32"] or len(totals) != 1:
                raise ValueError(f"missing MPI or final-energy marker: {name}")
            lines = log.splitlines()
            forces = _number_rows(lines, "TOTAL-FORCE (eV/Angstrom)", 5, 2)
            stress_kbar = _number_rows(lines, "TOTAL-STRESS (KBAR)", 3, 2)
            energy = float(totals[0])
            stress = -stress_kbar * 0.1 * GPa
            if (abs(energy - float(row["energy_eV"])) > 1e-6
                    or abs(float(np.max(np.linalg.norm(forces, axis=1)))
                           - float(row["max_atomic_force_eV_per_A"])) > 1e-6
                    or not np.allclose(stress, row["stress_eV_per_A3"], atol=1e-7, rtol=0)):
                raise ValueError(f"raw energy, force, or stress differs from ASE: {name}")
            measured = energy - float(old_audit["reference_energy_eV_per_BTO"])
            predicted = float(predictor.ev(q1[i], q2[j]))
            sample = {
                "grid_index_q1_q2": [i, j], "q1": float(q1[i]), "q2": float(q2[j]),
                "energy_minus_C_eV_per_BTO": measured,
                "old_9x9_cubic_prediction_eV_per_BTO": predicted,
                "prediction_minus_DFT_meV_per_BTO": 1000 * (predicted - measured),
                "source": "new_raw_audited_DFT", "log_sha256": row["log_sha256"],
            }
            combined[key] = sample
            samples.append(sample)
    if len(seen) != 289 or len(samples) != 208:
        raise ValueError("BTO nested 17x17 measured grid is incomplete")
    errors = np.array([row["prediction_minus_DFT_meV_per_BTO"] for row in samples])
    return {
        "status": "BTO_frozen_soft_mode_17x17_raw_audited",
        "claim_limit": "frozen cubic-cell dual-soft-mode cut; not a conditionally relaxed PES or T-to-C barrier",
        "n_reused_measured_points": 81, "n_new_raw_audited_statics": 208,
        "q1": q1.tolist(), "q2": q2.tolist(),
        "energy_minus_C_eV_per_BTO": [
            [combined[(i, j)]["energy_minus_C_eV_per_BTO"] for j in range(17)]
            for i in range(17)
        ],
        "prospective_max_abs_error_meV_per_BTO": float(np.max(np.abs(errors))),
        "prospective_rms_error_meV_per_BTO": float(np.sqrt(np.mean(errors**2))),
        "new_samples": samples,
        "source_sha256": {
            "plan": sha256(plan_path), "old_manifest": sha256(old_manifest_path),
            "old_audit": sha256(old_audit_path), "auditor": sha256(Path(__file__)),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("points-root", "work-root", "old-manifest", "old-audit", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    report = audit(args.points_root, args.work_root, args.old_manifest, args.old_audit)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: report[key] for key in (
        "status", "n_new_raw_audited_statics",
        "prospective_max_abs_error_meV_per_BTO",
    )}))


if __name__ == "__main__":
    main()

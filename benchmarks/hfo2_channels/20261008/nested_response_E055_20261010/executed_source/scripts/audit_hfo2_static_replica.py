"""Audit a same-input static replica; never edit or rerun its source image."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import sys

import numpy as np

from vcneb.abacus import _minimal_abacus_results


INPUT_FILES = (
    "INPUT", "KPT", "STRU", "Hf.upf", "O.upf",
    "Hf_gga_10au_100Ry_4s2p2d1f.orb", "O_gga_10au_100Ry_2s2p1d.orb",
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audited_results(directory: Path) -> dict:
    log = directory / "OUT.ABACUS/running_scf.log"
    body = log.read_text(encoding="utf-8", errors="replace")
    if "charge density convergence is achieved" not in body:
        raise ValueError(f"no electronic convergence evidence: {directory}")
    ranks = re.findall(r"\bDSIZE\s*=\s*(\d+)", body)
    if ranks != ["32"]:
        raise ValueError(f"expected one genuine DSIZE=32 calculation, got {ranks}")
    for name in ("abacus.out", "abacus.err"):
        path = directory / name
        if path.exists() and "PMI server not found" in path.read_text(encoding="utf-8", errors="replace"):
            raise ValueError("MPI singleton fallback detected")
    results = _minimal_abacus_results(log)
    if (results["forces"].shape != (12, 3) or results["stress"].shape != (6,)
            or not np.isfinite(results["energy"])
            or not np.isfinite(results["forces"]).all()
            or not np.isfinite(results["stress"]).all()):
        raise ValueError("nonfinite or incomplete Hf4O8 energy/forces/stress")
    return results


def audit_replica(source: Path, replica: Path) -> dict:
    hashes = {
        name: {"source": sha256(source / name), "replica": sha256(replica / name)}
        for name in INPUT_FILES
    }
    changed = [name for name, pair in hashes.items() if pair["source"] != pair["replica"]]
    if changed:
        raise ValueError(f"immutable input contract changed: {changed}")
    initial = audited_results(source)
    current = audited_results(replica)
    energy_delta = float(current["energy"] - initial["energy"])
    force_error = float(np.max(np.abs(current["forces"] - initial["forces"])))
    stress_error = float(np.max(np.abs(current["stress"] - initial["stress"]))) * 1602.176634
    passed = abs(energy_delta) <= 1e-5 and force_error <= 1e-4 and stress_error <= 0.02
    return {
        "schema_version": 1,
        "status": "passed" if passed else "review_required",
        "scope": "same-geometry same-input static reproducibility; not a pathway or TS certificate",
        "source_directory": str(source), "replica_directory": str(replica),
        "input_sha256": hashes,
        "source_log_sha256": sha256(source / "OUT.ABACUS/running_scf.log"),
        "replica_log_sha256": sha256(replica / "OUT.ABACUS/running_scf.log"),
        "energy_source_eV_cell": float(initial["energy"]),
        "energy_replica_eV_cell": float(current["energy"]),
        "energy_delta_eV_cell": energy_delta,
        "energy_delta_meV_fu": energy_delta * 1000 / 4,
        "maximum_force_difference_eV_A": force_error,
        "maximum_stress_difference_kbar": stress_error,
        "tolerances": {"energy_eV_cell": 1e-5, "force_eV_A": 1e-4, "stress_kbar": 0.02},
        "n_atoms": 12, "formula_units": 4, "real_mpi_ranks": 32,
        "results": {key: value.tolist() if isinstance(value, np.ndarray) else value
                    for key, value in current.items()},
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--replica", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        result = audit_replica(args.source, args.replica)
    except (OSError, ValueError) as error:
        result = {"status": "failed", "error": str(error), "source_directory": str(args.source)}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in ("status", "energy_delta_meV_fu") if key in result}))
    return 0 if result["status"] == "passed" else 1


if __name__ == "__main__":
    sys.exit(main())

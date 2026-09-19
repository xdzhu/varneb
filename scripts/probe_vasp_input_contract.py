"""Probe VASP initialization on exact regression structures, never a NEB run.

The default NELM=1 initialization probe is NOT a converged static result.
--full-scf separately verifies finite converged outputs and repeated input writes.
Neither mode automatically populates a production cache. Execution requires --run.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from ase.io import read
from vcneb import read_chain_trajectory
from vcneb.calculator import classify_calculator_failure
from vcneb.vasp import attach_vasp_calculators, default_vasp_command
from vcneb.vasp_contract import DEFAULT_VASP_SYMMETRY_PARAMETERS


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", required=True)
    sources = parser.add_mutually_exclusive_group(required=True)
    sources.add_argument("--structures", nargs="+")
    sources.add_argument("--trajectory", help="latest complete chain; never regenerate interpolation")
    parser.add_argument("--n-images", type=int, default=7, help="total images for --trajectory")
    parser.add_argument("--trajectory-index", type=int,
                        help="evaluate one zero-based image from the latest complete chain")
    parser.add_argument("--workdir", required=True)
    parser.add_argument("--symprec", type=float, default=DEFAULT_VASP_SYMMETRY_PARAMETERS["symprec"])
    parser.add_argument("--vasp-bin", default="vasp_std")
    parser.add_argument("--ncores", type=int, default=1)
    parser.add_argument("--timeout", type=float, default=120)
    parser.add_argument("--run", action="store_true")
    parser.add_argument("--full-scf", action="store_true", help="require converged energy/forces/stress instead of NELM=1")
    args = parser.parse_args()
    workdir = Path(args.workdir).resolve()
    if workdir.exists():
        raise FileExistsError(f"probe directory already exists: {workdir}")
    if args.trajectory:
        images = read_chain_trajectory(args.trajectory, n_images=args.n_images)
        labels = [f"{Path(args.trajectory).resolve()}#image={index}" for index in range(len(images))]
        if args.trajectory_index is not None:
            if not 0 <= args.trajectory_index < len(images):
                parser.error("--trajectory-index is outside the complete chain")
            labels = [labels[args.trajectory_index]]
            images = [images[args.trajectory_index]]
    else:
        if args.trajectory_index is not None:
            parser.error("--trajectory-index requires --trajectory")
        images = [read(path, format="vasp") for path in args.structures]
        labels = [str(Path(path).resolve()) for path in args.structures]
    command = default_vasp_command(args.ncores, args.vasp_bin)
    source_sha256 = {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
                     for name in ("vcneb/vasp.py", "vcneb/vasp_contract.py", "scripts/probe_vasp_input_contract.py")}
    overrides = {"xc": "PBE", "pp": "PBE", "isym": -1, "symprec": args.symprec}
    if not args.full_scf:
        overrides["nelm"] = 1
    attach_vasp_calculators(
        images, source_dir=args.source_dir, workdir=workdir, command=command,
        overrides=overrides,
        minimum_distance=1e-6,
    )
    # Validate and serialize every structure before allocating any MPI process.
    for image in images:
        image.calc.write_input(image)
    records = []
    for index, (image, structure) in enumerate(zip(images, labels)):
        directory = Path(image.calc.directory)
        record = {"structure": structure, "directory": str(directory)}
        if not args.run:
            record["status"] = "input_contract_checked_without_dft"
        else:
            try:
                with (directory / "vasp.out").open("w") as stdout, (directory / "vasp.err").open("w") as stderr:
                    with subprocess.Popen(command, shell=True, cwd=directory, stdout=stdout,
                                          stderr=stderr, start_new_session=os.name != "nt") as process:
                        try:
                            returncode = process.wait(timeout=args.timeout)
                        except subprocess.TimeoutExpired:
                            if os.name != "nt":
                                os.killpg(process.pid, signal.SIGKILL)
                            else:
                                process.kill()
                            process.wait()
                            raise
                oszicar = directory / "OSZICAR"
                electronic_step = oszicar.is_file() and bool(re.search(
                    r"^\s*(?:DAV|RMM|CG|DMP|SDA):", oszicar.read_text(errors="replace"), re.MULTILINE,
                ))
                category = classify_calculator_failure(
                    RuntimeError(f"VASP initialization return code {returncode}"),
                    diagnostic_paths=[directory / "vasp.out", directory / "vasp.err", directory / "OUTCAR"],
                )
                passed = returncode == 0 and electronic_step
                record.update(status="initialization_passed" if passed else "failed",
                              returncode=returncode, electronic_step_started=electronic_step,
                              failure_category=None if passed else category)
                if passed and args.full_scf:
                    image.calc.atoms = image.copy()
                    try:
                        image.calc.read_results()
                        import numpy as np
                        energy = float(image.calc.results["energy"])
                        forces = np.asarray(image.calc.results["forces"], dtype=float)
                        stress = np.asarray(image.calc.results["stress"], dtype=float)
                        finite = (np.isfinite(energy) and forces.shape == (len(image), 3)
                                  and stress.shape == (6,) and np.isfinite(forces).all()
                                  and np.isfinite(stress).all())
                        if not image.calc.converged or not finite:
                            raise RuntimeError("VASP SCF did not converge or returned nonfinite/invalid outputs")
                        # A second input write after ASE parses outputs must still
                        # match the immutable parameter contract; no second SCF.
                        image.calc.write_input(image)
                        record.update(status="converged_static_passed", energy_eV=energy,
                                      maximum_force_eV_per_A=float(np.linalg.norm(forces, axis=1).max()),
                                      repeated_input_write_passed=True)
                    except Exception as exc:
                        record.update(status="failed", error=str(exc), failure_category=classify_calculator_failure(exc))
            except subprocess.TimeoutExpired:
                record.update(status="failed", failure_category="timeout")
        records.append(record)
        print(f"[{record['status']}] image {index}: {structure}", flush=True)
    payload = {
        "mode": "converged_static_contract_probe" if args.full_scf else "initialization_probe_only_not_converged_scf",
        "nelm": None if args.full_scf else 1,
        "isym": -1, "symprec": args.symprec, "command": command, "records": records,
        "all_passed": all(record["status"] != "failed" for record in records),
        "production_cache_eligible": False,
        "source_sha256": source_sha256,
    }
    (workdir / "probe_summary.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return 0 if payload["all_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())

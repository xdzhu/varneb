"""Stage a bounded, same-contract GaN energy--stress derivative check.

Only the three volume-changing joint-coordinate directions are displaced at
0.01 and 0.005 Å. This checks the existing 0.02-Å discrepancy without
changing ENCUT, pseudopotentials, k points, or the production VCNEB chain.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path

import numpy as np
from ase.io import read, write

from scripts.prepare_gan_600eV_ts_hessian import geometry_preflight, same_geometry
from vcneb.joint_curvature import JointCurvatureCoordinates


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE = ROOT / "_varneb_work/gan_600eV_ts_newton_canary_20260928"
DEFAULT_AUDIT = ROOT / "benchmarks/numerical_integrity/gan_600eV_ts_newton_canary_20260928.json"
INPUT_SHA256 = {
    "INCAR": "83ba34d4b14ff4ea641bd74dec26e462ea1cc24b575db79a07138ccdfd552772",
    "KPOINTS": "b5215f3e7608c27f49d45e8129d3edca2cfaa3ba88b8c389d4b52604715712ee",
    "POTCAR": "f94781ce6cf9b9454c093353320c5e80a2a0aceea6fb8353b33b01ac37b95168",
}
CENTER_OUTCAR_SHA256 = "66e9ddf2bf5ee1f3f510bf71b8ec8dea2e05ae3d8be60e9593e741a87dabcfe5"
CELL_SCALE_A = 3.39827143300501
STEPS_A = (0.01, 0.005)
AXES = (12, 13, 14)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prepare(source: Path, center_audit: Path, output: Path) -> dict:
    source, center_audit, output = map(Path, (source, center_audit, output))
    if output.exists():
        raise FileExistsError(output)
    if sha256(source / "OUTCAR") != CENTER_OUTCAR_SHA256:
        raise ValueError("GaN candidate center OUTCAR changed")
    if any(sha256(source / name) != digest for name, digest in INPUT_SHA256.items()):
        raise ValueError("GaN 600-eV electronic input contract changed")
    audit = json.loads(center_audit.read_text(encoding="utf-8"))
    if (audit.get("status") != "GaN_600eV_one_step_static_audited_not_TS_certificate"
            or audit.get("OUTCAR_sha256") != CENTER_OUTCAR_SHA256
            or audit.get("input_contract", {}).get("original_input_sha256") != INPUT_SHA256):
        raise ValueError("GaN candidate center lacks its original raw audit")
    raw = (source / "OUTCAR").read_text(encoding="utf-8", errors="replace")
    if ("aborting loop because EDIFF is reached" not in raw
            or "General timing and accounting informations for this job" not in raw):
        raise ValueError("GaN candidate center has no completed static output")
    center = read(source / "OUTCAR")
    chart = JointCurvatureCoordinates(center, CELL_SCALE_A)
    if chart.size != 18 or center.get_chemical_symbols() != ["Ga", "Ga", "N", "N"]:
        raise ValueError("unexpected GaN candidate joint coordinates")

    output.mkdir(parents=True)
    cases_root = output / "cases"
    cases_root.mkdir()
    center_poscar = output / "center_POSCAR"
    write(center_poscar, center, format="vasp", direct=True, sort=False)
    if not same_geometry(center, read(center_poscar, format="vasp")):
        raise ValueError("center POSCAR roundtrip changed geometry")
    cases = []
    for step in STEPS_A:
        for axis in AXES:
            for sign, side in ((-1, "minus"), (1, "plus")):
                displacement = np.zeros(chart.size)
                displacement[axis] = sign * step
                atoms = chart.displaced(displacement)
                preflight = geometry_preflight(atoms)
                if preflight["empirical_near_symmetry_warning"]:
                    raise ValueError(f"near-symmetry input at {step}/{axis}/{side}")
                name = f"h{step:.3f}_axis{axis:02d}_{side}"
                case = cases_root / name
                case.mkdir()
                write(case / "POSCAR", atoms, format="vasp", direct=True,
                      sort=False)
                if not same_geometry(atoms, read(case / "POSCAR", format="vasp")):
                    raise ValueError(f"POSCAR roundtrip changed {name}")
                for filename in INPUT_SHA256:
                    shutil.copy2(source / filename, case / filename)
                hashes = {filename: sha256(case / filename)
                          for filename in ("POSCAR", *INPUT_SHA256)}
                (case / "sha256.inputs.json").write_text(
                    json.dumps(hashes, indent=2) + "\n", encoding="utf-8"
                )
                cases.append({"name": name, "step_A": step, "axis": axis,
                              "sign": sign, **preflight, "input_sha256": hashes})
    manifest = {
        "status": "inputs_finalized_no_DFT",
        "purpose": "GaN_600eV_normal_strain_derivative_step_test",
        "pressure_GPa": 45.7,
        "cell_scale_A": CELL_SCALE_A,
        "steps_A": list(STEPS_A), "axes": list(AXES), "n_cases": len(cases),
        "electronic_input_sha256": INPUT_SHA256,
        "center_OUTCAR_sha256": CENTER_OUTCAR_SHA256,
        "center_POSCAR_sha256": sha256(center_poscar),
        "source_sha256": {"center_audit": sha256(center_audit),
                          "preparer": sha256(Path(__file__))},
        "cases": cases,
        "limitations": [
            "Derivative consistency test only; not a new TS or VCNEB barrier.",
            "The finite plane-wave basis may change with strain at fixed ENCUT.",
        ],
    }
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--center-audit", type=Path, default=DEFAULT_AUDIT)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    manifest = prepare(args.source, args.center_audit, args.output)
    print(json.dumps({"status": manifest["status"],
                      "n_cases": manifest["n_cases"]}))


if __name__ == "__main__":
    main()

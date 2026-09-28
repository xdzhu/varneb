"""Stage a same-contract static check of one completed 600-eV GaN basin run.

Native VASP variable-cell relaxation can retain a history-dependent plane-wave
basis. The static point re-evaluates its final geometry with the original
600-eV VCNEB image contract; it is not another geometry optimization.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil

import numpy as np
from ase.io import read

from scripts.audit_gan_600eV_ts_basin_pilot import energy_triples, ga_n_coordination
from scripts.audit_gan_ts_basin_followups import residual_stress_kbar
from scripts.prepare_gan_600eV_ts_hessian import same_geometry, geometry_preflight
from scripts.prepare_gan_600eV_ts_newton_canary import PRODUCTION_INPUT_SHA256, sha256


PURPOSE = "GaN_45p7_600eV_relaxed_geometry_same_contract_static_not_basin_certificate"


def prepare(relax_case: Path, static_template: Path, output: Path) -> dict:
    if output.exists():
        raise FileExistsError(output)
    name = relax_case.name
    if name not in ("grid_um_vz", "grid_up_vz"):
        raise ValueError("expected one signed 600-eV GaN basin case")
    raw = (relax_case / "OUTCAR").read_text(encoding="utf-8", errors="replace")
    triples = energy_triples(raw)
    if ("General timing and accounting informations" not in raw
            or "PSTRESS=  457.0" not in raw
            or len(triples) < 10 or len(triples) > 100
            or raw.count("aborting loop because EDIFF is reached") != len(triples)):
        raise ValueError("native relaxation is not complete and SCF-audited")
    final = read(relax_case / "OUTCAR")
    contcar = read(relax_case / "CONTCAR", format="vasp")
    if (final.get_chemical_symbols() != ["Ga", "Ga", "N", "N"]
            or not same_geometry(final, contcar)
            or not np.isfinite(final.get_forces()).all()
            or not np.isfinite(final.get_stress(voigt=False)).all()):
        raise ValueError("native final geometry, atom order, force or stress invalid")
    if any(sha256(static_template / key) != digest
           for key, digest in PRODUCTION_INPUT_SHA256.items()):
        raise ValueError("static template differs from original 600-eV path contract")
    incar = (static_template / "INCAR").read_bytes()
    if (b"\r" in incar or b"PSTRESS" in incar
            or b" ENCUT = 600.000000\n" not in incar
            or b" IBRION = -1\n" not in incar
            or b" NSW = 0\n" not in incar
            or b" ISIF = 2\n" not in incar):
        raise ValueError("static INCAR is not a zero-ionic-step image evaluation")
    geometry = geometry_preflight(contcar)
    output.mkdir(parents=True)
    case = output / "case"
    case.mkdir()
    shutil.copy2(relax_case / "CONTCAR", case / "POSCAR")
    for key in PRODUCTION_INPUT_SHA256:
        shutil.copy2(static_template / key, case / key)
    if not same_geometry(read(case / "POSCAR", format="vasp"), final):
        raise ValueError("static POSCAR differs from last evaluated relaxation geometry")
    hashes = {key: sha256(case / key) for key in ("POSCAR", *PRODUCTION_INPUT_SHA256)}
    (case / "sha256.inputs.json").write_text(json.dumps(hashes, indent=2) + "\n", encoding="utf-8")
    record = {
        "status": "inputs_finalized_no_DFT", "purpose": PURPOSE,
        "case": name, "pressure_for_postprocessing_GPa": 45.7,
        "electronic_contract": {"ENCUT_eV": 600, "kpoints": "Gamma 8x8x6",
                                "PAW": "Ga_d+N", "ISYM": -1, "SYMPREC": 1e-4,
                                "static_input_sha256": PRODUCTION_INPUT_SHA256},
        "native_final_ionic_records": len(triples),
        "native_final_enthalpy_eV_per_cell": triples[-1][1],
        "native_final_max_force_eV_per_A": float(np.linalg.norm(final.get_forces(), axis=1).max()),
        "native_final_stress_residual_kbar": residual_stress_kbar(final.get_stress(voigt=False), 45.7),
        "native_final_GaN_coordination_2p4A": ga_n_coordination(final),
        "geometry": geometry, "input_sha256": hashes,
        "source_sha256": {"native_OUTCAR": sha256(relax_case / "OUTCAR"),
                          "native_CONTCAR": sha256(relax_case / "CONTCAR"),
                          "static_template_INCAR": sha256(static_template / "INCAR"),
                          "preparer": sha256(Path(__file__))},
        "limitation": "Same-geometry static comparison only; phase identity and full-variable-cell TS still require independent tests.",
    }
    (output / "manifest.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    return record


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ("relax-case", "static-template", "output"):
        parser.add_argument("--" + key, type=Path, required=True)
    args = parser.parse_args()
    result = prepare(args.relax_case, args.static_template, args.output)
    print(json.dumps({"status": result["status"], "case": result["case"],
                      "native_stress_residual_kbar": result["native_final_stress_residual_kbar"]}))


if __name__ == "__main__":
    main()

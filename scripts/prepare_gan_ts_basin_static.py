"""Stage a static 1000 eV GaN basin check from an audited relaxed cell.

The static DFT contract matches the near-TS and signed-mode calculations.
PSTRESS is absent in a static point; downstream analysis adds P*V exactly once.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil

from ase.io import read

from scripts.prepare_gan_ts_basin_continuation import geometry_metrics
from scripts.prepare_gan_ts_newton_probe import same_geometry, sha256


def prepare(endpoint: str, relax_root: Path, relax_audit_path: Path,
            structural_match_path: Path, static_template_root: Path,
            work_root: Path) -> dict:
    if work_root.exists():
        raise FileExistsError(work_root)
    if endpoint not in ("B4", "B1"):
        raise ValueError("only mapped GaN B4/B1 basins are supported")
    relax_audit = json.loads(relax_audit_path.read_text(encoding="utf-8"))
    match = json.loads(structural_match_path.read_text(encoding="utf-8"))
    if (not relax_audit.get("force_pass_0p02eV_per_A")
            or not relax_audit.get("stress_pass_1kbar")
            or not relax_audit.get("contcar_is_last_evaluated_geometry")
            or match.get("status")
            != "GaN_basin_candidate_endpoint_structural_comparison_only"
            or match.get("endpoint") != endpoint):
        raise ValueError("basin force, stress, final geometry or phase gate failed")
    source = relax_root / "case"
    if (sha256(source / "OUTCAR") != relax_audit["source_sha256"]["OUTCAR"]
            or sha256(source / "CONTCAR") != relax_audit["source_sha256"]["CONTCAR"]
            or match["source_sha256"]["candidate"]
            != relax_audit["source_sha256"]["CONTCAR"]
            or not same_geometry(read(source / "CONTCAR", format="vasp"),
                                 read(source / "OUTCAR"), tolerance=2e-5)):
        raise ValueError("audited relaxed source changed")
    static_source = static_template_root / "cases" / (
        "negative_mode_plus" if endpoint == "B4" else "negative_mode_minus"
    )
    incar = (static_source / "INCAR").read_bytes()
    if (b"\r" in incar or b"PSTRESS" in incar
            or b"ENCUT = 1000.000000\n" not in incar
            or b"SYMPREC = 1.00e-04\n" not in incar
            or b"ISYM = -1\n" not in incar
            or b"IBRION = -1\n" not in incar
            or b"NSW = 0\n" not in incar
            or (source / "KPOINTS").read_bytes() != (static_source / "KPOINTS").read_bytes()
            or sha256(source / "POTCAR") != sha256(static_source / "POTCAR")):
        raise ValueError("static template differs from 1000 eV GaN diagnostic contract")
    geometry = geometry_metrics(read(source / "CONTCAR", format="vasp"))
    destination = work_root / "case"
    destination.mkdir(parents=True)
    shutil.copy2(source / "CONTCAR", destination / "POSCAR")
    for filename in ("INCAR", "KPOINTS", "POTCAR"):
        shutil.copy2(static_source / filename, destination / filename)
    hashes = {filename: sha256(destination / filename)
              for filename in ("POSCAR", "INCAR", "KPOINTS", "POTCAR")}
    (destination / "sha256.inputs.json").write_text(
        json.dumps(hashes, indent=2) + "\n", encoding="utf-8",
    )
    result = {
        "purpose": "GaN_1000eV_audited_basin_static_E_plus_PV_not_TS_certificate",
        "status": "inputs_finalized_no_DFT", "endpoint": endpoint,
        "pressure_for_postprocessing_GPa": 45.7,
        "calculator": {"code": "VASP 6.3.2", "ENCUT_eV": 1000,
                       "kmesh": [8, 8, 6], "ISYM": -1, "SYMPREC": 1e-4},
        "initial_geometry_preflight": geometry,
        "input_sha256": hashes,
        "source_sha256": {"relax_manifest": sha256(relax_root / "manifest.json"),
                          "relax_audit": sha256(relax_audit_path),
                          "structural_match": sha256(structural_match_path),
                          "relaxed_CONTCAR": sha256(source / "CONTCAR"),
                          "static_template_INCAR": sha256(static_source / "INCAR"),
                          "preparer": sha256(Path(__file__))},
        "limitations": [
            "Static VASP free E must be combined with 45.7 GPa times its audited volume exactly once.",
            "A B1 or B4 basin check does not by itself certify both connections or TS index one.",
        ],
    }
    (work_root / "manifest.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8",
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--endpoint", choices=("B4", "B1"), required=True)
    for name in ("relax-root", "relax-audit", "structural-match",
                 "static-template-root", "work-root"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    result = prepare(args.endpoint, args.relax_root, args.relax_audit,
                     args.structural_match, args.static_template_root,
                     args.work_root)
    print(json.dumps({"status": result["status"], "endpoint": result["endpoint"]}))


if __name__ == "__main__":
    main()

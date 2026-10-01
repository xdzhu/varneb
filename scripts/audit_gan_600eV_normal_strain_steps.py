"""Audit fixed-600-eV GaN normal-strain finite-difference statics.

The comparison is descriptive: it does not infer a unique cause of
energy--stress disagreement or certify a variable-cell transition state.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from ase.io import read
from ase.units import GPa

from scripts.audit_gan_600eV_energy_stress_consistency import (
    basis_metadata, simpson_average,
)
from scripts.audit_gan_600eV_basis_history import parse_outcar
from scripts.prepare_gan_600eV_normal_strain_steps import (
    AXES, CENTER_OUTCAR_SHA256, INPUT_SHA256, STEPS_A,
)
from vcneb.joint_curvature import JointCurvatureCoordinates


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_WORK = ROOT / "paper/VARNEB_CPC/evidence/gan_600eV_normal_strain_steps_20261001"
DEFAULT_CENTER = DEFAULT_WORK / "center_OUTCAR"
DEFAULT_PRIOR = ROOT / "benchmarks/numerical_integrity/gan_600eV_energy_stress_consistency_20260928.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit(work: Path = DEFAULT_WORK, center_outcar: Path = DEFAULT_CENTER,
          prior_report: Path = DEFAULT_PRIOR) -> dict:
    work, center_outcar, prior_report = map(Path, (work, center_outcar, prior_report))
    manifest_path = work / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if (manifest.get("status") != "inputs_finalized_no_DFT"
            or manifest.get("purpose") != "GaN_600eV_normal_strain_derivative_step_test"
            or manifest.get("pressure_GPa") != 45.7
            or manifest.get("steps_A") != list(STEPS_A)
            or manifest.get("axes") != list(AXES)
            or manifest.get("n_cases") != 12
            or manifest.get("electronic_input_sha256") != INPUT_SHA256
            or manifest.get("center_OUTCAR_sha256") != CENTER_OUTCAR_SHA256
            or sha256(center_outcar) != CENTER_OUTCAR_SHA256):
        raise ValueError("normal-strain source or 600-eV contract changed")
    center_poscar = work / "center_POSCAR"
    if sha256(center_poscar) != manifest.get("center_POSCAR_sha256"):
        raise ValueError("normal-strain center POSCAR changed")
    prior = json.loads(prior_report.read_text(encoding="utf-8"))
    if (prior.get("status")
            != "GaN_600eV_raw_energy_stress_discrepancy_localized_to_volume_strain"
            or prior.get("step_A") != 0.02
            or prior.get("pressure_GPa") != 45.7):
        raise ValueError("prior 0.02-A gradient comparison changed")
    old = {record["axis"]: record for record in prior["axes"] if record["axis"] in AXES}
    if set(old) != set(AXES):
        raise ValueError("prior 0.02-A strain axes incomplete")

    center = read(center_outcar)
    chart = JointCurvatureCoordinates(read(center_poscar, format="vasp"),
                                      manifest["cell_scale_A"])
    pressure = 45.7 * GPa

    def gradient(atoms):
        return chart.enthalpy_gradient(
            atoms, atoms.get_forces(), atoms.get_stress(voigt=False), pressure
        )

    center_gradient = gradient(center)
    center_basis = basis_metadata(center_outcar)
    center_counts = parse_outcar(
        center_outcar.read_text(encoding="utf-8", errors="replace")
    )["kpoint_plane_wave_counts"]
    if len(center_counts) != 384:
        raise ValueError("candidate center has an unexpected number of k points")
    cases = {}
    for item in manifest["cases"]:
        key = (item["step_A"], item["axis"], item["sign"])
        if key in cases or key[0] not in STEPS_A or key[1] not in AXES or key[2] not in (-1, 1):
            raise ValueError("duplicate or unexpected normal-strain case")
        directory = work / "cases" / item["name"]
        expected = item["input_sha256"]
        if (set(expected) != {"INCAR", "KPOINTS", "POSCAR", "POTCAR"}
                or json.loads((directory / "sha256.inputs.json").read_text()) != expected
                or any(sha256(directory / name) != digest for name, digest in expected.items()
                       if name != "POTCAR")
                or any(expected.get(name) != digest for name, digest in INPUT_SHA256.items())):
            raise ValueError(f"normal-strain input changed: {item['name']}")
        # POTCAR is license-controlled and intentionally omitted from the
        # public evidence archive. The Slurm worker verified its bytes.
        if (directory / "POTCAR").exists() and sha256(directory / "POTCAR") != expected["POTCAR"]:
            raise ValueError(f"normal-strain POTCAR changed: {item['name']}")
        outcar = directory / "OUTCAR"
        raw = outcar.read_text(encoding="utf-8", errors="replace")
        if ("aborting loop because EDIFF is reached" not in raw
                or "General timing and accounting informations for this job" not in raw):
            raise ValueError(f"normal-strain static incomplete: {item['name']}")
        counts = parse_outcar(raw)["kpoint_plane_wave_counts"]
        if len(counts) != len(center_counts):
            raise ValueError(f"normal-strain k-point count changed: {item['name']}")
        atoms = read(outcar)
        input_atoms = read(directory / "POSCAR", format="vasp")
        if (atoms.get_chemical_symbols() != ["Ga", "Ga", "N", "N"]
                or np.max(np.abs(atoms.cell.array - input_atoms.cell.array)) > 1e-5
                or np.max(np.abs(atoms.positions - input_atoms.positions)) > 1e-4):
            raise ValueError(f"normal-strain evaluated geometry changed: {item['name']}")
        cases[key] = {
            "atoms": atoms, "gradient": gradient(atoms),
            "basis": basis_metadata(outcar), "OUTCAR_sha256": sha256(outcar),
            "input_sha256": expected, "kpoint_plane_wave_counts": counts,
        }
    if len(cases) != 12:
        raise ValueError("expected all 12 independently evaluated statics")

    rows = []
    for step in STEPS_A:
        for axis in AXES:
            minus, plus = cases[(step, axis, -1)], cases[(step, axis, 1)]
            a, b = minus["atoms"], plus["atoms"]
            h_a = a.get_potential_energy() + pressure * a.get_volume()
            h_b = b.get_potential_energy() + pressure * b.get_volume()
            secant = (h_b - h_a) / (2 * step)
            integrated = simpson_average(
                minus["gradient"][axis], center_gradient[axis],
                plus["gradient"][axis],
            )
            volume_secant = (b.get_volume() - a.get_volume()) / (2 * step)
            residual = secant - integrated
            minus_counts = minus["kpoint_plane_wave_counts"]
            plus_counts = plus["kpoint_plane_wave_counts"]
            rows.append({
                "step_A": step, "axis": axis,
                "energy_enthalpy_secant_eV_per_A": float(secant),
                "simpson_force_stress_gradient_eV_per_A": float(integrated),
                "secant_minus_simpson_eV_per_A": float(residual),
                "observed_minus_integrated_pressure_GPa": float(
                    -residual / volume_secant / GPa
                ),
                "maximum_plane_waves": {
                    "minus": minus["basis"]["maximum_plane_waves"],
                    "center": center_basis["maximum_plane_waves"],
                    "plus": plus["basis"]["maximum_plane_waves"],
                },
                "kpoint_plane_wave_counts": {
                    "n_kpoints": len(center_counts),
                    "sum_minus_center_plus": [
                        sum(minus_counts), sum(center_counts), sum(plus_counts),
                    ],
                    "n_changed_minus_vs_center": sum(
                        a != b for a, b in zip(minus_counts, center_counts)
                    ),
                    "n_changed_plus_vs_center": sum(
                        a != b for a, b in zip(plus_counts, center_counts)
                    ),
                    "n_changed_minus_vs_plus": sum(
                        a != b for a, b in zip(minus_counts, plus_counts)
                    ),
                },
                "fft_grids_unchanged": all(
                    case["basis"][key] == center_basis[key]
                    for case in (minus, plus)
                    for key in ("coarse_fft_grid", "fine_fft_grid")
                ),
                "OUTCAR_sha256": {"minus": minus["OUTCAR_sha256"],
                                  "plus": plus["OUTCAR_sha256"]},
            })
    comparison = []
    for axis in AXES:
        series = [old[axis]["secant_minus_simpson_eV_per_A"]]
        series += [next(row["secant_minus_simpson_eV_per_A"] for row in rows
                        if row["axis"] == axis and row["step_A"] == step)
                   for step in STEPS_A]
        comparison.append({"axis": axis, "steps_A": [0.02, *STEPS_A],
                           "residual_eV_per_A": series})
    return {
        "status": "GaN_600eV_normal_strain_step_dependence_raw_audited",
        "purpose": manifest["purpose"],
        "pressure_GPa": 45.7,
        "electronic_input_sha256": INPUT_SHA256,
        "n_new_statics": len(cases),
        "POTCAR_locally_reverified": all(
            (work / "cases" / item["name"] / "POTCAR").exists()
            for item in manifest["cases"]
        ),
        "rows": rows, "comparison_to_prior_0p02_A": comparison,
        "limitations": [
            "A step-size trend does not uniquely identify finite-basis/Pulay effects.",
            "This is a derivative-consistency check, not a stationary TS or new barrier.",
            "POTCAR is not redistributed; the archived hash was verified on hf at each worker launch.",
        ],
        "source_sha256": {
            "manifest": sha256(manifest_path),
            "center_OUTCAR": sha256(center_outcar),
            "prior_report": sha256(prior_report),
            "auditor": sha256(Path(__file__)),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work", type=Path, default=DEFAULT_WORK)
    parser.add_argument("--center-outcar", type=Path, default=DEFAULT_CENTER)
    parser.add_argument("--prior-report", type=Path, default=DEFAULT_PRIOR)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = audit(args.work, args.center_outcar, args.prior_report)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": result["status"],
                      "comparison": result["comparison_to_prior_0p02_A"]}))


if __name__ == "__main__":
    main()

"""Audit GaN 600-eV energy/force/stress derivatives without changing inputs.

The 0.02-A finite-difference Hessian uses 36 independent VASP statics.  This
script compares their enthalpy secants with Simpson-integrated generalized
gradients, which separates ordinary finite-step curvature from a genuine
energy-versus-gradient inconsistency.  It never runs VASP.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import numpy as np
from ase.io import read
from ase.units import GPa

from vcneb.joint_curvature import JointCurvatureCoordinates


ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "_varneb_work/gan_600eV_ts_hessian_20260928"
CENTER = ROOT / "_varneb_work/gan_600eV_ts_newton_canary_20260928"
HESSIAN_AUDIT = ROOT / "benchmarks/numerical_integrity/gan_600eV_ts_hessian_0p02_20260928/audit.json"
OUTPUT = ROOT / "benchmarks/numerical_integrity/gan_600eV_energy_stress_consistency_20260928.json"
INPUT_SHA256 = {
    "INCAR": "83ba34d4b14ff4ea641bd74dec26e462ea1cc24b575db79a07138ccdfd552772",
    "KPOINTS": "b5215f3e7608c27f49d45e8129d3edca2cfaa3ba88b8c389d4b52604715712ee",
    "POTCAR": "f94781ce6cf9b9454c093353320c5e80a2a0aceea6fb8353b33b01ac37b95168",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def basis_metadata(outcar: Path) -> dict:
    content = outcar.read_text(encoding="utf-8", errors="replace")
    count = re.search(r"maximum number of plane-waves:\s*(\d+)", content)
    coarse = re.search(
        r"dimension x,y,z NGX\s*=\s*(\d+) NGY\s*=\s*(\d+) NGZ\s*=\s*(\d+)",
        content,
    )
    fine = re.search(
        r"dimension x,y,z NGXF\s*=\s*(\d+) NGYF\s*=\s*(\d+) NGZF\s*=\s*(\d+)",
        content,
    )
    if not (count and coarse and fine):
        raise ValueError(f"missing plane-wave or FFT metadata: {outcar}")
    return {
        "maximum_plane_waves": int(count.group(1)),
        "coarse_fft_grid": [int(value) for value in coarse.groups()],
        "fine_fft_grid": [int(value) for value in fine.groups()],
    }


def simpson_average(minus: float, center: float, plus: float) -> float:
    """Mean derivative over a symmetric interval from three samples."""

    return (minus + 4.0 * center + plus) / 6.0


def audit(work: Path, center_work: Path, prior_audit: Path) -> dict:
    manifest_path = work / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    previous = json.loads(prior_audit.read_text(encoding="utf-8"))
    center_outcar = center_work / "OUTCAR"
    if (
        manifest.get("purpose") != "GaN_600eV_candidate_joint_hessian_0p02_A"
        or manifest.get("step_A") != 0.02
        or manifest.get("pressure_GPa") != 45.7
        or manifest.get("source_input_sha256") != INPUT_SHA256
        or manifest.get("source_sha256", {}).get("center_OUTCAR") != sha256(center_outcar)
        or previous.get("source_sha256", {}).get("manifest") != sha256(manifest_path)
        or previous.get("source_sha256", {}).get("center_OUTCAR") != sha256(center_outcar)
        or previous.get("n_static_displacements") != 36
    ):
        raise ValueError("GaN 600-eV source or prior raw audit changed")
    for name, expected in INPUT_SHA256.items():
        if sha256(center_work / name) != expected:
            raise ValueError(f"center {name} changed from the 600-eV contract")
    center_poscar = work / "center_POSCAR"
    if sha256(center_poscar) != manifest.get("center_POSCAR_sha256"):
        raise ValueError("Hessian reference structure changed")

    center = read(center_outcar)
    chart = JointCurvatureCoordinates(read(center_poscar, format="vasp"), manifest["cell_scale_A"])
    pressure = manifest["pressure_GPa"] * GPa
    step = manifest["step_A"]

    def gradient(atoms):
        return chart.enthalpy_gradient(
            atoms, atoms.get_forces(), atoms.get_stress(voigt=False), pressure
        )

    center_gradient = gradient(center)
    center_basis = basis_metadata(center_outcar)
    audited_cases = {record["case"]: record for record in previous["cases"]}
    rows = []
    for axis in range(chart.size):
        signed = []
        for side in ("minus", "plus"):
            name = f"axis{axis:02d}_{side}"
            directory = work / "cases" / name
            record = audited_cases[name]
            staged = manifest["cases"][2 * axis + (0 if side == "plus" else 1)]
            if staged["name"] != name:
                raise ValueError(f"manifest displacement order changed at {name}")
            expected_inputs = staged["input_sha256"]
            if expected_inputs != record["input_sha256"]:
                raise ValueError(f"input audit disagrees for {name}")
            for filename, expected in expected_inputs.items():
                if sha256(directory / filename) != expected:
                    raise ValueError(f"{name}/{filename} changed after input audit")
            outcar = directory / "OUTCAR"
            if sha256(outcar) != record["outcar_sha256"]:
                raise ValueError(f"{name}/OUTCAR changed after raw audit")
            atoms = read(outcar)
            signed.append((atoms, gradient(atoms), basis_metadata(outcar)))
        (minus, gm, bm), (plus, gp, bp) = signed
        h_minus = minus.get_potential_energy() + pressure * minus.get_volume()
        h_plus = plus.get_potential_energy() + pressure * plus.get_volume()
        numerical = (h_plus - h_minus) / (2.0 * step)
        integrated = simpson_average(gm[axis], center_gradient[axis], gp[axis])
        group = "atomic" if axis < 12 else "diagonal_strain" if axis < 15 else "shear_strain"
        row = {
            "axis": axis,
            "group": group,
            "energy_enthalpy_secant_eV_per_A": float(numerical),
            "center_force_stress_gradient_eV_per_A": float(center_gradient[axis]),
            "simpson_force_stress_gradient_eV_per_A": float(integrated),
            "secant_minus_simpson_eV_per_A": float(numerical - integrated),
            "volume_secant_A2": float((plus.get_volume() - minus.get_volume()) / (2.0 * step)),
            "maximum_plane_waves": {
                "minus": bm["maximum_plane_waves"],
                "center": center_basis["maximum_plane_waves"],
                "plus": bp["maximum_plane_waves"],
            },
            "fft_grids_unchanged": (
                bm["coarse_fft_grid"] == center_basis["coarse_fft_grid"] == bp["coarse_fft_grid"]
                and bm["fine_fft_grid"] == center_basis["fine_fft_grid"] == bp["fine_fft_grid"]
            ),
        }
        if 12 <= axis < 15:
            component = axis - 12
            volume_per_q = center.get_volume() / chart.cell_scale_A
            for atoms, derived in ((minus, gm), (center, center_gradient), (plus, gp)):
                independent = (
                    atoms.get_stress(voigt=False)[component, component] + pressure
                ) * volume_per_q
                if not np.isclose(independent, derived[axis], atol=1e-9, rtol=0):
                    raise ValueError(f"diagonal strain {axis} stress transform failed")
            row["observed_minus_integrated_pressure_GPa"] = float(
                -(numerical - integrated) / row["volume_secant_A2"] / GPa
            )
        rows.append(row)

    max_residual = {
        group: max(abs(row["secant_minus_simpson_eV_per_A"]) for row in rows if row["group"] == group)
        for group in ("atomic", "diagonal_strain", "shear_strain")
    }
    if not all(row["fft_grids_unchanged"] for row in rows):
        raise ValueError("FFT mesh changed within a displacement pair")
    return {
        "status": "GaN_600eV_raw_energy_stress_discrepancy_localized_to_volume_strain",
        "electronic_contract": "VASP/PBE/Ga_d+N/600 eV/Gamma 8x8x6; original INCAR/KPOINTS/POTCAR hashes",
        "pressure_GPa": manifest["pressure_GPa"],
        "step_A": step,
        "method": "central enthalpy secant minus Simpson-integrated generalized gradient",
        "units": {
            "gradient": "eV/Angstrom",
            "volume_secant": "Angstrom^2",
            "pressure_difference": "GPa",
        },
        "max_abs_residual_eV_per_A_by_group": max_residual,
        "center_basis": center_basis,
        "axes": rows,
        "interpretation": (
            "The mismatch is concentrated in volume-changing normal strains; atomic and "
            "volume-conserving shear directions are much closer after finite-step integration. "
            "The plane-wave count changes while the FFT grids stay fixed. This is compatible "
            "with finite-basis/Pulay discontinuities, but the present data do not uniquely prove "
            "that mechanism or certify a transition state. Do not change ENCUT or mix protocols."
        ),
        "source_sha256": {
            "manifest": sha256(manifest_path),
            "prior_hessian_audit": sha256(prior_audit),
            "center_POSCAR": sha256(center_poscar),
            "center_OUTCAR": sha256(center_outcar),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work", type=Path, default=WORK)
    parser.add_argument("--center-work", type=Path, default=CENTER)
    parser.add_argument("--prior-audit", type=Path, default=HESSIAN_AUDIT)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = audit(args.work, args.center_work, args.prior_audit)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": result["status"], "max_residuals": result["max_abs_residual_eV_per_A_by_group"]}))


if __name__ == "__main__":
    main()

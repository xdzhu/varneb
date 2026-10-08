"""BFGS a prepared endpoint with a fixed substrate and audited open traction.

The factory remains calculator-independent. The seed manifest prescribes the
mechanical ensemble; it does not override any calculator parameter. A fresh
workdir is mandatory, and no free-cell evaluation is reused as a seed result.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import tempfile

import numpy as np
from ase.io import read, write

from scripts.audit_hfo2_static_replica import sha256
from scripts.prepare_hfo2_clamped_endpoints import geometry_guard
from scripts.relax_ase_endpoint import EndpointBFGS, _json_object, _load_symbol
from vcneb import (ClampedPlaneFilter, attach_image_calculators,
                   clamped_plane_vcneb_boundary, endpoint_structure_record,
                   validate_image_calculators)


EV_A3_TO_KBAR = 1602.176634
GPA_PER_EV_A3 = EV_A3_TO_KBAR / 10


def stress_report(atoms, boundary, pressure_gpa=0.0):
    stress = np.asarray(atoms.get_stress(voigt=False), dtype=float)
    pressure = pressure_gpa / GPA_PER_EV_A3
    open_traction = boundary.open_traction(stress, pressure=pressure)
    full_traction = (stress + pressure * np.eye(3)) @ boundary.normal
    normal_stress = float(np.dot(full_traction, boundary.normal))
    return {
        "stress_eV_per_A3": stress.tolist(),
        "open_traction_kbar": (open_traction * EV_A3_TO_KBAR).tolist(),
        "open_traction_norm_kbar": float(np.linalg.norm(open_traction) * EV_A3_TO_KBAR),
        "normal_residual_kbar": normal_stress * EV_A3_TO_KBAR,
        "shear_traction_norm_kbar": float(np.linalg.norm(
            full_traction - normal_stress * boundary.normal) * EV_A3_TO_KBAR),
        "max_abs_full_stress_kbar_information_only": float(np.abs(stress).max() * EV_A3_TO_KBAR),
    }


class ClampedEndpointBFGS(EndpointBFGS):
    """Use physical atomic forces and only the boundary's released stresses."""

    def __init__(self, *args, boundary, stress_kbar, **kwargs):
        if not np.isfinite(stress_kbar) or stress_kbar <= 0:
            raise ValueError("a finite positive open-stress threshold is required")
        self.boundary = boundary
        super().__init__(*args, stress_kbar=stress_kbar, **kwargs)

    def _physical_converged(self):
        self.boundary.validate_images([self.endpoint_atoms])
        forces = np.asarray(self.endpoint_atoms.get_forces(), dtype=float)
        if forces.shape != (len(self.endpoint_atoms), 3) or not np.all(np.isfinite(forces)):
            raise ValueError("finite atomic forces required")
        return bool(np.linalg.norm(forces, axis=1).max() < self.fmax
                    and stress_report(self.endpoint_atoms, self.boundary, self.pressure_gpa)
                    ["open_traction_norm_kbar"] < self.stress_kbar)


def load_seed(manifest_path):
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest["format_version"] != 1 or manifest["status"] != "geometry_seed_not_relaxed":
        raise ValueError("prepared, unrelaxed endpoint seed manifest required")
    filename = manifest["seed_file"]
    if not isinstance(filename, str) or Path(filename).name != filename or "/" in filename or "\\" in filename:
        raise ValueError("seed_file must name a sibling file, not another namespace")
    structure = manifest_path.parent / filename
    if sha256(structure) != manifest["seed_sha256"]:
        raise ValueError("seed geometry checksum changed")
    atoms = read(structure, format="vasp")
    if atoms.get_chemical_symbols() != manifest["ordered_species"]:
        raise ValueError("ordered seed identity changed")
    if atoms.constraints:
        raise ValueError("all atomic freedoms must be released")
    boundary = clamped_plane_vcneb_boundary(len(atoms), np.array(manifest["reference_cell_A"]),
                                            allow_tilt=manifest["allow_tilt"])
    boundary.validate_images([atoms])
    geometry_guard([atoms], minimum_distance_A=manifest["minimum_distance_guard_A"])
    if not np.isfinite(manifest["cell_scale_A"]) or manifest["cell_scale_A"] <= 0:
        raise ValueError("manifest cell_scale_A must be finite and positive")
    if not np.isfinite(manifest["pressure_gpa"]):
        raise ValueError("manifest pressure must be finite")
    return atoms, boundary, manifest


def write_json_atomic(path, value):
    descriptor, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(value, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def relax(atoms, boundary, manifest, workdir, *, factory, fmax=.03, stress_kbar=2.,
          steps=80, maxstep=.02, seed_manifest_sha256=None):
    """Fresh bounded endpoint optimization; preserves every raw calculator call."""
    if (not np.isfinite(fmax) or fmax <= 0 or not np.isfinite(maxstep) or maxstep <= 0
            or not isinstance(steps, int) or isinstance(steps, bool) or steps < 1
            or not np.isfinite(stress_kbar) or stress_kbar <= 0):
        raise ValueError("positive finite optimizer gates and an integer step cap required")
    boundary.validate_images([atoms])

    def guard(images):
        geometry_guard(images, minimum_distance_A=manifest["minimum_distance_guard_A"])

    guard([atoms])
    workdir.mkdir(parents=True, exist_ok=False)
    attach_image_calculators([atoms], workdir=workdir / "calculator", factory=factory)
    validate_image_calculators([atoms], require_stress=True, require_variable_cell=True,
                              require_directory=True, require_unique_directories=True)
    target = ClampedPlaneFilter(atoms, boundary, cell_scale_A=manifest["cell_scale_A"],
                               pressure_eV_per_A3=manifest["pressure_gpa"] / GPA_PER_EV_A3,
                               candidate_validator=guard)
    optimizer = ClampedEndpointBFGS(target, endpoint_atoms=atoms, boundary=boundary,
                                    pressure_gpa=manifest["pressure_gpa"], stress_kbar=stress_kbar,
                                    maxstep=maxstep, logfile=str(workdir / "relax.log"),
                                    trajectory=str(workdir / "relax.traj"))

    def snapshot():
        boundary.validate_images([atoms])
        energy = float(atoms.get_potential_energy())
        forces = np.asarray(atoms.get_forces(), dtype=float)
        generalized = np.asarray(target.get_forces(), dtype=float)
        if not np.isfinite(energy) or not np.isfinite(forces).all() or not np.isfinite(generalized).all():
            raise ValueError("nonfinite endpoint evaluation")
        return {
            "optimizer": "BFGS(ClampedPlaneFilter)", "optimizer_steps": optimizer.nsteps,
            "phase_label_seed_not_certified": manifest["phase_label"],
            "strain": manifest["strain"], "external_pressure_gpa": manifest["pressure_gpa"],
            "potential_energy_eV": energy,
            "enthalpy_eV": energy + manifest["pressure_gpa"] / GPA_PER_EV_A3 * atoms.get_volume(),
            "max_atomic_force_eV_per_A": float(np.linalg.norm(forces, axis=1).max()),
            "max_generalized_force_eV_per_A": float(np.linalg.norm(generalized, axis=1).max()),
            "fmax_target_eV_per_A": fmax, "open_stress_target_kbar": stress_kbar,
            "steps_requested": steps, "maxstep": maxstep, "cell_scale_A": manifest["cell_scale_A"],
            "cell_A": atoms.cell.array.tolist(), "volume_A3": float(atoms.get_volume()),
            "substrate_plane_drift_A": float(np.max(np.abs(atoms.cell[:2] - boundary.reference_cell[:2]))),
            "allow_tilt": boundary.allow_tilt, "endpoint": endpoint_structure_record(atoms),
            "seed_manifest_sha256": seed_manifest_sha256,
            "physical_contract_sha256": manifest.get("physical_contract_sha256"),
            "phase_and_variant_gate": "pending_independent_audit_no_symmetry_forcing",
            **stress_report(atoms, boundary, manifest["pressure_gpa"]),
        }

    def checkpoint():
        report = {"status": "running", **snapshot()}
        write(workdir / "CONTCAR.current", atoms, format="vasp", direct=True, vasp5=True, sort=False)
        write_json_atomic(workdir / "endpoint_relax_state.json", report)

    optimizer.attach(checkpoint)
    converged = bool(optimizer.run(fmax=fmax, steps=steps))
    report = {"status": "completed" if converged else "step_limit", "converged": converged, **snapshot()}
    write(workdir / "CONTCAR", atoms, format="vasp", direct=True, vasp5=True, sort=False)
    write_json_atomic(workdir / "endpoint_relax_summary.json", report)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed-manifest", type=Path, required=True)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--factory", required=True, help="module:attribute ASE factory")
    parser.add_argument("--parameters", required=True)
    parser.add_argument("--command", default=None)
    parser.add_argument("--fmax", type=float, default=.03)
    parser.add_argument("--stress-kbar", type=float, default=2.)
    parser.add_argument("--steps", type=int, default=80)
    parser.add_argument("--maxstep", type=float, default=.02)
    args = parser.parse_args()
    atoms, boundary, manifest = load_seed(args.seed_manifest)
    kwargs = {"parameters": _json_object(args.parameters)}
    if args.command is not None:
        kwargs["command"] = args.command
    factory = _load_symbol(args.factory)(**kwargs)
    report = relax(atoms, boundary, manifest, args.workdir, factory=factory,
                   fmax=args.fmax, stress_kbar=args.stress_kbar, steps=args.steps,
                   maxstep=args.maxstep, seed_manifest_sha256=sha256(args.seed_manifest))
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if report["converged"] else 1


if __name__ == "__main__":
    raise SystemExit(main())

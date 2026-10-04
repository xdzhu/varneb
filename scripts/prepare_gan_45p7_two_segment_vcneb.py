"""Prepare the two fixed-center GaN VCNEB halves without running DFT.

The 29-image production chain and three completed 600-eV statics are immutable
inputs.  The one-step-refined image 15 is a *candidate*, not a certified TS.
Licensed POTCAR bytes are copied only within the remote calculation directory.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil

import numpy as np
from ase.io import read, write
from ase.units import GPa

from scripts.audit_gan_45p7_vasp_endpoints import (
    EXPECTED_CHAIN_SHA256,
    EXPECTED_SOURCE_SHA256,
)
from scripts.prepare_gan_600eV_ts_newton_canary import PRODUCTION_INPUT_SHA256
from vcneb import endpoint_structure_record, validate_path_geometry
from vcneb.vasp import prepare_vasp_static_parameters, vasp_input_fingerprints


PRESSURE_GPA = 45.7
CENTER_INDEX = 15
TOTAL_IMAGES = 29
EXPECTED_CENTER_OUTCAR_SHA256 = (
    "66e9ddf2bf5ee1f3f510bf71b8ec8dea2e05ae3d8be60e9593e741a87dabcfe5"
)
EXPECTED_CENTER_ENTHALPY_EV_PER_CELL = -11.496211487167525


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _same_periodic_geometry(a, b, *, cell_tolerance_A=2e-5,
                            position_tolerance_A=2e-5) -> bool:
    if a.get_chemical_symbols() != b.get_chemical_symbols():
        return False
    if not np.array_equal(a.pbc, b.pbc):
        return False
    if np.max(np.abs(a.cell.array - b.cell.array)) > cell_tolerance_A:
        return False
    delta = a.get_scaled_positions(wrap=False) - b.get_scaled_positions(wrap=False)
    delta[:, a.pbc] -= np.rint(delta[:, a.pbc])
    return bool(np.max(np.linalg.norm(delta @ b.cell.array, axis=1)) <= position_tolerance_A)


def split_frames(frames, center):
    """Preserve the original ordering and make one shared fixed center."""
    if len(frames) != TOTAL_IMAGES:
        raise ValueError("expected exactly 29 archived production images")
    if center.get_chemical_symbols() != frames[CENTER_INDEX].get_chemical_symbols():
        raise ValueError("candidate atom order differs from image 15")
    if not np.array_equal(center.pbc, frames[CENTER_INDEX].pbc):
        raise ValueError("candidate PBC differs from image 15")
    local_center = center.copy()
    scaled = local_center.get_scaled_positions(wrap=False)
    reference = frames[CENTER_INDEX].get_scaled_positions(wrap=False)
    scaled[:, local_center.pbc] += np.rint((reference - scaled)[:, local_center.pbc])
    local_center.set_scaled_positions(scaled)
    if np.max(np.abs(local_center.cell.array - frames[CENTER_INDEX].cell.array)) > 0.03:
        raise ValueError("candidate cell is not a local image-15 refinement")
    delta = (local_center.get_scaled_positions(wrap=False)
             - frames[CENTER_INDEX].get_scaled_positions(wrap=False))
    if np.max(np.linalg.norm(delta @ local_center.cell.array, axis=1)) > 0.03:
        raise ValueError("candidate atoms are not a local image-15 refinement")
    left = [item.copy() for item in frames[:CENTER_INDEX + 1]]
    right = [item.copy() for item in frames[CENTER_INDEX:]]
    left[-1] = local_center.copy()
    right[0] = local_center.copy()
    return left, right


def _raw_static(directory: Path, expected_outcar_hash: str):
    for name, digest in PRODUCTION_INPUT_SHA256.items():
        if sha256(directory / name) != digest:
            raise ValueError(f"{directory}/{name} differs from the 600-eV input contract")
    outcar = directory / "OUTCAR"
    if sha256(outcar) != expected_outcar_hash:
        raise ValueError(f"raw static hash changed: {outcar}")
    content = outcar.read_text(encoding="utf-8", errors="replace")
    if ("aborting loop because EDIFF is reached" not in content
            or "General timing and accounting informations for this job" not in content):
        raise ValueError(f"incomplete or unconverged electronic static: {outcar}")
    return read(outcar)


def _static_summary(raw, target, *, role: str, n_images: int,
                    parameters: dict, fingerprints: dict, outcar: Path) -> dict:
    if not _same_periodic_geometry(raw, target):
        raise ValueError(f"{role} cached static does not match the fixed endpoint")
    forces = np.asarray(raw.get_forces(), dtype=float)
    stress = np.asarray(raw.get_stress(voigt=True), dtype=float)
    energy = float(raw.get_potential_energy())
    if (forces.shape != (len(target), 3) or stress.shape != (6,)
            or not np.isfinite(energy) or not np.isfinite(forces).all()
            or not np.isfinite(stress).all()):
        raise ValueError(f"{role} cached static lacks complete energy/force/stress")
    return {
        "status": "completed",
        "execution_mode": f"fixed_{role}_endpoint_static_scf",
        "evaluated_image_index": 0 if role == "initial" else n_images - 1,
        "n_images": n_images,
        "endpoint_structures": {role: endpoint_structure_record(target)},
        "licensed_input_fingerprints": fingerprints,
        "calculator_parameters": parameters,
        "potential_energy_eV": energy,
        "forces_eV_per_A": forces.tolist(),
        "stress_eV_per_A3_voigt": stress.tolist(),
        "max_force_eV_per_A": float(np.linalg.norm(forces, axis=1).max()),
        "raw_OUTCAR_sha256": sha256(outcar),
        "provenance": "completed same-contract raw static, repackaged for fixed segment endpoint",
    }


def prepare(*, trajectory: Path, center_static: Path, initial_static: Path,
            final_static: Path, input_dir: Path, output: Path) -> dict:
    if output.exists():
        raise FileExistsError(output)
    if sha256(trajectory) != EXPECTED_CHAIN_SHA256:
        raise ValueError("archived 29-image production chain changed")
    for name in ("INCAR", "KPOINTS", "POTCAR"):
        if not (input_dir / name).is_file():
            raise FileNotFoundError(input_dir / name)
    if sha256(input_dir / "POTCAR") != PRODUCTION_INPUT_SHA256["POTCAR"]:
        raise ValueError("nominal VASP input does not use the original Ga_d+N POTCAR")
    params, _ = prepare_vasp_static_parameters(
        input_dir, overrides={"xc": "PBE", "pp": "PBE", "isym": -1, "symprec": 1e-4}
    )
    if (float(params.get("encut", 0)) != 600 or list(params.get("kpts", ())) != [8, 8, 6]
            or params.get("isym") != -1 or float(params.get("symprec", 0)) != 1e-4
            or float(params.get("ediff", 0)) != 1e-7
            or params.get("setups") != {"Ga": "_d", "N": ""}):
        raise ValueError("effective VASP static parameters changed from production")

    frames = read(trajectory, index=":")
    b4 = _raw_static(initial_static, EXPECTED_SOURCE_SHA256["initial_00_OUTCAR"])
    b1 = _raw_static(final_static, EXPECTED_SOURCE_SHA256["final_28_OUTCAR"])
    peak = _raw_static(center_static, EXPECTED_CENTER_OUTCAR_SHA256)
    if (not _same_periodic_geometry(b4, frames[0])
            or not _same_periodic_geometry(b1, frames[-1])):
        raise ValueError("archived endpoint statics differ from the production chain")
    center_h = float(peak.get_potential_energy()
                     + PRESSURE_GPA * GPa * peak.get_volume())
    if abs(center_h - EXPECTED_CENTER_ENTHALPY_EV_PER_CELL) > 1e-6:
        raise ValueError("candidate enthalpy differs from the audited static")
    left, right = split_frames(frames, peak)
    for segment in (left, right):
        validate_path_geometry(segment, minimum_distance=1.4,
                               maximum_deformation=0.50)

    output.mkdir(parents=True)
    shared = output / "shared_center.vasp"
    write(shared, left[-1], format="vasp", direct=True, sort=False, vasp5=True)
    canonical_center = read(shared, format="vasp")
    if not _same_periodic_geometry(canonical_center, peak):
        raise ValueError("shared center POSCAR roundtrip changed geometry")
    left[-1] = canonical_center.copy()
    right[0] = canonical_center.copy()
    sources = {"B4": (b4, initial_static / "OUTCAR"),
               "center": (peak, center_static / "OUTCAR"),
               "B1": (b1, final_static / "OUTCAR")}
    manifest = {
        "status": "two_segment_preflight_passed_no_DFT",
        "pressure_GPa": PRESSURE_GPA,
        "fixed_center": "one-step-refined image 15, not a certified stationary TS",
        "electronic_contract": "VASP/PBE/Ga_d+N/600 eV/Gamma 8x8x6/ISYM=-1/SYMPREC=1e-4",
        "source_chain_sha256": EXPECTED_CHAIN_SHA256,
        "center_OUTCAR_sha256": EXPECTED_CENTER_OUTCAR_SHA256,
        "shared_center_POSCAR_sha256": sha256(shared),
        "segments": {},
    }
    for name, indices, images, endpoint_names in (
        ("left", [0, CENTER_INDEX], left, ("B4", "center")),
        ("right", [CENTER_INDEX, TOTAL_IMAGES - 1], right, ("center", "B1")),
    ):
        root = output / name
        start_dir, stop_dir = root / "initial", root / "final"
        start_dir.mkdir(parents=True)
        stop_dir.mkdir()
        write(start_dir / "CONTCAR", images[0], format="vasp", direct=True,
              sort=False, vasp5=True)
        write(stop_dir / "CONTCAR", images[-1], format="vasp", direct=True,
              sort=False, vasp5=True)
        for input_name in ("INCAR", "KPOINTS", "POTCAR"):
            shutil.copy2(input_dir / input_name, start_dir / input_name)
        images[0] = read(start_dir / "CONTCAR", format="vasp")
        images[-1] = read(stop_dir / "CONTCAR", format="vasp")
        if not (_same_periodic_geometry(images[0], sources[endpoint_names[0]][0])
                and _same_periodic_geometry(images[-1], sources[endpoint_names[1]][0])):
            raise ValueError(f"{name} endpoint POSCAR roundtrip changed a static geometry")
        write(root / "seed.traj", images)
        fingerprints = vasp_input_fingerprints(start_dir)
        for role, image, source_name in (
            ("initial", images[0], endpoint_names[0]),
            ("final", images[-1], endpoint_names[1]),
        ):
            raw, outcar = sources[source_name]
            payload = _static_summary(raw, image, role=role, n_images=len(images),
                                      parameters=params, fingerprints=fingerprints,
                                      outcar=outcar)
            (root / f"static_{role}.json").write_text(
                json.dumps(payload, indent=2) + "\n", encoding="utf-8"
            )
        manifest["segments"][name] = {
            "source_image_indices_inclusive": indices,
            "total_images": len(images),
            "active_interior_images": len(images) - 2,
            "seed_trajectory_sha256": sha256(root / "seed.traj"),
            "endpoint_names": endpoint_names,
            "endpoint_structures": {
                "initial": endpoint_structure_record(images[0]),
                "final": endpoint_structure_record(images[-1]),
            },
        }
    if (manifest["segments"]["left"]["endpoint_structures"]["final"]["sha256"]
            != manifest["segments"]["right"]["endpoint_structures"]["initial"]["sha256"]):
        raise ValueError("two segments do not share the same exact center")
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n",
                                            encoding="utf-8")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("trajectory", "center-static", "initial-static", "final-static",
                 "input-dir", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    report = prepare(trajectory=args.trajectory, center_static=args.center_static,
                     initial_static=args.initial_static, final_static=args.final_static,
                     input_dir=args.input_dir, output=args.output)
    print(json.dumps({"status": report["status"],
                      "segments": {key: value["active_interior_images"]
                                   for key, value in report["segments"].items()}}))


if __name__ == "__main__":
    main()

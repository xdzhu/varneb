"""Run VC-NEB with VASP through ASE.

This script is intentionally close to ``run_NEB/run_ase_cineb.sh`` but allows
the cell to vary along the band.  It reads relaxed endpoint structures and
copies INCAR/KPOINTS/POTCAR settings from the initial endpoint directory.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np
from ase.io import read, write

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from vcneb import (
    Mode,
    VCNEB,
    build_mode_basis,
    interpolate_vcneb,
    endpoint_structure_record,
    path_geometry_diagnostics,
    mode_guided_path,
    read_chain_trajectory,
    run_vcneb,
    validate_image_calculators,
)
from vcneb.executor import ThreadedCalculatorExecutor
from vcneb.vasp import (
    attach_vasp_calculators,
    default_vasp_command,
    prepare_vasp_static_parameters,
    vasp_input_fingerprints,
)


def read_endpoint(directory: Path):
    for name in ["CONTCAR", "CONTCAR.vasp", "POSCAR"]:
        path = directory / name
        if path.exists():
            atoms = read(path)
            print(f"[OK] read {path}")
            return atoms
    raise FileNotFoundError(f"No endpoint structure found in {directory}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    root = Path(__file__).resolve().parents[1]
    parser.add_argument("--initial", default=str(root / "initial_state" / "relax"))
    parser.add_argument("--final", default=str(root / "final_state" / "relax"))
    parser.add_argument("--workdir", default=str(root / "run_VCNEB" / "run"))
    parser.add_argument("--n-images", type=int, default=7, help="Total images including endpoints")
    parser.add_argument("--fmax", type=float, default=0.10, help="NEB convergence threshold in eV/A")
    parser.add_argument("--steps", type=int, default=300)
    parser.add_argument("--k", type=float, default=0.10)
    parser.add_argument("--pressure-gpa", type=float, default=0.0)
    parser.add_argument(
        "--optimizer",
        choices=["FIRE", "BFGS", "LBFGS", "BFGSLineSearch"],
        default="FIRE",
    )
    parser.add_argument("--line-search-retries", type=int, default=0)
    parser.add_argument("--line-search-retry-factor", type=float, default=0.5)
    parser.add_argument("--vasp-bin", default=os.environ.get("VASP_BIN", "vasp_std"))
    parser.add_argument("--ncores", type=int, default=int(os.environ.get("NP", "8")))
    parser.add_argument("--image-workers", type=int, default=0)
    parser.add_argument("--image-retries", type=int, default=0)
    parser.add_argument("--image-manifest", default=None)
    parser.add_argument("--image-cache-dir", default=None)
    parser.add_argument("--image-cache-namespace", default=None)
    parser.add_argument("--validate-only", action="store_true")
    parser.add_argument(
        "--static-only",
        action="store_true",
        help="Evaluate fixed initial endpoint 00 once, then write vasp_static_summary.json; never optimize a path",
    )
    parser.add_argument("--mic", action="store_true")
    parser.add_argument("--cell-interpolation", choices=["linear", "log_strain"], default="linear")
    parser.add_argument("--mapping", choices=["identity", "auto"], default="identity")
    parser.add_argument("--align-translation", action="store_true")
    parser.add_argument("--minimum-distance", type=float, default=None)
    parser.add_argument("--maximum-deformation", type=float, default=None)
    parser.add_argument("--mode", default=None, help="Optional JSON/NPZ/text atomic-plus-cell mode")
    parser.add_argument("--mode-guided", action="store_true", help="Apply the mode to the initial path")
    parser.add_argument("--mode-amplitude", type=float, default=0.25)
    parser.add_argument("--mode-envelope", choices=["sin", "bell", "linear"], default="sin")
    parser.add_argument(
        "--constraint-mode",
        choices=["none", "subspace", "projected"],
        default="none",
        help="Optional strict mode-subspace or projected-update constraint",
    )
    parser.add_argument("--mode-cell-scale", type=float, default=None)
    parser.add_argument("--mode-mass-weighted-input", action="store_true")
    parser.add_argument("--mode-remove-translation", action="store_true")
    parser.add_argument("--no-climb", action="store_true")
    parser.add_argument("--resume", action="store_true", help="Resume from the latest complete chain in vcneb.traj")
    parser.add_argument("--resume-trajectory", default=None, help="Trajectory to resume from; defaults to workdir/vcneb.traj")
    return parser.parse_args()


def _write_json_atomic(path: Path, payload: dict) -> None:
    def json_default(value):
        if isinstance(value, Path):
            return str(value)
        tolist = getattr(value, "tolist", None)
        if callable(tolist):
            return tolist()
        item = getattr(value, "item", None)
        if callable(item):
            return item()
        raise TypeError(f"cannot serialize {type(value).__name__} in VCNEB provenance")

    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=path.parent, delete=False) as handle:
        json.dump(payload, handle, indent=2, sort_keys=True, default=json_default)
        handle.write("\n")
        temporary = Path(handle.name)
    temporary.replace(path)


def _git_revision() -> str | None:
    declared = os.environ.get("VCNEB_GIT_REVISION")
    if declared:
        return declared
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=ROOT,
            check=True,
            capture_output=True,
            text=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return None
    return result.stdout.strip() or None


def main() -> None:
    args = parse_args()
    if args.image_workers < 0 or args.image_retries < 0:
        raise ValueError("--image-workers and --image-retries must be non-negative")
    if args.validate_only and args.static_only:
        raise ValueError("--validate-only and --static-only are mutually exclusive")
    initial_dir = Path(args.initial).resolve()
    final_dir = Path(args.final).resolve()
    workdir = Path(args.workdir).resolve()
    workdir.mkdir(parents=True, exist_ok=True)

    initial = read_endpoint(initial_dir)
    final = read_endpoint(final_dir)
    mode = Mode.from_file(args.mode, n_atoms=len(initial)) if args.mode else None
    if args.mode_guided and mode is None:
        raise ValueError("--mode-guided requires --mode")
    if args.constraint_mode != "none" and mode is None:
        raise ValueError("--constraint-mode requires --mode")
    if args.mode_mass_weighted_input and mode is None:
        raise ValueError("--mode-mass-weighted-input requires --mode")
    if args.mode_remove_translation and mode is None:
        raise ValueError("--mode-remove-translation requires --mode")
    mode_masses = (
        initial.get_masses()
        if mode is not None and (args.mode_mass_weighted_input or args.mode_remove_translation)
        else None
    )

    traj_path = workdir / "vcneb.traj"
    resume_path = Path(args.resume_trajectory).resolve() if args.resume_trajectory else traj_path
    if args.resume:
        images = read_chain_trajectory(resume_path, n_images=args.n_images)
        print(f"[OK] resumed latest complete {args.n_images}-image chain from {resume_path}")
    else:
        path_kwargs = dict(
            n_images=args.n_images,
            align_cells=True,
            mic=args.mic,
            cell_interpolation=args.cell_interpolation,
            mapping=None if args.mapping == "identity" else "auto",
            align_translation=args.align_translation,
            minimum_distance=args.minimum_distance,
            maximum_deformation=args.maximum_deformation,
        )
        if args.mode_guided:
            images = mode_guided_path(
                initial,
                final,
                mode=mode,
                amplitude=args.mode_amplitude,
                envelope=args.mode_envelope,
                masses=mode_masses,
                mass_weighted_input=args.mode_mass_weighted_input,
                remove_translation=args.mode_remove_translation,
                **path_kwargs,
            )
        else:
            images = interpolate_vcneb(initial, final, **path_kwargs)
    mode_basis = None
    if mode is not None and args.constraint_mode == "subspace":
        write(workdir / "initial-vcneb-unprojected.traj", images)
    if mode is not None and args.constraint_mode != "none":
        mode_basis = build_mode_basis(
            mode,
            initial,
            cell_scale=args.mode_cell_scale,
            masses=mode_masses,
            mass_weighted_input=args.mode_mass_weighted_input,
            remove_translation=args.mode_remove_translation,
        )
        if args.constraint_mode == "subspace":
            VCNEB(images, mode_basis=mode_basis, constraint_mode="subspace", k=args.k)
    write(workdir / "initial-vcneb.traj", images)

    command = default_vasp_command(args.ncores, args.vasp_bin)
    static_parameters, _ = prepare_vasp_static_parameters(
        initial_dir,
        overrides={"xc": "PBE", "pp": "PBE"},
    )
    attach_vasp_calculators(
        images,
        source_dir=initial_dir,
        workdir=workdir,
        command=command,
        overrides={"xc": "PBE", "pp": "PBE"},
    )
    reports = validate_image_calculators(
        images,
        require_stress=True,
        require_variable_cell=True,
        require_directory=True,
        require_unique_directories=True,
    )
    image_manifest = (
        Path(args.image_manifest).resolve()
        if args.image_manifest
        else workdir / "image_worker_manifest.jsonl"
    )
    metadata = {
        "git_revision": _git_revision(),
        "command_line": sys.argv,
        "calculator": "ASE VASP",
        "n_images": args.n_images,
        "n_interior_images": max(0, args.n_images - 2),
        "endpoint_evaluation_policy": (
            "fixed_initial_endpoint_static_scf"
            if args.static_only
            else "fixed_cached_once" if args.image_workers else "ASE_calculator_cache"
        ),
        "image_workers": args.image_workers,
        "image_retries": args.image_retries,
        "cell_interpolation": args.cell_interpolation,
        "mapping": args.mapping,
        "align_translation": args.align_translation,
        "fmax_target_eV_per_A": args.fmax,
        "calculator_parameters": static_parameters,
        "licensed_input_fingerprints": vasp_input_fingerprints(initial_dir),
        "calculator_reports": [report.to_dict() for report in reports],
        "initial_path_geometry": path_geometry_diagnostics(images),
        "endpoint_structures": {
            "initial": endpoint_structure_record(initial),
            "final": endpoint_structure_record(final),
        },
    }
    _write_json_atomic(workdir / "vcneb_preflight.json", {"status": "ok", **metadata})
    if args.validate_only:
        print(f"[OK] VASP VCNEB preflight passed; report={workdir / 'vcneb_preflight.json'}")
        return
    if args.static_only:
        # This establishes a real electronic baseline without turning the
        # fixed endpoint into an NEB worker or permitting ionic/cell updates.
        static_image = images[0]
        forces = np.asarray(static_image.get_forces(), dtype=float)
        stress = np.asarray(static_image.get_stress(), dtype=float)
        summary = {
            "status": "completed",
            **metadata,
            "execution_mode": "fixed_initial_endpoint_static_scf",
            "evaluated_image_index": 0,
            "potential_energy_eV": float(static_image.get_potential_energy()),
            "forces_eV_per_A": forces.tolist(),
            "stress_eV_per_A3_voigt": stress.tolist(),
            "max_force_eV_per_A": float(np.linalg.norm(forces, axis=1).max()),
        }
        _write_json_atomic(workdir / "vasp_static_summary.json", summary)
        print(f"[DONE] VASP fixed-endpoint static SCF; workdir={workdir}")
        return

    executor = (
        ThreadedCalculatorExecutor(
            args.image_workers,
            max_retries=args.image_retries,
            manifest_path=image_manifest,
            cache_dir=args.image_cache_dir,
            cache_namespace=args.image_cache_namespace,
        )
        if args.image_workers
        else None
    )

    chain, _ = run_vcneb(
        images,
        pressure_gpa=args.pressure_gpa,
        k=args.k,
        climb=not args.no_climb,
        mic=args.mic,
        image_executor=executor,
        mode_basis=mode_basis,
        constraint_mode=None if args.constraint_mode == "none" else args.constraint_mode,
        optimizer=args.optimizer,
        line_search_retries=args.line_search_retries,
        line_search_retry_factor=args.line_search_retry_factor,
        fmax=args.fmax,
        steps=args.steps,
        logfile=workdir / "vcneb.opt.log",
        trajectory=traj_path,
        trajectory_mode="a" if args.resume and resume_path == traj_path and traj_path.exists() else "w",
        snapshot_dir=workdir / "snapshots",
        failure_report=workdir / "vcneb_failure.json",
    )

    for image_index, image in enumerate(chain.images):
        write(workdir / f"{image_index:02d}" / "POSCAR.final", image, format="vasp", direct=True, vasp5=True)
    chain.plot_band(workdir / "vcneb_barrier.png")
    barrier, delta = chain.barrier()
    summary = {
        "status": "completed",
        **metadata,
        "barrier_enthalpy_eV": barrier,
        "reaction_enthalpy_eV": delta,
        "final_max_generalized_force_eV_per_A": chain.gradient_norm(-chain.get_forces()),
        "image_enthalpies_eV": [float(value) for value in chain.enthalpies],
        "path_diagnostics": chain.path_diagnostics(),
        "saddle_diagnostics": chain.saddle_diagnostics(),
        "image_manifest": str(image_manifest) if executor is not None else None,
    }
    _write_json_atomic(workdir / "vcneb_summary.json", summary)
    with open(workdir / "vcneb_summary.txt", "w", encoding="utf-8") as handle:
        handle.write(f"Forward barrier (enthalpy) = {barrier:.8f} eV\n")
        handle.write(f"Reaction enthalpy          = {delta:.8f} eV\n")
        handle.write("Image enthalpies (eV)      = " + " ".join(f"{e:.8f}" for e in chain.enthalpies) + "\n")
    print(f"[DONE] barrier={barrier:.6f} eV delta={delta:.6f} eV workdir={workdir}")


if __name__ == "__main__":
    main()

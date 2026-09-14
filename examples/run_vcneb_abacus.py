"""Run VC-NEB with ABACUS through a calculator factory.

The default path uses ``ase.calculators.abacus.Abacus`` when available.  If
your production environment uses a different ABACUS ASE adapter, replace
``make_ase_abacus_factory(...)`` with a custom factory that returns one
calculator per image directory.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

from ase.io import read, write

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from vcneb import (
    Mode,
    build_mode_basis,
    interpolate_vcneb,
    mode_guided_path,
    path_geometry_diagnostics,
    read_chain_trajectory,
    run_vcneb,
    validate_image_calculators,
)
from vcneb.abacus import attach_abacus_calculators, make_ase_abacus_factory
from vcneb.executor import ThreadedCalculatorExecutor


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--initial", required=True, help="ASE-readable initial endpoint")
    parser.add_argument("--final", required=True, help="ASE-readable final endpoint")
    parser.add_argument("--workdir", default="abacus_vcneb_run")
    parser.add_argument(
        "--n-images",
        type=int,
        default=7,
        help="Total images including the two fixed endpoints (7 means 5 interior images)",
    )
    parser.add_argument("--fmax", type=float, default=0.05)
    parser.add_argument("--steps", type=int, default=300)
    parser.add_argument("--k", type=float, default=0.10)
    parser.add_argument("--pressure-gpa", type=float, default=0.0)
    parser.add_argument(
        "--command",
        default=os.environ.get("ABACUS_COMMAND"),
        help="ABACUS launcher command; defaults to ABACUS_COMMAND when set",
    )
    parser.add_argument("--pseudo-dir", default=os.environ.get("ABACUS_PP_PATH"))
    parser.add_argument("--basis-dir", default=os.environ.get("ABACUS_ORBITAL_PATH"))
    parser.add_argument("--pp", action="append", default=[], metavar="SPECIES=FILE")
    parser.add_argument("--basis", action="append", default=[], metavar="SPECIES=FILE")
    parser.add_argument("--ecutwfc", type=float, default=None)
    parser.add_argument("--scf-thr", type=float, default=None)
    parser.add_argument("--scf-nmax", type=int, default=None)
    parser.add_argument("--mixing-beta", type=float, default=None)
    parser.add_argument("--kpts", type=int, nargs=3, default=None, metavar=("NX", "NY", "NZ"))
    parser.add_argument(
        "--optimizer",
        choices=["FIRE", "BFGS", "LBFGS", "BFGSLineSearch"],
        default="FIRE",
    )
    parser.add_argument(
        "--maxstep",
        type=float,
        default=None,
        help="Optional ASE optimizer maxstep in Angstrom-like extended coordinates",
    )
    parser.add_argument("--mic", action="store_true")
    parser.add_argument(
        "--cell-interpolation",
        choices=["linear", "log_strain"],
        default="linear",
        help="Cell deformation path used for initial images",
    )
    parser.add_argument(
        "--mapping",
        choices=["identity", "auto"],
        default="identity",
        help="Endpoint atom mapping strategy",
    )
    parser.add_argument(
        "--align-translation",
        action="store_true",
        help="Optimize a common periodic endpoint translation with atom mapping",
    )
    parser.add_argument(
        "--minimum-distance",
        type=float,
        default=None,
        help="Reject initial paths whose periodic atom separation is below this Angstrom threshold",
    )
    parser.add_argument(
        "--maximum-deformation",
        type=float,
        default=None,
        help="Reject initial paths whose Frobenius deformation exceeds this threshold",
    )
    parser.add_argument(
        "--mode",
        default=None,
        help="Optional JSON/NPZ/text atomic-plus-cell mode used for path guidance or constraints",
    )
    parser.add_argument(
        "--mode-guided",
        action="store_true",
        help="Apply the supplied mode only to the initial path (subsequent VC-NEB is unconstrained unless requested)",
    )
    parser.add_argument(
        "--mode-amplitude",
        type=float,
        default=0.25,
        help="Mode-guided initial-path amplitude in the extended coordinate metric",
    )
    parser.add_argument(
        "--mode-envelope",
        choices=["sin", "bell", "linear"],
        default="sin",
        help="Endpoint-zero envelope for --mode-guided",
    )
    parser.add_argument(
        "--constraint-mode",
        choices=["none", "subspace", "projected"],
        default="none",
        help="Optional strict mode-subspace or projected-update constraint",
    )
    parser.add_argument(
        "--mode-cell-scale",
        type=float,
        default=None,
        help="Cell scale used when building a mode basis; defaults to reference-cell volume root",
    )
    parser.add_argument(
        "--mode-mass-weighted-input",
        action="store_true",
        help="Interpret supplied atomic mode components as mass-weighted eigenvector values",
    )
    parser.add_argument(
        "--mode-remove-translation",
        action="store_true",
        help="Remove the (mass-weighted, when available) translational component before normalization",
    )
    parser.add_argument("--no-climb", action="store_true")
    parser.add_argument(
        "--climb-after",
        type=int,
        default=None,
        help="Enable CI after this many completed ordinary-NEB optimizer steps",
    )
    parser.add_argument("--resume", action="store_true", help="Resume from the latest complete chain in vcneb.traj")
    parser.add_argument("--resume-trajectory", default=None, help="Trajectory to resume from; defaults to workdir/vcneb.traj")
    parser.add_argument(
        "--resume-step",
        type=int,
        default=None,
        help="complete chain snapshot index to resume (negative counts from the end; default latest)",
    )
    parser.add_argument(
        "--validate-only",
        action="store_true",
        help="Build the path and calculator directories, run preflight, write a report, and skip DFT",
    )
    parser.add_argument(
        "--image-workers",
        type=int,
        default=0,
        help="concurrently evaluate this many image calculators; 0 keeps serial evaluation",
    )
    parser.add_argument(
        "--image-retries",
        type=int,
        default=0,
        help="retry a failed image calculator this many times when image-workers is enabled",
    )
    parser.add_argument(
        "--image-manifest",
        default=None,
        help="append one JSONL record per controller image-evaluation batch",
    )
    parser.add_argument(
        "--image-cache-dir",
        default=None,
        help="optional durable cache directory for exact image energy/force/stress evaluations",
    )
    parser.add_argument(
        "--image-cache-namespace",
        default=None,
        help="calculator-parameter namespace recorded with --image-cache-dir",
    )
    parser.add_argument(
        "--line-search-retries",
        type=int,
        default=0,
        help="bounded retries for ASE BFGSLineSearch failures only",
    )
    parser.add_argument(
        "--line-search-retry-factor",
        type=float,
        default=0.5,
        help="step-cap factor applied to each line-search retry",
    )
    return parser.parse_args()


def parse_species_files(values: list[str], option: str) -> dict[str, str]:
    result = {}
    for value in values:
        if "=" not in value:
            raise ValueError(f"{option} expects SPECIES=FILE, got {value!r}")
        species, filename = value.split("=", 1)
        if not species or not filename:
            raise ValueError(f"{option} expects non-empty SPECIES and FILE, got {value!r}")
        result[species] = filename
    return result


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


def _write_json_atomic(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, ensure_ascii=False)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass


def _run_metadata(args: argparse.Namespace, workdir: Path) -> dict:
    slurm_keys = (
        "SLURM_JOB_ID",
        "SLURM_JOB_NAME",
        "SLURM_JOB_NODELIST",
        "SLURM_NTASKS",
        "SLURM_CPUS_PER_TASK",
        "SLURM_MEM_PER_NODE",
        "SLURM_JOB_PARTITION",
    )
    return {
        "git_revision": _git_revision(),
        "workdir": str(workdir),
        "command_line": sys.argv,
        "slurm": {key: os.environ[key] for key in slurm_keys if os.environ.get(key)},
        "resume": bool(args.resume),
        "resume_trajectory": str(args.resume_trajectory) if args.resume_trajectory else None,
        "resume_step": args.resume_step,
        "image_workers": int(args.image_workers),
        "image_retries": int(args.image_retries),
        "image_manifest": str(args.image_manifest) if args.image_manifest else None,
        "image_cache_dir": str(args.image_cache_dir) if args.image_cache_dir else None,
        "image_cache_namespace": args.image_cache_namespace,
        "mode": str(args.mode) if args.mode else None,
        "mode_guided": bool(args.mode_guided),
        "mode_amplitude": float(args.mode_amplitude),
        "mode_envelope": args.mode_envelope,
        "constraint_mode": args.constraint_mode,
        "mode_cell_scale": args.mode_cell_scale,
        "mode_mass_weighted_input": bool(args.mode_mass_weighted_input),
        "mode_remove_translation": bool(args.mode_remove_translation),
    }


def main() -> None:
    args = parse_args()
    if args.image_workers < 0:
        raise ValueError("--image-workers must be non-negative")
    if args.image_retries < 0:
        raise ValueError("--image-retries must be non-negative")
    workdir = Path(args.workdir).resolve()
    workdir.mkdir(parents=True, exist_ok=True)

    initial = read(args.initial)
    final = read(args.final)
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
        resume_step = -1 if args.resume_step is None else args.resume_step
        images = read_chain_trajectory(resume_path, n_images=args.n_images, step=resume_step)
        label = "latest complete" if resume_step == -1 else f"complete step {resume_step}"
        print(f"[OK] resumed {label} {args.n_images}-image chain from {resume_path}")
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
    initial_trajectory = workdir / "initial-vcneb.traj"
    if not args.resume or not initial_trajectory.exists():
        write(initial_trajectory, images)
    initial_geometry = path_geometry_diagnostics(images)
    mode_basis = None
    if mode is not None and args.constraint_mode != "none":
        mode_basis = build_mode_basis(
            mode,
            initial,
            cell_scale=args.mode_cell_scale,
            masses=mode_masses,
            mass_weighted_input=args.mode_mass_weighted_input,
            remove_translation=args.mode_remove_translation,
        )

    parameters = {
        "calculation": "scf",
        "basis_type": "lcao",
        "dft_functional": "pbe",
        "cal_force": 1,
        "cal_stress": 1,
        "out_stru": 1,
    }
    if args.pseudo_dir:
        parameters["pseudo_dir"] = args.pseudo_dir
    if args.basis_dir:
        parameters["basis_dir"] = args.basis_dir
    pp = parse_species_files(args.pp, "--pp")
    basis = parse_species_files(args.basis, "--basis")
    if pp:
        parameters["pp"] = pp
    if basis:
        parameters["basis"] = basis
    optional_parameters = {
        "ecutwfc": args.ecutwfc,
        "scf_thr": args.scf_thr,
        "scf_nmax": args.scf_nmax,
        "mixing_beta": args.mixing_beta,
        "kpts": args.kpts,
    }
    parameters.update({key: value for key, value in optional_parameters.items() if value is not None})
    factory = make_ase_abacus_factory(parameters=parameters, command=args.command)
    attach_abacus_calculators(images, workdir=workdir, factory=factory)
    calculator_reports = validate_image_calculators(
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
    metadata = _run_metadata(args, workdir)
    metadata["image_manifest"] = str(image_manifest) if args.image_workers else None
    metadata["n_interior_images"] = max(0, int(args.n_images) - 2)
    metadata["endpoint_evaluation_policy"] = (
        "fixed_cached_once" if args.image_workers else "ASE_calculator_cache"
    )
    preflight = {
        **metadata,
        "status": "ok",
        "n_images": args.n_images,
        "n_interior_images": max(0, int(args.n_images) - 2),
        "initial_path_geometry": initial_geometry,
        "calculator_reports": [report.to_dict() for report in calculator_reports],
        "calculator_parameters": parameters,
    }
    _write_json_atomic(workdir / "vcneb_preflight.json", preflight)
    if args.validate_only:
        print(f"[OK] calculator preflight passed; report={workdir / 'vcneb_preflight.json'}")
        return

    image_executor = (
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
        climb_after=args.climb_after,
        mode_basis=mode_basis,
        constraint_mode=None if args.constraint_mode == "none" else args.constraint_mode,
        image_executor=image_executor,
        optimizer=args.optimizer,
        optimizer_kwargs={} if args.maxstep is None else {"maxstep": args.maxstep},
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
    max_force = chain.gradient_norm(-chain.get_forces())
    diagnostics = chain.path_diagnostics()
    saddle = chain.saddle_diagnostics()
    summary = {
        **metadata,
        "status": "completed",
        "workdir": str(workdir),
        "n_images": args.n_images,
        "n_interior_images": max(0, int(args.n_images) - 2),
        "endpoint_evaluation_policy": metadata["endpoint_evaluation_policy"],
        "cell_interpolation": args.cell_interpolation,
        "mapping": args.mapping,
        "align_translation": args.align_translation,
        "minimum_distance_threshold_A": args.minimum_distance,
        "maximum_deformation_threshold": args.maximum_deformation,
        "climb_after": args.climb_after,
        "mode": str(Path(args.mode).resolve()) if args.mode else None,
        "mode_guided": args.mode_guided,
        "mode_amplitude": args.mode_amplitude,
        "mode_envelope": args.mode_envelope,
        "constraint_mode": args.constraint_mode,
        "mode_cell_scale": args.mode_cell_scale,
        "mode_basis_columns": int(mode_basis.shape[1]) if mode_basis is not None else 0,
        "initial_path_geometry": initial_geometry,
        "initial_path_metadata": images[0].info.get("vcneb_path_metadata", {}),
        "optimizer": args.optimizer,
        "line_search_retries": args.line_search_retries,
        "line_search_retry_factor": args.line_search_retry_factor,
        "steps_requested": args.steps,
        "fmax_target_eV_per_A": args.fmax,
        "final_max_generalized_force_eV_per_A": max_force,
        "barrier_enthalpy_eV": barrier,
        "reaction_enthalpy_eV": delta,
        "image_enthalpies_eV": [float(value) for value in chain.enthalpies],
        "image_manifest": str(image_manifest) if image_executor is not None else None,
        "calculator": "ASE ABACUS",
        "calculator_parameters": parameters,
        "calculator_reports": [report.to_dict() for report in calculator_reports],
        "stress_required": True,
        "saddle_diagnostics": saddle,
        "path_diagnostics": diagnostics,
    }
    _write_json_atomic(workdir / "vcneb_summary.json", summary)
    with (workdir / "vcneb_summary.txt").open("w", encoding="utf-8") as handle:
        handle.write(f"Forward barrier (enthalpy) = {barrier:.8f} eV\n")
        handle.write(f"Reaction enthalpy          = {delta:.8f} eV\n")
        handle.write(f"Final max generalized force = {max_force:.8f} eV/A\n")
        handle.write(
            "Highest image diagnostics    = "
            + json.dumps(saddle, ensure_ascii=False, sort_keys=True)
            + "\n"
        )
        handle.write("Per-image physical diagnostics = vcneb_summary.json[path_diagnostics]\n")
        handle.write("Image enthalpies (eV)      = " + " ".join(f"{value:.8f}" for value in chain.enthalpies) + "\n")
    print(
        f"[DONE] barrier={barrier:.6f} eV delta={delta:.6f} eV "
        f"max_force={max_force:.6f} eV/A workdir={workdir}"
    )


if __name__ == "__main__":
    main()

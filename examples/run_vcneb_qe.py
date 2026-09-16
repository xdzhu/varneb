"""Run ordinary or climbing VCNEB with Quantum ESPRESSO through ASE.

All images are static ``pw.x`` evaluations.  QE must never perform ``relax``
or ``vc-relax`` itself: atomic and cell updates belong exclusively to VARNEB.
Use ``--validate-only`` to construct the 7-total-image BTO layout and verify
calculator capabilities without launching DFT.
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
    interpolate_vcneb,
    endpoint_structure_record,
    path_geometry_diagnostics,
    run_vcneb,
    validate_image_calculators,
)
from vcneb.executor import ThreadedCalculatorExecutor
from vcneb.qe import (
    attach_qe_calculators,
    make_ase_espresso_factory,
    load_approved_qe_pseudopotential_manifest,
    static_qe_input_data,
    validate_qe_pseudopotentials,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--initial", required=True, help="ASE-readable initial endpoint")
    parser.add_argument("--final", required=True, help="ASE-readable final endpoint")
    parser.add_argument("--workdir", default="qe_vcneb_run")
    parser.add_argument(
        "--n-images", type=int, default=7,
        help="Total images including fixed endpoints (7 means five interior workers)",
    )
    parser.add_argument("--fmax", type=float, default=0.10, help="NEB convergence threshold in eV/A")
    parser.add_argument("--steps", type=int, default=300)
    parser.add_argument("--k", type=float, default=0.10)
    parser.add_argument("--pressure-gpa", type=float, default=0.0)
    parser.add_argument("--command", default=os.environ.get("QE_COMMAND"), help="QE launcher, e.g. 'srun pw.x'")
    parser.add_argument("--pseudo-dir", default=os.environ.get("ESPRESSO_PSEUDO"))
    parser.add_argument("--pp", action="append", default=[], metavar="SPECIES=FILE")
    parser.add_argument("--pp-manifest", default=None, help="Approved JSON manifest that pins QE UPF filenames and MD5")
    parser.add_argument("--ecutwfc", type=float, default=100.0, help="QE wavefunction cutoff in Ry")
    parser.add_argument("--ecutrho", type=float, default=None, help="QE charge-density cutoff in Ry")
    parser.add_argument("--scf-thr", type=float, default=1e-8, help="QE electron convergence threshold in Ry")
    parser.add_argument("--kpts", type=int, nargs=3, default=[4, 4, 4], metavar=("NX", "NY", "NZ"))
    parser.add_argument("--smearing", default=None, help="Optional QE smearing name")
    parser.add_argument("--degauss", type=float, default=None, help="QE degauss in Ry; requires --smearing")
    parser.add_argument("--optimizer", choices=["FIRE", "BFGS", "LBFGS", "BFGSLineSearch"], default="FIRE")
    parser.add_argument("--maxstep", type=float, default=None)
    parser.add_argument("--no-climb", action="store_true", help="Use ordinary NEB without CI")
    parser.add_argument("--climb-after", type=int, default=None)
    parser.add_argument("--mic", action="store_true")
    parser.add_argument("--cell-interpolation", choices=["linear", "log_strain"], default="linear")
    parser.add_argument("--mapping", choices=["identity", "auto"], default="identity")
    parser.add_argument("--align-translation", action="store_true")
    parser.add_argument("--minimum-distance", type=float, default=None)
    parser.add_argument("--maximum-deformation", type=float, default=None)
    parser.add_argument("--image-workers", type=int, default=0)
    parser.add_argument("--image-retries", type=int, default=0)
    parser.add_argument("--image-manifest", default=None)
    parser.add_argument("--image-cache-dir", default=None)
    parser.add_argument("--image-cache-namespace", default=None)
    parser.add_argument("--validate-only", action="store_true")
    return parser.parse_args()


def parse_species_files(values: list[str]) -> dict[str, str]:
    result: dict[str, str] = {}
    for value in values:
        if "=" not in value:
            raise ValueError(f"--pp expects SPECIES=FILE, got {value!r}")
        species, filename = value.split("=", 1)
        if not species or not filename or species in result:
            raise ValueError(f"invalid or duplicate --pp value {value!r}")
        result[species] = filename
    return result


def build_qe_parameters(
    args: argparse.Namespace,
    symbols: set[str],
    *,
    pseudopotentials: dict[str, str] | None = None,
) -> dict:
    pseudopotentials = parse_species_files(args.pp) if pseudopotentials is None else dict(pseudopotentials)
    missing = sorted(symbols - set(pseudopotentials))
    if missing:
        raise ValueError("missing --pp mappings for " + ", ".join(missing))
    if args.ecutrho is not None and args.ecutrho < args.ecutwfc:
        raise ValueError("--ecutrho must be at least --ecutwfc")
    if (args.smearing is None) != (args.degauss is None):
        raise ValueError("--smearing and --degauss must be supplied together")
    system = {"ecutwfc": args.ecutwfc}
    if args.ecutrho is not None:
        system["ecutrho"] = args.ecutrho
    if args.smearing is not None:
        system.update({"occupations": "smearing", "smearing": args.smearing, "degauss": args.degauss})
    parameters = {
        "pseudopotentials": pseudopotentials,
        "kpts": tuple(args.kpts),
        "input_data": {
            "system": system,
            "electrons": {"conv_thr": args.scf_thr},
        },
    }
    # Store exactly the normalized static input sent to the factory in the
    # preflight/summary provenance, rather than only applying it internally.
    parameters["input_data"] = static_qe_input_data(parameters["input_data"])
    return parameters


def _git_revision() -> str | None:
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def _write_json_atomic(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=path.parent, delete=False) as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
        temporary = Path(handle.name)
    temporary.replace(path)


def main() -> None:
    args = parse_args()
    if args.image_workers < 0 or args.image_retries < 0:
        raise ValueError("--image-workers and --image-retries must be non-negative")
    if args.command is None or args.pseudo_dir is None:
        raise ValueError("QE VCNEB requires --command and --pseudo-dir (or QE_COMMAND/ESPRESSO_PSEUDO)")
    initial, final = read(args.initial), read(args.final)
    if set(initial.get_chemical_symbols()) != set(final.get_chemical_symbols()):
        raise ValueError("endpoint species sets differ")
    workdir = Path(args.workdir).resolve()
    workdir.mkdir(parents=True, exist_ok=True)
    images = interpolate_vcneb(
        initial,
        final,
        args.n_images,
        align_cells=True,
        mic=args.mic,
        cell_interpolation=args.cell_interpolation,
        mapping=None if args.mapping == "identity" else "auto",
        align_translation=args.align_translation,
        minimum_distance=args.minimum_distance,
        maximum_deformation=args.maximum_deformation,
    )
    write(workdir / "initial-vcneb.traj", images)
    symbols = set(initial.get_chemical_symbols())
    if args.pp_manifest and args.pp:
        raise ValueError("pass either --pp-manifest or one or more --pp mappings, not both")
    pp_manifest = (
        load_approved_qe_pseudopotential_manifest(args.pp_manifest, symbols)
        if args.pp_manifest
        else None
    )
    parameters = build_qe_parameters(
        args,
        symbols,
        pseudopotentials=None if pp_manifest is None else pp_manifest["pseudopotentials"],
    )
    pseudopotential_reports = validate_qe_pseudopotentials(
        args.pseudo_dir,
        parameters["pseudopotentials"],
        expected_md5=None if pp_manifest is None else pp_manifest["expected_md5"],
    )
    factory = make_ase_espresso_factory(
        parameters=parameters,
        command=args.command,
        pseudo_dir=args.pseudo_dir,
    )
    attach_qe_calculators(images, workdir=workdir, factory=factory)
    reports = validate_image_calculators(
        images,
        require_stress=True,
        require_variable_cell=True,
        require_directory=True,
        require_unique_directories=True,
    )
    image_manifest = Path(args.image_manifest).resolve() if args.image_manifest else workdir / "image_worker_manifest.jsonl"
    metadata = {
        "git_revision": _git_revision(),
        "command_line": sys.argv,
        "calculator": "ASE Quantum ESPRESSO",
        "n_images": args.n_images,
        "n_interior_images": max(0, args.n_images - 2),
        "endpoint_evaluation_policy": "fixed_cached_once" if args.image_workers else "ASE_calculator_cache",
        "image_workers": args.image_workers,
        "image_retries": args.image_retries,
        "cell_interpolation": args.cell_interpolation,
        "mapping": args.mapping,
        "align_translation": args.align_translation,
        "fmax_target_eV_per_A": args.fmax,
        "calculator_parameters": parameters,
        "pseudopotential_reports": pseudopotential_reports,
        "pseudopotential_manifest": pp_manifest,
        "calculator_reports": [report.to_dict() for report in reports],
        "initial_path_geometry": path_geometry_diagnostics(images),
        "endpoint_structures": {
            "initial": endpoint_structure_record(initial),
            "final": endpoint_structure_record(final),
        },
    }
    _write_json_atomic(workdir / "vcneb_preflight.json", {"status": "ok", **metadata})
    if args.validate_only:
        print(f"[OK] QE VCNEB preflight passed; report={workdir / 'vcneb_preflight.json'}")
        return
    executor = ThreadedCalculatorExecutor(
        args.image_workers,
        max_retries=args.image_retries,
        manifest_path=image_manifest,
        cache_dir=args.image_cache_dir,
        cache_namespace=args.image_cache_namespace,
    ) if args.image_workers else None
    chain, _ = run_vcneb(
        images,
        pressure_gpa=args.pressure_gpa,
        k=args.k,
        climb=not args.no_climb,
        climb_after=args.climb_after,
        image_executor=executor,
        optimizer=args.optimizer,
        optimizer_kwargs={} if args.maxstep is None else {"maxstep": args.maxstep},
        fmax=args.fmax,
        steps=args.steps,
        logfile=workdir / "vcneb.opt.log",
        trajectory=workdir / "vcneb.traj",
        snapshot_dir=workdir / "snapshots",
        failure_report=workdir / "vcneb_failure.json",
    )
    barrier, reaction = chain.barrier()
    summary = {
        "status": "completed",
        **metadata,
        "barrier_enthalpy_eV": barrier,
        "reaction_enthalpy_eV": reaction,
        "final_max_generalized_force_eV_per_A": chain.gradient_norm(-chain.get_forces()),
        "image_enthalpies_eV": [float(value) for value in chain.enthalpies],
        "path_diagnostics": chain.path_diagnostics(),
        "saddle_diagnostics": chain.saddle_diagnostics(),
        "image_manifest": str(image_manifest) if executor is not None else None,
    }
    _write_json_atomic(workdir / "vcneb_summary.json", summary)
    for index, image in enumerate(chain.images):
        write(workdir / f"{index:02d}" / "POSCAR.final", image, format="vasp", direct=True, vasp5=True)
    print(f"[DONE] barrier={barrier:.6f} eV reaction={reaction:.6f} eV workdir={workdir}")


if __name__ == "__main__":
    main()

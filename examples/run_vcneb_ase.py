"""Run a material VCNEB path through any ASE calculator or VARNEB factory.

This is the calculator-agnostic production driver.  Use ``--calculator`` for
an ASE calculator class with JSON-serializable constructor parameters, or
``--factory`` for a VARNEB/specialized ASE factory (for example CP2K, QE,
ABINIT, or LAMMPS) and pass backend-only keyword arguments through
``--factory-kwargs``.  Endpoints are fixed and evaluated once by the VCNEB
controller; only interior images are sent to the threaded executor.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib
import json
import math
from pathlib import Path
import sys

import numpy as np
from ase.io import read, write

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from vcneb import (  # noqa: E402
    VCNEB,
    attach_image_calculators,
    endpoint_structure_record,
    interpolate_vcneb,
    load_vcneb_subspace_artifact,
    make_ase_calculator_factory,
    run_vcneb,
    validate_static_endpoint_identity,
    validate_path_geometry,
    validate_image_calculators,
)
from vcneb.executor import ThreadedCalculatorExecutor  # noqa: E402
from vcneb.mode_subspace import _vcneb_x  # noqa: E402


def _load_symbol(spec: str):
    if ":" not in spec:
        raise ValueError("symbol must use module:attribute form")
    module_name, attribute = spec.split(":", 1)
    if not module_name or not attribute:
        raise ValueError(f"invalid symbol specification {spec!r}")
    return getattr(importlib.import_module(module_name), attribute)


def _json_object(path: str | None) -> dict:
    if path is None:
        return {}
    value = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"JSON object expected in {path}")
    return value


def _validate_static_endpoint_identity(summary_path: str | Path, images) -> dict:
    """Require the static gate to describe the effective VCNEB endpoints.

    Atom mapping and translation alignment can change the structures that the
    calculator actually sees.  A static calculation of the source files is
    therefore not a valid gate unless its endpoint hashes match the mapped,
    aligned band endpoints exactly.
    """

    summary_path = Path(summary_path).resolve()
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    report = validate_static_endpoint_identity(summary, images)
    report["summary"] = str(summary_path)
    return report


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--initial", required=True)
    parser.add_argument("--final", required=True)
    seed_group = parser.add_mutually_exclusive_group()
    seed_group.add_argument(
        "--resume-snapshot",
        default=None,
        help="complete chain snapshot to resume from in a new work directory",
    )
    seed_group.add_argument(
        "--initial-chain",
        help="unconstrained starting chain to project into a declared mode subspace",
    )
    parser.add_argument(
        "--subspace-artifact",
        help="audited global VCNEB mode-subspace JSON, independent of calculator backend",
    )
    parser.add_argument(
        "--subspace-artifact-sha256",
        help="pin the exact subspace artifact bytes; required for a production subspace run",
    )
    parser.add_argument("--workdir", required=True)
    parser.add_argument("--n-images", type=int, default=7)
    parser.add_argument("--fmax", type=float, default=0.10)
    parser.add_argument("--steps", type=int, default=300)
    parser.add_argument("--k", type=float, default=0.20)
    parser.add_argument("--pressure-gpa", type=float, default=0.0)
    parser.add_argument("--optimizer", default="FIRE")
    parser.add_argument("--maxstep", type=float, default=None)
    parser.add_argument("--image-workers", type=int, default=0)
    parser.add_argument("--image-retries", type=int, default=0)
    parser.add_argument(
        "--candidate-step-retries",
        type=int,
        default=0,
        help="FIRE geometry backtracking retries before an electronic evaluation",
    )
    parser.add_argument("--command", default=None)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--calculator", help="ASE class/factory as module:attribute")
    group.add_argument("--factory", help="VARNEB factory as module:attribute")
    parser.add_argument("--parameters", default=None, help="JSON calculator parameter object")
    parser.add_argument("--factory-kwargs", default=None, help="JSON factory-only keyword object")
    parser.add_argument(
        "--endpoint-static-summary",
        default=None,
        help="static endpoint summary whose structure hashes must match the effective band endpoints",
    )
    parser.add_argument("--no-climb", action="store_true")
    parser.add_argument("--climb-after", type=int, default=None)
    parser.add_argument("--mic", action="store_true")
    parser.add_argument("--cell-interpolation", choices=("linear", "log_strain"), default="log_strain")
    parser.add_argument("--mapping", choices=("identity", "auto"), default="auto")
    parser.add_argument(
        "--no-align-cells",
        action="store_true",
        help="preserve the endpoint cell frames instead of removing a global rotation",
    )
    parser.add_argument("--align-translation", action="store_true")
    parser.add_argument("--minimum-distance", type=float, default=None)
    parser.add_argument("--maximum-deformation", type=float, default=None)
    parser.add_argument(
        "--minimum-endpoint-separation",
        type=float,
        default=None,
        help="reject a nearly identical endpoint pair before DFT (extended-coordinate Angstrom)",
    )
    parser.add_argument(
        "--maximum-cell-step",
        type=float,
        default=None,
        help="maximum relative Frobenius cell deformation per optimizer step",
    )
    parser.add_argument("--validate-only", action="store_true")
    parser.add_argument(
        "--static-only",
        action="store_true",
        help="evaluate the fixed endpoints once and write ase_static_summary.json",
    )
    parser.add_argument(
        "--static-endpoints",
        choices=("both", "initial", "final"),
        default="both",
        help="endpoints evaluated by --static-only (default: both)",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.n_images < 3 or args.image_workers < 0 or args.image_retries < 0:
        raise ValueError("n_images must be >=3 and worker/retry counts non-negative")
    if not args.validate_only and args.calculator is None and args.factory is None:
        raise ValueError("a calculator or factory is required unless --validate-only is set")
    if args.subspace_artifact_sha256 is not None and args.subspace_artifact is None:
        raise ValueError("--subspace-artifact-sha256 requires --subspace-artifact")
    if args.subspace_artifact is not None and not args.validate_only and args.subspace_artifact_sha256 is None:
        raise ValueError("production subspace runs require --subspace-artifact-sha256")
    if args.no_align_cells and args.cell_interpolation != "linear":
        raise ValueError("--no-align-cells requires --cell-interpolation linear")
    initial = read(args.initial)
    final = read(args.final)
    workdir = Path(args.workdir).resolve()
    workdir.mkdir(parents=True, exist_ok=True)
    supplied_chain = args.resume_snapshot or args.initial_chain
    if supplied_chain is None:
        images = interpolate_vcneb(
            initial,
            final,
            args.n_images,
            align_cells=not args.no_align_cells,
            mic=args.mic,
            cell_interpolation=args.cell_interpolation,
            mapping=None if args.mapping == "identity" else "auto",
            align_translation=args.align_translation,
            minimum_distance=args.minimum_distance,
            maximum_deformation=args.maximum_deformation,
        )
    else:
        images = read(supplied_chain, index=":")
        if len(images) != args.n_images:
            raise ValueError(
                f"supplied chain contains {len(images)} images; expected {args.n_images}"
            )
        for label, expected, actual in (
            ("initial", endpoint_structure_record(initial), endpoint_structure_record(images[0])),
            ("final", endpoint_structure_record(final), endpoint_structure_record(images[-1])),
        ):
            if expected["sha256"] != actual["sha256"]:
                raise ValueError(f"supplied chain {label} endpoint does not match the requested endpoint")
    geometry = validate_path_geometry(
        images,
        minimum_distance=args.minimum_distance,
        maximum_deformation=args.maximum_deformation,
        minimum_endpoint_separation=args.minimum_endpoint_separation,
    )
    endpoint_identity_gate = (
        _validate_static_endpoint_identity(args.endpoint_static_summary, images)
        if args.endpoint_static_summary is not None
        else None
    )
    mode_basis = None
    mode_scale = None
    subspace_record = None
    if args.subspace_artifact is not None:
        artifact_path = Path(args.subspace_artifact).resolve()
        artifact_sha256 = hashlib.sha256(artifact_path.read_bytes()).hexdigest()
        if args.subspace_artifact_sha256 is not None and (
            len(args.subspace_artifact_sha256) != 64
            or args.subspace_artifact_sha256.lower() != artifact_sha256
        ):
            raise ValueError("subspace artifact SHA256 does not match the pinned digest")
        mode_basis, artifact = load_vcneb_subspace_artifact(
            artifact_path, images[0], images[-1],
        )
        mode_scale = float(artifact["cell_scale_A"])
        unprojected = [image.copy() for image in images]
        projected = [image.copy() for image in images]
        VCNEB(
            projected, cell_scale=mode_scale, mode_basis=mode_basis,
            constraint_mode="subspace", k=args.k, climb=False,
        )
        first_cell = images[0].cell.array
        before = [_vcneb_x(image, first_cell, mode_scale) for image in unprojected]
        displacements = [float(np.linalg.norm(
            _vcneb_x(image, first_cell, mode_scale) - original
        )) for image, original in zip(projected, before)]
        if args.resume_snapshot is not None and max(displacements) > 1e-8:
            raise ValueError("resume snapshot is outside the declared global mode subspace")
        if args.resume_snapshot is None:
            raw_path = workdir / "initial-vcneb-unprojected.traj"
            if raw_path.exists():
                raise FileExistsError(f"refusing to replace raw initial chain: {raw_path}")
            write(raw_path, unprojected)
        images = projected
        geometry = validate_path_geometry(
            images,
            minimum_distance=args.minimum_distance,
            maximum_deformation=args.maximum_deformation,
            minimum_endpoint_separation=args.minimum_endpoint_separation,
        )
        subspace_record = {
            "artifact": str(artifact_path),
            "artifact_sha256": artifact_sha256,
            "artifact_sha256_pinned": args.subspace_artifact_sha256 is not None,
            "subspace_kind": artifact["subspace_kind"],
            "reference_id": artifact["reference_id"],
            "n_directions": artifact["n_directions"],
            "cell_scale_A": mode_scale,
            "initial_projection_displacements_vcneb_A": displacements,
            "source_sha256": artifact["source_sha256"],
        }
    write(workdir / "initial-vcneb.traj", images)

    parameters = _json_object(args.parameters)
    metadata = {
        "status": "preflight",
        "driver": "examples/run_vcneb_ase.py",
        "calculator": args.calculator,
        "factory": args.factory,
        "calculator_parameters": parameters,
        "factory_kwargs": _json_object(args.factory_kwargs),
        "resume_snapshot": args.resume_snapshot,
        "initial_chain": args.initial_chain,
        "n_images": args.n_images,
        "n_interior_images": args.n_images - 2,
        "external_pressure_gpa": args.pressure_gpa,
        "fmax_target_eV_per_A": args.fmax,
        "candidate_step_retries": args.candidate_step_retries,
        "maximum_cell_step": args.maximum_cell_step,
        "endpoint_evaluation_policy": "fixed_cached_once",
        "cell_interpolation": args.cell_interpolation,
        "mapping": args.mapping,
        "align_cells": not args.no_align_cells,
        "align_translation": args.align_translation,
        "endpoint_structures": {
            "initial": endpoint_structure_record(images[0]),
            "final": endpoint_structure_record(images[-1]),
        },
        "endpoint_static_identity_gate": endpoint_identity_gate,
        "initial_path_geometry": geometry,
        "mode_subspace": subspace_record,
        "calculator_validation": "not_instantiated_validate_only" if args.validate_only else "factory_configuration_only",
        "calculator_reports": [],
    }
    (workdir / "vcneb_preflight.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    if args.validate_only:
        print(f"[OK] ASE VCNEB preflight passed: {workdir / 'vcneb_preflight.json'}")
        return

    if args.factory:
        builder = _load_symbol(args.factory)
        factory_kwargs = _json_object(args.factory_kwargs)
        factory_kwargs["parameters"] = parameters
        if args.command is not None:
            factory_kwargs["command"] = args.command
        factory = builder(**factory_kwargs)
    else:
        calculator = _load_symbol(args.calculator)
        factory = make_ase_calculator_factory(
            calculator,
            parameters=parameters,
            command=args.command,
        )

    # Some ASE calculators (notably CP2K) start a persistent external process
    # in their constructor.  Instantiate image calculators only inside an
    # actual scheduler run, never during a geometry-only validation on a login
    # node.  Production still validates stress support and isolated image
    # directories before the first electronic-structure evaluation.
    attach_image_calculators(images, workdir=workdir, factory=factory)
    reports = validate_image_calculators(
        images,
        require_stress=True,
        require_variable_cell=True,
        require_directory=True,
        require_unique_directories=True,
    )
    metadata["calculator_validation"] = "runtime_instantiated"
    metadata["calculator_reports"] = [report.to_dict() for report in reports]
    (workdir / "vcneb_preflight.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    if args.static_only:
        endpoint_indices = {
            "both": (0, len(images) - 1),
            "initial": (0,),
            "final": (len(images) - 1,),
        }[args.static_endpoints]
        endpoint_results = []
        for index in endpoint_indices:
            image = images[index]
            energy = float(image.get_potential_energy())
            forces = image.get_forces()
            stress = image.get_stress(voigt=False)
            if not all(math.isfinite(float(value)) for value in forces.ravel()):
                raise RuntimeError(f"non-finite endpoint forces at image {index}")
            if not all(math.isfinite(float(value)) for value in stress.ravel()):
                raise RuntimeError(f"non-finite endpoint stress at image {index}")
            endpoint_results.append(
                {
                    "index": index,
                    "label": "initial" if index == 0 else "final",
                    "energy_eV": energy,
                    "max_force_eV_per_A": max(
                        math.sqrt(sum(float(value) ** 2 for value in row))
                        for row in forces
                    ),
                    "stress_eV_per_A3": stress.tolist(),
                    "n_atoms": len(image),
                }
            )
        static_summary = {
            **metadata,
            "status": "static_completed",
            "static_endpoints": endpoint_results,
            "endpoint_evaluation_policy": "fixed_cached_once",
        }
        (workdir / "ase_static_summary.json").write_text(
            json.dumps(static_summary, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        print(json.dumps(static_summary, indent=2, sort_keys=True))
        return

    executor = (
        ThreadedCalculatorExecutor(
            args.image_workers,
            max_retries=args.image_retries,
            manifest_path=workdir / "image_worker_manifest.jsonl",
            cache_dir=workdir / "image_cache",
            cache_namespace=json.dumps(
                {"calculator": args.calculator, "factory": args.factory, "parameters": parameters},
                sort_keys=True,
            ),
        )
        if args.image_workers
        else None
    )
    chain, _ = run_vcneb(
        images,
        pressure_gpa=args.pressure_gpa,
        cell_scale=mode_scale,
        k=args.k,
        climb=not args.no_climb,
        climb_after=args.climb_after,
        mode_basis=mode_basis,
        constraint_mode="subspace" if mode_basis is not None else None,
        image_executor=executor,
        optimizer=args.optimizer,
        optimizer_kwargs={} if args.maxstep is None else {"maxstep": args.maxstep},
        candidate_step_retries=args.candidate_step_retries,
        maximum_cell_step=args.maximum_cell_step,
        fmax=args.fmax,
        steps=args.steps,
        logfile=workdir / "vcneb.opt.log",
        trajectory=workdir / "vcneb.traj",
        snapshot_dir=workdir / "snapshots",
        candidate_step_manifest=workdir / "candidate_step_manifest.jsonl",
        failure_report=workdir / "vcneb_failure.json",
    )
    barrier, reaction = chain.barrier()
    final_force = chain.gradient_norm(-chain.get_forces())
    converged = bool(final_force <= args.fmax)
    summary = {
        **metadata,
        "status": "converged" if converged else "max_steps_reached",
        "converged": converged,
        "termination": "force_threshold" if converged else "step_limit",
        "barrier_enthalpy_eV": barrier,
        "reaction_enthalpy_eV": reaction,
        "final_max_generalized_force_eV_per_A": final_force,
        "image_enthalpies_eV": [float(value) for value in chain.enthalpies],
        "path_diagnostics": chain.path_diagnostics(),
        "saddle_diagnostics": chain.saddle_diagnostics(),
    }
    (workdir / "vcneb_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

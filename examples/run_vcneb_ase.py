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
import importlib
import json
import math
from pathlib import Path
import sys

from ase.io import read, write

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from vcneb import (  # noqa: E402
    attach_image_calculators,
    endpoint_structure_record,
    interpolate_vcneb,
    make_ase_calculator_factory,
    path_geometry_diagnostics,
    run_vcneb,
    validate_image_calculators,
)
from vcneb.executor import ThreadedCalculatorExecutor  # noqa: E402


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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--initial", required=True)
    parser.add_argument("--final", required=True)
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
    parser.add_argument("--command", default=None)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--calculator", help="ASE class/factory as module:attribute")
    group.add_argument("--factory", help="VARNEB factory as module:attribute")
    parser.add_argument("--parameters", default=None, help="JSON calculator parameter object")
    parser.add_argument("--factory-kwargs", default=None, help="JSON factory-only keyword object")
    parser.add_argument("--no-climb", action="store_true")
    parser.add_argument("--climb-after", type=int, default=None)
    parser.add_argument("--mic", action="store_true")
    parser.add_argument("--cell-interpolation", choices=("linear", "log_strain"), default="log_strain")
    parser.add_argument("--mapping", choices=("identity", "auto"), default="auto")
    parser.add_argument("--align-translation", action="store_true")
    parser.add_argument("--minimum-distance", type=float, default=None)
    parser.add_argument("--maximum-deformation", type=float, default=None)
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
    initial = read(args.initial)
    final = read(args.final)
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

    parameters = _json_object(args.parameters)
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
    attach_image_calculators(images, workdir=workdir, factory=factory)
    reports = validate_image_calculators(
        images,
        require_stress=True,
        require_variable_cell=True,
        require_directory=True,
        require_unique_directories=True,
    )
    metadata = {
        "status": "preflight",
        "driver": "examples/run_vcneb_ase.py",
        "calculator": args.calculator,
        "factory": args.factory,
        "calculator_parameters": parameters,
        "factory_kwargs": _json_object(args.factory_kwargs),
        "n_images": args.n_images,
        "n_interior_images": args.n_images - 2,
        "fmax_target_eV_per_A": args.fmax,
        "endpoint_evaluation_policy": "fixed_cached_once",
        "cell_interpolation": args.cell_interpolation,
        "mapping": args.mapping,
        "align_translation": args.align_translation,
        "endpoint_structures": {
            "initial": endpoint_structure_record(images[0]),
            "final": endpoint_structure_record(images[-1]),
        },
        "initial_path_geometry": path_geometry_diagnostics(
            images,
            minimum_distance=args.minimum_distance,
            maximum_deformation=args.maximum_deformation,
        ),
        "calculator_reports": [report.to_dict() for report in reports],
    }
    (workdir / "vcneb_preflight.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    if args.validate_only:
        print(f"[OK] ASE VCNEB preflight passed: {workdir / 'vcneb_preflight.json'}")
        return

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
                    "max_force_eV_per_A": float(abs(forces).max()),
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
        **metadata,
        "status": "completed",
        "barrier_enthalpy_eV": barrier,
        "reaction_enthalpy_eV": reaction,
        "final_max_generalized_force_eV_per_A": chain.gradient_norm(-chain.get_forces()),
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

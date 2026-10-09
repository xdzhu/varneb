"""Exercise the actual clamped path CLI on HF geometry, with zero calculators."""

import argparse
from contextlib import redirect_stdout
import hashlib
from importlib.metadata import version
import io
import json
from pathlib import Path
import sys


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    source, output = args.source.resolve(), args.output.resolve()
    if not source.is_dir() or output.exists() or output.is_relative_to(source):
        raise ValueError("immutable source and a fresh separate output namespace required")
    sys.path.insert(0, str(source))
    from vcneb import material_runner
    from ase.io import read, write
    if not Path(material_runner.__file__).resolve().is_relative_to(source):
        raise ValueError("imports must come from the declared tested source")

    def forbidden(_):
        raise AssertionError("geometry preflight must not load a calculator")

    seeds = source / "benchmarks/hfo2_channels/20261008/clamped_endpoint_seeds"
    master_path = seeds / "clamped_seed_manifest.json"
    master_hash = digest(master_path)
    master = json.loads(master_path.read_text(encoding="utf-8"))
    if (len(master["seeds"]) != 10 or master["strain_conditions"] != [0., .01]
            or not master["allow_tilt"] or master["holdout_strain_reserved_not_generated"] != .005):
        raise ValueError("finite registered seed matrix changed")
    rows = []
    for index, entry in enumerate(master["seeds"]):
        manifest_path = seeds / entry["manifest_file"]
        if digest(manifest_path) != entry["manifest_sha256"]:
            raise ValueError("registered seed manifest changed")
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        final = manifest_path.parent / manifest["seed_file"]
        reference = manifest_path.parent.parent / "T/POSCAR.seed"
        if digest(final) != entry["seed_sha256"]:
            raise ValueError("registered geometry changed")
        arguments = ["--initial", str(reference), "--final", str(final),
                     "--clamped-plane-reference", str(reference), "--clamped-allow-tilt", "true",
                     "--cell-scale", str(master["cell_scale_A"]),
                     "--workdir", str(output / f"seed_{index:02d}"), "--n-images", "3",
                     "--no-align-cells", "--cell-interpolation", "linear", "--mapping", "identity",
                     "--validate-only"]
        with redirect_stdout(io.StringIO()):
            material_runner.main(arguments, symbol_loader=forbidden)
        report = json.loads((output / f"seed_{index:02d}/vcneb_preflight.json").read_text())
        if (report["calculator_reports"] or not report["requires_stress"]
                or report["mechanical_boundary"]["cell_dofs"] != 3
                or report["calculator_validation"] != "not_instantiated_validate_only"):
            raise ValueError("clamped preflight contract not honored")
        # A changed internal plane must be rejected, not projected or sent to DFT.
        chain = read(output / f"seed_{index:02d}/initial-vcneb.traj", index=":")
        chain[1].cell[0, 0] += .001
        raw_path = output / f"bad_raw_{index:02d}.traj"
        write(raw_path, chain)
        bad_arguments = arguments.copy()
        bad_arguments[bad_arguments.index("--workdir") + 1] = str(output / f"reject_{index:02d}")
        bad_arguments += ["--resume-snapshot", str(raw_path)]
        try:
            with redirect_stdout(io.StringIO()):
                material_runner.main(bad_arguments, symbol_loader=forbidden)
        except ValueError as error:
            if "substrate vectors" not in str(error):
                raise
            rejection = str(error)
        else:
            raise AssertionError("incompatible raw image was accepted")
        if (output / f"reject_{index:02d}/initial-vcneb.traj").exists():
            raise AssertionError("rejected geometry was persisted as a valid seed")
        if digest(manifest_path) != entry["manifest_sha256"] or digest(final) != entry["seed_sha256"]:
            raise ValueError("preflight modified scientific inputs")
        rows.append({"phase_label": entry["phase_label"], "strain": entry["strain"],
                     "seed_sha256": entry["seed_sha256"], "open_cell_dofs": 3,
                     "valid_preflight": True, "bad_internal_plane_rejected": rejection})
    if digest(master_path) != master_hash:
        raise ValueError("source master manifest changed")
    print(json.dumps({
        "status": "clamped_path_entry_geometry_only_passed_not_G2_material_results",
        "source": str(source), "output": str(output), "python": sys.executable,
        "packages": {name: version(name) for name in ("numpy", "ase")},
        "master_manifest_sha256": master_hash,
        "source_hashes": {str(path.relative_to(source)): digest(path) for path in (
            Path(material_runner.__file__), source / "vcneb/config.py", source / "vcneb/cli.py",
            source / "vcneb/epitaxial_boundary.py", Path(__file__).resolve(),
        )},
        "seeds": rows, "calculator_symbols_loaded": 0, "DFT_calls": 0,
        "job_submissions": 0, "G1_or_G2_gate_certified": False,
    }, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()

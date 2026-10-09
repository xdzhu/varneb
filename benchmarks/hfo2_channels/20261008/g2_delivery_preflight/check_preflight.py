"""Read-only G2 delivery check: no calculator evaluation or job submission."""

import argparse
from importlib.metadata import version
import json
from pathlib import Path
import subprocess
import sys

import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--probe-root", type=Path, required=True,
                        help="reserved nonexistent path; must remain uncreated")
    args = parser.parse_args()
    source, probe = args.source.resolve(), args.probe_root.resolve()
    if not source.is_dir() or probe.exists():
        raise ValueError("existing immutable source and nonexistent probe path required")
    sys.path.insert(0, str(source))

    import examples.hfo2_fixed_input_factory as transport
    import scripts.relax_clamped_ase_endpoint as entry
    from scripts.audit_hfo2_static_replica import sha256

    for module in (transport, entry):
        if not Path(module.__file__).resolve().is_relative_to(source):
            raise ValueError("refuse imports outside the declared archive")
    case = source / "benchmarks/hfo2_channels/20261008"
    seeds = case / "clamped_endpoint_seeds"
    master_path = seeds / "clamped_seed_manifest.json"
    master_before = sha256(master_path)
    master = json.loads(master_path.read_text(encoding="utf-8"))
    if (master["n_geometry_seeds"] != 10 or len(master["seeds"]) != 10
            or master["strain_conditions"] != [0., .01]
            or master["holdout_strain_reserved_not_generated"] != .005
            or master["physical_contract_sha256"] != transport.CONTRACT):
        raise ValueError("the registered finite seed set or calculator contract changed")
    parameters_path = case / "M_endpoint_factory_parameters.json"
    parameters = json.loads(parameters_path.read_text(encoding="utf-8"))
    physical_source = Path(parameters["source_directory"])
    physical_hashes = {name: sha256(physical_source / name) for name in transport.CONTRACT}
    if physical_hashes != transport.CONTRACT:
        raise ValueError("physical source contract changed")
    command = ("mpirun -np 32 /public/home/iai806/apprepo/abacus/"
               "v3.10.0LTS-intelmpi2025/app/bin/abacus")
    factory = transport.make_factory(parameters=parameters, command=command)
    rows, planes = [], {}
    for index, record in enumerate(master["seeds"]):
        manifest_path = seeds / record["manifest_file"]
        if not manifest_path.resolve().is_relative_to(seeds):
            raise ValueError("seed manifest must stay in the registered namespace")
        if sha256(manifest_path) != record["manifest_sha256"]:
            raise ValueError("seed manifest changed")
        atoms, boundary, manifest = entry.load_seed(manifest_path)
        if (atoms.calc is not None or manifest["pressure_gpa"] != 0.
                or not boundary.allow_tilt
                or manifest["physical_contract_sha256"] != transport.CONTRACT):
            raise ValueError("fresh, tilt-released P=0 geometry seed required")
        plane = atoms.cell.array[:2].copy()
        planes.setdefault(record["strain"], plane)
        np.testing.assert_allclose(plane, planes[record["strain"]], atol=1e-14, rtol=0)
        calculator = factory(index, atoms, probe / f"endpoint_{index:02d}")
        if calculator.results or calculator.next_call != 0 or probe.exists():
            raise ValueError("constructor must not evaluate or create a work directory")
        if sha256(manifest_path) != record["manifest_sha256"]:
            raise ValueError("preflight changed a source manifest")
        rows.append({"phase_label": manifest["phase_label"], "strain": manifest["strain"],
                     "seed_sha256": manifest["seed_sha256"],
                     "manifest_sha256": record["manifest_sha256"],
                     "substrate_drift_A": float(np.abs(plane - boundary.reference_cell[:2]).max()),
                     "calculator_evaluated": False})
    np.testing.assert_allclose(planes[.01], 1.01 * planes[0.], atol=1e-12, rtol=0)
    job_script = source / "cluster/hf_hfo2_clamped_endpoint_20261008.slurm"
    syntax = subprocess.run(["bash", "-n", str(job_script)], capture_output=True,
                            text=True, check=True)
    if sha256(master_path) != master_before or probe.exists():
        raise ValueError("preflight unexpectedly changed its inputs or reserved path")
    report = {"status": "read_only_delivery_preflight_passed_not_G2_material_passed",
              "source": str(source), "python": sys.executable,
              "packages": {name: version(name) for name in ("numpy", "ase")},
              "validator_sha256": sha256(Path(__file__)),
              "master_manifest_sha256": master_before,
              "parameters_sha256": sha256(parameters_path),
              "physical_source": str(physical_source), "physical_source_sha256": physical_hashes,
              "source_hashes": {str(path.relative_to(source)): sha256(path) for path in
                                (Path(transport.__file__), Path(entry.__file__), job_script)},
              "shell_syntax_exit_code": syntax.returncode, "command_not_executed": command,
              "seeds": rows, "new_DFT_calls": 0, "job_submissions": 0,
              "probe_directory_created": False, "G1_or_G2_material_gate_certified": False}
    print(json.dumps(report, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()

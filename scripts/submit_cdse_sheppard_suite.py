"""Submit six independently tracked CdSe jobs, never touching GaN jobs."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile


def atomic_json(path, data):
    with tempfile.NamedTemporaryFile("w", dir=path.parent, delete=False) as handle:
        json.dump(data, handle, indent=2)
        handle.write("\n")
        temporary = Path(handle.name)
    temporary.replace(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", required=True, type=Path)
    parser.add_argument("--run", action="store_true")
    parser.add_argument("--phase-candidates", type=Path,
                        help="submit only the four recovery check/production jobs; never repeat completed relaxations")
    parser.add_argument("--lattice-probe", type=Path)
    parser.add_argument("--routes", nargs="+", choices=("cell_mapping", "atomic_mapping"),
                        default=["cell_mapping", "atomic_mapping"])
    args = parser.parse_args()
    repo = args.repo.resolve()
    manifest = repo / "launch.json"
    if manifest.exists():
        raise FileExistsError("launch manifest exists; inspect IDs, do not resubmit")
    for label in ("rs", "wz"):
        directory = args.phase_candidates / label if args.phase_candidates else repo / "cases" / "endpoint_seeds" / label
        names = ("CONTCAR", "INCAR", "KPOINTS", "POTCAR", "phase_candidate_audit.json") if args.phase_candidates else ("POSCAR", "INCAR", "KPOINTS", "POTCAR", "seed.json")
        for name in names:
            path = directory / name
            if not path.is_file():
                raise FileNotFoundError(path)
    if not args.run:
        print("Prepared seeds found; no submission without --run")
        return
    (repo / "logs").mkdir(exist_ok=True)
    payload = {"status": "submitting", "repo": str(repo), "jobs": {},
               "endpoint_recovery": "bounded phase candidates, full-SCF raw force/virial gate; original relaxations not repeated" if args.phase_candidates else None,
               "source_sha256": {str(p.relative_to(repo)): hashlib.sha256(p.read_bytes()).hexdigest()
                   for base in ("vcneb", "scripts", "cluster", "examples")
                   for p in sorted((repo / base).rglob("*")) if p.suffix in (".py", ".slurm")},
               "resources": {"endpoint_or_check": {"nodes": 1, "mpi": 32},
                             "each_production": {"nodes": 2, "allocated_tasks": 194, "workers": 6,
                                                 "mpi_per_worker": 32, "n_interiors": 15, "endpoint_workers": 0}},
               "physical_contract": {"xc": "PBE", "encut_eV": 455, "kmesh": [10,10,10],
                                     "kpoint_centering": "Monkhorst-Pack", "pressure_gpa": 0.,
                                     "fmax_eV_per_A": .10, "CI": False, "SYMPREC": 1e-4, "ISYM": -1}}
    # Exclusive claim precedes every sbatch, including partial failures.
    with manifest.open("x") as handle:
        json.dump(payload, handle, indent=2)

    def submit(role, script, values, dependencies=()):
        values = {"SUITE_REPO": str(repo), **values}
        if args.phase_candidates:
            values["ENDPOINT_CANDIDATES"] = str(args.phase_candidates.resolve())
        if args.lattice_probe:
            values["LATTICE_PROBE"] = str(args.lattice_probe.resolve(strict=True))
        command = ["sbatch", "--parsable", f"--job-name={role}",
                   f"--output={repo}/logs/%x-%j.out", f"--error={repo}/logs/%x-%j.err",
                   "--export=ALL,RUN_DFT=1," + ",".join(f"{k}={v}" for k,v in values.items())]
        if dependencies:
            command.append("--dependency=afterok:" + ":".join(dependencies))
        command.append(str(repo / "cluster" / script))
        result = subprocess.run(command, check=True, text=True, capture_output=True, cwd=repo)
        job_id = result.stdout.strip().split(";")[0]
        if not job_id.isdigit():
            raise ValueError(f"invalid job ID: {job_id!r}")
        payload["jobs"][role] = {"job_id": job_id, "afterok": list(dependencies), "command": command}
        atomic_json(manifest, payload)
        print(f"{role}: {job_id}; dependencies={dependencies}", flush=True)
        return job_id

    try:
        endpoints = () if args.phase_candidates else tuple(submit(f"cdse-{label}-bfgs", "hf_cdse_sheppard_endpoint.slurm", {"ENDPOINT": label}) for label in ("rs", "wz"))
        for route in dict.fromkeys(args.routes):
            check = submit(f"cdse-{route}-check", "hf_cdse_sheppard_validate.slurm", {"ROUTE": route}, endpoints)
            submit(f"cdse-{route}-vcneb", "hf_cdse_sheppard_vcneb.slurm", {"ROUTE": route}, (check,))
        payload["status"] = "submitted"
    except Exception as error:
        payload["status"] = "partial_submission_requires_inspection"
        payload["error"] = str(error)
        raise
    finally:
        atomic_json(manifest, payload)


if __name__ == "__main__":
    main()

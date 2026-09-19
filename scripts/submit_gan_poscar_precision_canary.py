"""Submit one immutable full-SCF serialization canary; never a production path."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", required=True, type=Path)
    parser.add_argument("--trajectory", required=True, type=Path)
    parser.add_argument("--source-input", required=True, type=Path)
    parser.add_argument("--run", action="store_true")
    args = parser.parse_args()
    repo = args.repo.resolve()
    manifest = repo / "canary_launch.json"
    if manifest.exists():
        raise FileExistsError("canary manifest exists; do not resubmit")
    for path in (args.trajectory, args.source_input / "INCAR",
                 args.source_input / "KPOINTS", args.source_input / "POTCAR"):
        if not path.is_file():
            raise FileNotFoundError(path)
    workdir = repo / "canary_image7_poscar17"
    payload = {
        "status": "prepared_not_submitted", "purpose": "one image full-SCF serialization canary",
        "production_continuation": False, "nodes": 1, "mpi_tasks": 32,
        "trajectory": str(args.trajectory.resolve()), "trajectory_sha256": hashlib.sha256(args.trajectory.read_bytes()).hexdigest(),
        "source_input": str(args.source_input.resolve()), "workdir": str(workdir),
        "physical_contract": {"vasp": "6.3.2", "xc": "PBE", "encut_eV": 600,
                              "kmesh": [8,8,6], "gamma": True, "isym": -1, "symprec": 1e-4},
        "source_sha256": {str(path.relative_to(repo)): hashlib.sha256(path.read_bytes()).hexdigest()
                          for base in ("vcneb", "scripts", "cluster")
                          for path in sorted((repo/base).rglob("*")) if path.suffix in (".py", ".slurm")},
    }
    if args.run:
        (repo / "logs").mkdir(exist_ok=True)
        command = ["sbatch", "--parsable", f"--output={repo}/logs/%x-%j.out",
                   f"--error={repo}/logs/%x-%j.err",
                   "--export=ALL,RUN_DFT=1," + ",".join((
                       f"SOURCE_REPO={repo}", f"OLD_TRAJECTORY={args.trajectory.resolve()}",
                       f"SOURCE_INPUT={args.source_input.resolve()}", f"CANARY_WORKDIR={workdir}")),
                   str(repo / "cluster/hf_gan_poscar_precision_canary.slurm")]
        result = subprocess.run(command, cwd=repo, text=True, capture_output=True, check=True)
        job_id = result.stdout.strip().split(";")[0]
        if not job_id.isdigit():
            raise ValueError(f"invalid Slurm job id: {job_id!r}")
        payload.update(status="submitted", job_id=job_id, command=command)
    with manifest.open("x", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)
        handle.write("\n")
    print(json.dumps({key: payload[key] for key in ("status", "production_continuation")
                      if key in payload} | ({"job_id": payload["job_id"]} if "job_id" in payload else {})))


if __name__ == "__main__":
    main()

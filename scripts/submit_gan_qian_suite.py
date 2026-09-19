"""Submit the prepared GaN suite with explicit afterok dependencies.

Requires --run; records each returned ID and refuses duplicate submission.
"""
import argparse
import hashlib
import json
import os
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
    args = parser.parse_args()
    repo = args.repo.resolve()
    manifest = repo / "launch.json"
    if manifest.exists():
        raise FileExistsError("launch manifest exists; inspect recorded jobs instead of resubmitting")
    for path in ("cases/endpoint_seeds/b3/POSCAR", "cases/endpoint_seeds/b1/POSCAR",
                 "cases/b4_b1_hexagonal/initial.traj"):
        if not (repo / path).is_file():
            raise FileNotFoundError(repo / path)
    if not args.run:
        print("Prepared suite found; no submission without --run")
        return
    environment = dict(os.environ, VCNEB_GIT_REVISION="d8af06db-dirty-gan-qian-suite")
    (repo / "logs").mkdir(exist_ok=True)
    payload = {"status": "submitting", "jobs": {}, "repo": str(repo),
               "source_sha256": {str(path.relative_to(repo)): hashlib.sha256(path.read_bytes()).hexdigest()
                   for path in sorted(repo.glob("cluster/hf_gan_qian_*.slurm"))}}
    atomic_json(manifest, payload)

    def submit(role, script, values, dependencies=()):
        command = ["sbatch", "--parsable", f"--job-name={role}",
                   f"--output={repo}/logs/%x-%j.out", f"--error={repo}/logs/%x-%j.err",
                   "--export=ALL,RUN_DFT=1," + ",".join(f"{key}={value}" for key, value in values.items())]
        if dependencies:
            command.append("--dependency=afterok:" + ":".join(dependencies))
        command.append(str(repo / "cluster" / script))
        job_id = subprocess.run(command, check=True, text=True, capture_output=True,
                                cwd=repo, env=environment).stdout.strip().split(";")[0]
        if not job_id.isdigit():
            raise ValueError(f"invalid Slurm job ID: {job_id!r}")
        payload["jobs"][role] = {"job_id": job_id, "afterok": list(dependencies),
                                 "command": command}
        atomic_json(manifest, payload)
        print(f"{role}: {job_id}; afterok={','.join(dependencies) or 'none'}", flush=True)
        return job_id

    try:
        b3 = submit("gan-qian-b3-relax", "hf_gan_qian_endpoint.slurm", {"ENDPOINT": "b3"})
        b1 = submit("gan-qian-b1-relax", "hf_gan_qian_endpoint.slurm", {"ENDPOINT": "b1"})
        hexa = submit("gan-qian-hex-check", "hf_gan_qian_validate.slurm", {"CASE_NAME": "b4_b1_hexagonal"})
        submit("gan-qian-hex-vcneb", "hf_gan_qian_vcneb.slurm", {"CASE_NAME": "b4_b1_hexagonal"}, (hexa,))
        diag = submit("gan-qian-b3-check", "hf_gan_qian_validate.slurm", {"CASE_NAME": "b3_b1_diagonal"}, (b3, b1))
        submit("gan-qian-b3-vcneb", "hf_gan_qian_vcneb.slurm", {"CASE_NAME": "b3_b1_diagonal"}, (diag,))
        payload["status"] = "submitted"
    except Exception as error:
        payload["status"] = "partial_submission_requires_inspection"
        payload["error"] = str(error)
        raise
    finally:
        atomic_json(manifest, payload)


if __name__ == "__main__":
    main()

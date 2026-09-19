"""Submit one evidence-gated GaN hex continuation from an immutable chain."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", required=True, type=Path)
    parser.add_argument("--trajectory", required=True, type=Path)
    parser.add_argument("--case-root", required=True, type=Path)
    parser.add_argument("--native-probe", required=True, type=Path)
    parser.add_argument("--canary-job", required=True)
    parser.add_argument("--run", action="store_true")
    args = parser.parse_args()
    if not args.canary_job.isdigit():
        parser.error("--canary-job must be a Slurm job id")
    repo = args.repo.resolve()
    manifest = repo / "hex_resume_launch.json"
    workdir = repo / "vcneb_hex_n29_w9_mpi32_poscar17"
    if manifest.exists() or workdir.exists():
        raise FileExistsError("launch/work directory exists; do not resubmit")
    required = [args.trajectory, args.native_probe]
    for endpoint in ("initial", "final"):
        required.extend([args.case_root / "input" / endpoint / name for name in ("CONTCAR","INCAR","KPOINTS","POTCAR")])
        required.append(args.case_root / f"static_{endpoint}" / "vasp_static_summary.json")
    for path in required:
        if not path.is_file():
            raise FileNotFoundError(path)
    payload = {
        "status": "prepared_not_submitted", "canary_afterok": args.canary_job,
        "trajectory": str(args.trajectory.resolve()), "trajectory_sha256": digest(args.trajectory),
        "native_probe": str(args.native_probe.resolve()), "native_probe_sha256": digest(args.native_probe),
        "workdir": str(workdir), "resources": {"nodes":3,"allocated_tasks":291,"workers":9,"mpi_per_worker":32,"n_interiors":27,"endpoint_workers":0},
        "physical_contract": {"vasp":"6.3.2","xc":"PBE","encut_eV":600,"kmesh":[8,8,6],"gamma":True,
                              "pressure_gpa":45.7,"fmax_eV_per_A":0.10,"CI":False,"isym":-1,"symprec":1e-4},
        "serialization_contract": {"policy":"exact_binary64_roundtrip","significant_decimal_digits":17},
        "source_sha256": {str(path.relative_to(repo)): digest(path)
                          for base in ("vcneb","scripts","cluster","examples")
                          for path in sorted((repo/base).rglob("*")) if path.suffix in (".py",".slurm")},
    }
    if args.run:
        (repo / "logs").mkdir(exist_ok=True)
        exported = ",".join((f"SOURCE_REPO={repo}", f"OLD_TRAJECTORY={args.trajectory.resolve()}",
                             f"CASE_ROOT={args.case_root.resolve()}", f"NATIVE_PROBE={args.native_probe.resolve()}",
                             f"WORKDIR={workdir}"))
        command = ["sbatch","--parsable",f"--dependency=afterok:{args.canary_job}",
                   f"--output={repo}/logs/%x-%j.out",f"--error={repo}/logs/%x-%j.err",
                   "--export=ALL,RUN_DFT=1,"+exported,str(repo/"cluster/hf_gan_hex_poscar17_resume.slurm")]
        result = subprocess.run(command,cwd=repo,text=True,capture_output=True,check=True)
        job_id=result.stdout.strip().split(";")[0]
        if not job_id.isdigit(): raise ValueError(f"invalid Slurm job id: {job_id!r}")
        payload.update(status="submitted",job_id=job_id,command=command)
    with manifest.open("x",encoding="utf-8") as handle:
        json.dump(payload,handle,indent=2); handle.write("\n")
    print(json.dumps({k:payload[k] for k in ("status","canary_afterok")}|({"job_id":payload["job_id"]} if "job_id" in payload else {})))


if __name__ == "__main__": main()

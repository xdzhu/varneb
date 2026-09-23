#!/bin/bash
# Pressure-consistent GaN B4/B1 endpoint BFGS for the 45.7 GPa VCNEB path.
# One Slurm task owns eight CPUs; Intel MPI's mpirun starts eight ABACUS ranks.
#SBATCH -J varneb-gan-aba-end
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=08:00:00
#SBATCH -o /public/home/iai806/abacus/agent-runs/20260921-varneb-material/slurm-%x-%j.out
#SBATCH -e /public/home/iai806/abacus/agent-runs/20260921-varneb-material/slurm-%x-%j.err

set -euo pipefail
[[ "${RUN_DFT:-0}" == 1 ]] || { echo "set RUN_DFT=1" >&2; exit 2; }
phase=${PHASE:?set PHASE=B4 or B1}
case "${phase}" in B4|B1) ;; *) echo "invalid PHASE=${phase}" >&2; exit 2;; esac
repo=${REPO:-/public/home/iai806/abacus/agent-runs/20260921-varneb-material/varneb_abacus_marker_source}
endpoint_dir=${ENDPOINT_DIR:-/public/home/iai806/abacus/agent-runs/20260921-varneb-material/cases/gan/endpoints}
workdir=${WORKDIR:?set a fresh phase-specific WORKDIR}
python=${PYTHON:-/public/home/iai806/.conda/envs/icu/bin/python}
abacus=${ABACUS_BIN:-/public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/app/bin/abacus}
profile=${PROFILE:-examples/material_profiles/abacus_gan_pbe_dzp_100ry_hf.json}
[[ "${SLURM_JOB_NUM_NODES:-0}" == 1 && "${SLURM_NTASKS:-0}" == 1 \
   && "${SLURM_CPUS_PER_TASK:-0}" -ge 8 ]] || {
  echo "requires one Slurm task owning at least eight CPUs" >&2; exit 2;
}
[[ ! -e "${workdir}" ]] || { echo "refusing existing WORKDIR=${workdir}" >&2; exit 2; }
[[ -s "${endpoint_dir}/${phase}_POSCAR" && -s "${repo}/${profile}" ]] || {
  echo "missing endpoint or ABACUS material profile" >&2; exit 2;
}

module purge
source /public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/scripts/env.sh
unset I_MPI_PMI_LIBRARY I_MPI_HYDRA_BOOTSTRAP_EXEC_EXTRA_ARGS I_MPI_HYDRA_BOOTSTRAP
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONUNBUFFERED=1
cd "${repo}"
printf 'phase=%s pressure_GPa=45.7 job=%s node=%s\n' \
  "${phase}" "${SLURM_JOB_ID}" "${SLURM_JOB_NODELIST}"
sha256sum "${endpoint_dir}/${phase}_POSCAR" "${profile}" "${abacus}"

exec "${python}" scripts/relax_ase_endpoint.py \
  --structure "${endpoint_dir}/${phase}_POSCAR" --workdir "${workdir}" \
  --factory vcneb.abacus:make_ase_abacus_factory --parameters "${profile}" \
  --command "mpirun -np 8 ${abacus}" \
  --pressure-gpa 45.7 --fmax 0.02 --steps 150 --maxstep 0.05

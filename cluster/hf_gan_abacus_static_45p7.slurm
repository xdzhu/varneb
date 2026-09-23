#!/bin/bash
# Independent frozen-endpoint ABACUS static evaluations at 45.7 GPa.
# Stress gate checks residual stress relative to -45.7 GPa, not absolute stress.
#SBATCH -J varneb-gan-aba-static
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00
#SBATCH -o /public/home/iai806/abacus/agent-runs/20260921-varneb-material/slurm-%x-%j.out
#SBATCH -e /public/home/iai806/abacus/agent-runs/20260921-varneb-material/slurm-%x-%j.err

set -euo pipefail
[[ "${RUN_DFT:-0}" == 1 ]] || { echo "set RUN_DFT=1" >&2; exit 2; }
repo=${REPO:-/public/home/iai806/abacus/agent-runs/20260921-varneb-material/varneb_abacus_45p7_source}
endpoint_root=${ENDPOINT_ROOT:-/public/home/iai806/abacus/agent-runs/20260921-varneb-material/cases/gan/abacus_endpoint_45p7_mpi8_20260922}
workdir=${WORKDIR:?set fresh static WORKDIR}
python=${PYTHON:-/public/home/iai806/.conda/envs/icu/bin/python}
abacus=${ABACUS_BIN:-/public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/app/bin/abacus}
profile=examples/material_profiles/abacus_gan_pbe_dzp_100ry_hf.json
[[ "${SLURM_JOB_NUM_NODES:-0}" == 1 && "${SLURM_NTASKS:-0}" == 1 \
   && "${SLURM_CPUS_PER_TASK:-0}" -ge 8 ]] || {
  echo "requires one Slurm task owning at least eight CPUs" >&2; exit 2;
}
[[ ! -e "${workdir}" ]] || { echo "refusing existing WORKDIR=${workdir}" >&2; exit 2; }
[[ -s "${endpoint_root}/B4/CONTCAR" && -s "${endpoint_root}/B1/CONTCAR" \
   && -s "${repo}/${profile}" ]] || { echo "missing optimized endpoints or profile" >&2; exit 2; }

module purge
source /public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/scripts/env.sh
unset I_MPI_PMI_LIBRARY I_MPI_HYDRA_BOOTSTRAP_EXEC_EXTRA_ARGS I_MPI_HYDRA_BOOTSTRAP
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONUNBUFFERED=1
cd "${repo}"
sha256sum "${endpoint_root}/B4/CONTCAR" "${endpoint_root}/B1/CONTCAR" "${profile}"

"${python}" examples/run_vcneb_ase.py \
  --initial "${endpoint_root}/B4/CONTCAR" \
  --final "${endpoint_root}/B1/CONTCAR" \
  --workdir "${workdir}" --n-images 29 --pressure-gpa 45.7 \
  --factory vcneb.abacus:make_ase_abacus_factory --parameters "${profile}" \
  --command "mpirun -np 8 ${abacus}" \
  --static-only --mic --cell-interpolation log_strain --mapping auto \
  --align-translation --minimum-distance 1.4 --maximum-deformation 0.50

# 2 kbar = 0.2 GPa, 0.44% of the target pressure.  This bounds the
# pressure-baseline error well below the GaN path barrier scale.
"${python}" scripts/validate_ase_static_gate.py \
  "${workdir}/ase_static_summary.json" \
  --pressure-gpa 45.7 --fmax 0.10 --stress-kbar 2.0 \
  --output "${workdir}/endpoint_static_gate.json"

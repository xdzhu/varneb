#!/bin/bash
# Real, fixed-endpoint QE SCF baseline for the BTO cross-calculator case.
# It consumes one 32-MPI allocation and evaluates only endpoint 00. It is a
# cutoff/SCF baseline, not an NEB step and not a path calculation. Use a
# distinct WORKDIR for each cutoff point and retain every result.
#SBATCH -J varneb-bto-qe-static
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH -o /public/home/iai806/abacus/agent-runs/20260912-vcneb-hf/slurm-%x-%j.out
#SBATCH -e /public/home/iai806/abacus/agent-runs/20260912-vcneb-hf/slurm-%x-%j.err

set -euo pipefail

repo=${REPO:-/public/home/iai806/abacus/agent-runs/20260912-vcneb-hf}
python=${PYTHON:-/public/home/iai806/.conda/envs/icu/bin/python}
pseudo_dir=${ESPRESSO_PSEUDO:-}
pp_manifest=${QE_PP_MANIFEST:-}
endpoint_tag=${ENDPOINT_TAG:-pbe100_dzp10au}
ecutwfc=${ECUTWFC:-100}
ecutrho=${ECUTRHO:-600}
scf_thr=${SCF_THR:-1e-8}
workdir=${WORKDIR:-validation/batio3_cubic_to_tetragonal/qe_static_ecut${ecutwfc}_rho${ecutrho}}

if [[ "${RUN_DFT:-0}" != 1 ]]; then
  echo "Refusing a 32-rank allocation without RUN_DFT=1" >&2
  exit 2
fi
[[ -n "${pseudo_dir}" ]] || { echo "Set ESPRESSO_PSEUDO to the reviewed UPF directory" >&2; exit 2; }
[[ -n "${pp_manifest}" ]] || { echo "Set QE_PP_MANIFEST to an approved MD5-pinned UPF manifest" >&2; exit 2; }

module purge
: "${QE_ENV_SCRIPT:?Set QE_ENV_SCRIPT to the reviewed QE environment script}"
source "${QE_ENV_SCRIPT}"
export OMP_NUM_THREADS=1
cd "${repo}"

exec "${python}" examples/run_vcneb_qe.py \
  --initial "validation/batio3_cubic_to_tetragonal/endpoint_relax_tetragonal_${endpoint_tag}/CONTCAR" \
  --final "validation/batio3_cubic_to_tetragonal/endpoint_relax_cubic_${endpoint_tag}/CONTCAR" \
  --workdir "${workdir}" --n-images 7 --fmax 0.10 --steps 0 --k 0.2 \
  --image-workers 0 --command "srun --exclusive --nodes=1 --ntasks=32 pw.x" \
  --pseudo-dir "${pseudo_dir}" --pp-manifest "${pp_manifest}" \
  --ecutwfc "${ecutwfc}" --ecutrho "${ecutrho}" --scf-thr "${scf_thr}" \
  --kpts 4 4 4 --no-climb --mic --cell-interpolation log_strain --mapping auto \
  --align-translation --minimum-distance 1.6 --maximum-deformation 0.10 --static-only

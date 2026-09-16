#!/bin/bash
# BTO T->C QE VCNEB: one manager and five independent 32-MPI interior workers.
# This is the production allocation. Use the separate one-rank preflight
# template before submitting it; RUN_DFT=1 is deliberately mandatory here.
#SBATCH -J varneb-bto-qe
#SBATCH -p hfacnormal01
#SBATCH -N 2
#SBATCH --ntasks=160
#SBATCH --cpus-per-task=1
#SBATCH -o /public/home/iai806/abacus/agent-runs/20260912-vcneb-hf/slurm-%x-%j.out
#SBATCH -e /public/home/iai806/abacus/agent-runs/20260912-vcneb-hf/slurm-%x-%j.err

set -euo pipefail

repo=${REPO:-/public/home/iai806/abacus/agent-runs/20260912-vcneb-hf}
python=${PYTHON:-/public/home/iai806/.conda/envs/icu/bin/python}
qe_bin=${QE_BIN:-pw.x}
pseudo_dir=${ESPRESSO_PSEUDO:-}
n_images=${N_IMAGES:-7}
n_interior=$((n_images - 2))
image_workers=${IMAGE_WORKERS:-${n_interior}}
image_mpi=${IMAGE_MPI:-32}
steps=${STEPS:-300}
fmax=${FMAX:-0.10}
spring=${SPRING:-0.2}
optimizer=${OPTIMIZER:-FIRE}
maxstep=${MAXSTEP:-0.02}
ecutwfc=${ECUTWFC:-100}
ecutrho=${ECUTRHO:-600}
scf_thr=${SCF_THR:-1e-8}
kpt_x=${KPT_X:-4}
kpt_y=${KPT_Y:-4}
kpt_z=${KPT_Z:-4}
pp_ba=${PP_BA:-Ba.upf}
pp_ti=${PP_TI:-Ti.upf}
pp_o=${PP_O:-O.upf}
pp_manifest=${QE_PP_MANIFEST:-}
endpoint_identity_gate=${ENDPOINT_IDENTITY_GATE:-}
endpoint_tag=${ENDPOINT_TAG:-pbe100_dzp10au}
workdir=${WORKDIR:-validation/batio3_cubic_to_tetragonal/vcneb_tetragonal_to_cubic_n${n_images}_${optimizer,,}_qe_pbe}

if (( n_images != 7 )); then
  echo "This distributed BTO QE validation template is intentionally fixed at 7 total images" >&2
  exit 2
fi
if (( image_workers != n_interior || image_workers * image_mpi != 160 )); then
  echo "BTO QE validation requires five 32-MPI interior workers (5 x 32 = 160)" >&2
  exit 2
fi
if (( SLURM_NTASKS < image_workers * image_mpi )); then
  echo "allocated tasks are insufficient for five 32-MPI workers" >&2
  exit 2
fi
if [[ "${RUN_DFT:-0}" != 1 ]]; then
  echo "Refusing a 160-rank allocation without RUN_DFT=1; use hf_batio3_vcneb_qe_preflight.slurm first" >&2
  exit 2
fi
if [[ -z "${pseudo_dir}" ]]; then
  echo "Set ESPRESSO_PSEUDO to the reviewed UPF directory" >&2
  exit 2
fi
if [[ -z "${pp_manifest}" ]]; then
  echo "Set QE_PP_MANIFEST to an approved MD5-pinned UPF manifest" >&2
  exit 2
fi

module purge
: "${QE_ENV_SCRIPT:?Set QE_ENV_SCRIPT to the reviewed QE environment script}"
: "${ENDPOINT_IDENTITY_GATE:?Set ENDPOINT_IDENTITY_GATE to a passing endpoint_identity_gate.json}"
source "${QE_ENV_SCRIPT}"
export OMP_NUM_THREADS=1
export VCNEB_GIT_REVISION="${VCNEB_GIT_REVISION:-remote-sync-unknown}"
cd "${repo}"
"${python}" scripts/validate_endpoint_identity_gate.py --gate "${endpoint_identity_gate}"

exec "${python}" examples/run_vcneb_qe.py \
  --initial "validation/batio3_cubic_to_tetragonal/endpoint_relax_tetragonal_${endpoint_tag}/CONTCAR" \
  --final "validation/batio3_cubic_to_tetragonal/endpoint_relax_cubic_${endpoint_tag}/CONTCAR" \
  --workdir "${workdir}" --n-images "${n_images}" --fmax "${fmax}" \
  --steps "${steps}" --k "${spring}" --optimizer "${optimizer}" --maxstep "${maxstep}" \
  --image-workers "${image_workers}" \
  --command "srun --exclusive --nodes=1 --ntasks=${image_mpi} --ntasks-per-node=${image_mpi} ${qe_bin}" \
  --pseudo-dir "${pseudo_dir}" --pp-manifest "${pp_manifest}" \
  --ecutwfc "${ecutwfc}" --ecutrho "${ecutrho}" --scf-thr "${scf_thr}" \
  --kpts "${kpt_x}" "${kpt_y}" "${kpt_z}" \
  --no-climb --mic --cell-interpolation log_strain --mapping auto --align-translation \
  --minimum-distance 1.6 --maximum-deformation 0.10

#!/bin/bash
# Two model predictions were committed before these new DFT calculations.
# The original BTO 100 Ry / 10 au DZP / 4x4x4 electronic contract is unchanged.
#SBATCH -J varneb-bto-qy0-even-validate
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --array=0-1%2
#SBATCH --time=05:00:00
#SBATCH -o /public/home/iai806/abacus/agent-runs/20260928-varneb-bto-qy0-pilot/slurm-even-holdout-%A_%a.out
#SBATCH -e /public/home/iai806/abacus/agent-runs/20260928-varneb-bto-qy0-pilot/slurm-even-holdout-%A_%a.err

set -euo pipefail
[[ "${RUN_DFT:-0}" == 1 ]] || { echo 'set RUN_DFT=1 explicitly' >&2; exit 2; }
[[ "${SLURM_JOB_PARTITION:-}" == hfacnormal01 && "${SLURM_JOB_NUM_NODES:-0}" == 1 \
   && "${SLURM_NTASKS:-0}" == 1 && "${SLURM_CPUS_PER_TASK:-0}" -ge 32 ]] || {
  echo 'requires one 32-CPU Slurm task on hfacnormal01' >&2
  exit 2
}

root=/public/home/iai806/abacus/agent-runs/20260928-varneb-bto-qy0-pilot
code=${root}/code-holdouts
inputs=${root}/inputs
plan=${root}/bto_qy0_even_mode_two_holdout_plan_20260928.json
python=/public/home/iai806/.conda/envs/icu/bin/python
abacus=/public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/app/bin/abacus
pseudo_dir=/public/home/iai806/abacus/PSEUDO/ABACUS-orbitals/Dojo-NC-FR/Pseudopotential
basis_dir=/public/home/iai806/abacus/PSEUDO/ABACUS-orbitals/Dojo-NC-FR/Orb-DZP-10au
env_script=/public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/scripts/env.sh

case "${SLURM_ARRAY_TASK_ID:-}" in
  0) point=q045_q025; qz=0.45; qx=0.25; preflight_hash=85a5b4357e230597e128b8afa558e294aa0f0c0ac78ce40aa10cc3e88eb8b6f3 ;;
  1) point=q105_q025; qz=1.05; qx=0.25; preflight_hash=79f2dbdbc6de831b9f8c8c13eb971e4b5fde957b96d0ab6f1f0dd370ed82e396 ;;
  *) echo 'unexpected holdout array index' >&2; exit 2 ;;
esac
preflight=${inputs}/bto_qy0_${point}_model_preflight_2026-09-28.json
workdir=${root}/run-model-holdout-${point}
result=${point}_model_holdout_result.json

[[ ! -e "${workdir}" && -s "${plan}" ]] || {
  echo 'frozen plan required and holdout workdir must be unused' >&2
  exit 2
}
echo "${preflight_hash}  ${preflight}" | sha256sum -c -
echo "0cd1a509cec79dffbc353548d2d7fb13975966886fc6bdb14cb22841a3ecc07d  ${plan}" | sha256sum -c -
echo "23fa2b5230bf1904c759d79f6b373be855e5db20e6cf9183bdc2cda3b8eaca4b  ${code}/examples/run_bto_q1q2_conditional_abacus.py" | sha256sum -c -
echo "55bcf70c3af12b4a5d4eb445ae30f8461c07f18890afc34c05ee3b084b64164a  ${code}/vcneb/mode_surface.py" | sha256sum -c -
"${python}" - "${plan}" "${SLURM_ARRAY_TASK_ID}" "${qz}" "${qx}" <<'PY'
import json
from pathlib import Path
import sys

plan = json.loads(Path(sys.argv[1]).read_text())
index, qz, qx = int(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
point = plan['new_holdouts'][index]
if (plan.get('status') != 'two_model_predictions_frozen_before_new_DFT'
        or plan.get('energy_absolute_error_gate_meV_per_BTO') != 2.0
        or plan.get('full_structure_seed_error_gate_sqrt_amu_A') != 0.10
        or point['q_sqrt_amu_A'] != [qz, qx]):
    raise SystemExit('model holdout position or frozen gate changed')
PY

module purge
source "${env_script}"
unset I_MPI_PMI_LIBRARY I_MPI_HYDRA_BOOTSTRAP_EXEC_EXTRA_ARGS I_MPI_HYDRA_BOOTSTRAP
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
cd "${code}"
"${python}" examples/run_bto_q1q2_conditional_abacus.py \
  --preflight "${preflight}" \
  --report "${inputs}/bto_t_to_c_q1q2_strain_gradient_preflight_2026-09-25.json" \
  --reference "${inputs}/cubic_CONTCAR" \
  --force-constants "${inputs}/bto_cubic_gamma_force_constants.npz" \
  --phonopy-eigenpairs "${inputs}/bto_cubic_phonopy_gamma_eigenpairs.npz" \
  --gamma-provenance "${inputs}/bto_cubic_gamma_phonon_provenance.json" \
  --force-sets "${inputs}/bto_cubic_gamma_FORCE_SETS" \
  --eigenpairs-provenance "${inputs}/bto_cubic_phonopy_gamma_eigenpairs_provenance.json" \
  --grid-result-manifest "${inputs}/result_manifest_hf_27777454.json" \
  --workdir "${workdir}" \
  --pseudo-dir "${pseudo_dir}" --basis-dir "${basis_dir}" \
  --abacus-bin "${abacus}" --mpi-ranks 32 \
  --gradient-tolerance 0.003 --orthogonal-amplitude-bound 4.0 \
  --max-iterations 64 --max-new-evaluations 120 \
  --max-absolute-stress-kbar 2.0 \
  --result-name "${result}" --run-dft

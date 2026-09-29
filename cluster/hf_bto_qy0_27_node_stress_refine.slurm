#!/bin/bash
# Continue only the two previously gradient-converged Qy=0 cases whose
# raw stress was 2.45 kbar > the predeclared 2.0-kbar gate. No DFT settings change.
#SBATCH -J varneb-bto-qy0-stress
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --array=0-1%2
#SBATCH --time=03:00:00
#SBATCH -o /public/home/iai806/abacus/agent-runs/20260929-varneb-bto-qy0-grid27/slurm-stress-%A_%a.out
#SBATCH -e /public/home/iai806/abacus/agent-runs/20260929-varneb-bto-qy0-grid27/slurm-stress-%A_%a.err

set -euo pipefail
[[ "${RUN_DFT:-0}" == 1 ]] || exit 2
[[ "${SLURM_JOB_PARTITION:-}" == hfacnormal01 && "${SLURM_JOB_NUM_NODES:-0}" == 1 \
   && "${SLURM_NTASKS:-0}" == 1 && "${SLURM_CPUS_PER_TASK:-0}" -ge 32 ]] || exit 2
root=/public/home/iai806/abacus/agent-runs/20260929-varneb-bto-qy0-grid27
prior=/public/home/iai806/abacus/agent-runs/20260928-varneb-bto-qy0-pilot
code=${prior}/code-holdouts
inputs=${prior}/inputs
python=/public/home/iai806/.conda/envs/icu/bin/python
abacus=/public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/app/bin/abacus
pseudo_dir=/public/home/iai806/abacus/PSEUDO/ABACUS-orbitals/Dojo-NC-FR/Pseudopotential
basis_dir=/public/home/iai806/abacus/PSEUDO/ABACUS-orbitals/Dojo-NC-FR/Orb-DZP-10au
case "${SLURM_ARRAY_TASK_ID:-}" in
  0) point=qz01_qx02; prior_hash=e0de115364c965f939cfc20239ee2ded0efe0b9d7c526958c5736c990385613a ;;
  1) point=qz02_qx01; prior_hash=66406d0af8c5e5b00a7477a62438d440023f9310bd31821783b3d17fcaadd279 ;;
  *) exit 2 ;;
esac
workdir=${root}/run-${point}
preflight=${root}/preflights/${point}.json
prior_result=${workdir}/${point}_result.json
result=${point}_stress_refined_result.json
[[ -s "${prior_result}" && -s "${preflight}" && ! -e "${workdir}/${result}" ]] || exit 3
printf '%s  %s\n' "${prior_hash}" "${prior_result}" | sha256sum -c -
printf '%s  %s\n' '23fa2b5230bf1904c759d79f6b373be855e5db20e6cf9183bdc2cda3b8eaca4b' \
  "${code}/examples/run_bto_q1q2_conditional_abacus.py" | sha256sum -c -
"${python}" - "${prior_result}" "${preflight}" <<'PY'
import hashlib
import json
from pathlib import Path
import sys

summary = json.loads(Path(sys.argv[1]).read_text())
preflight = Path(sys.argv[2])
if (summary.get('status') != 'orthogonal_gradient_converged_stress_target_failed_curvature_unchecked'
        or summary.get('preflight_sha256') != hashlib.sha256(preflight.read_bytes()).hexdigest()
        or summary.get('stress_target_kbar') != 2.0
        or not 2.0 < summary.get('maximum_absolute_stress_kbar', 0) < 2.5):
    raise SystemExit('unexpected failed-stress warm-start contract')
PY
module purge
source /public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/scripts/env.sh
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
  --workdir "${workdir}" --resume --start-from-result "${prior_result}" \
  --pseudo-dir "${pseudo_dir}" --basis-dir "${basis_dir}" \
  --abacus-bin "${abacus}" --mpi-ranks 32 \
  --gradient-tolerance 0.0015 --orthogonal-amplitude-bound 4.0 \
  --max-iterations 64 --max-new-evaluations 120 \
  --max-absolute-stress-kbar 2.0 \
  --result-name "${result}" --run-dft

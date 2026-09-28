#!/bin/bash
# Resume two raw-audited low-Q BTO points from their complete caches.
# Tighten only the optimizer gradient; retain 2-kbar stress and all DFT inputs.
#SBATCH -J varneb-bto-qy0-lowq
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --array=0-1%2
#SBATCH --time=03:00:00
#SBATCH -o /public/home/iai806/abacus/agent-runs/20260928-varneb-bto-qy0-pilot/slurm-lowq-%A_%a.out
#SBATCH -e /public/home/iai806/abacus/agent-runs/20260928-varneb-bto-qy0-pilot/slurm-lowq-%A_%a.err

set -euo pipefail
[[ "${RUN_DFT:-0}" == 1 ]] || { echo 'set RUN_DFT=1 explicitly' >&2; exit 2; }
[[ "${SLURM_JOB_PARTITION:-}" == hfacnormal01 && "${SLURM_JOB_NUM_NODES:-0}" == 1 \
   && "${SLURM_NTASKS:-0}" == 1 && "${SLURM_CPUS_PER_TASK:-0}" -ge 32 ]] || {
  echo 'requires one 32-CPU Slurm task on hfacnormal01' >&2
  exit 2
}

root=/public/home/iai806/abacus/agent-runs/20260928-varneb-bto-qy0-pilot
code=${root}/code
inputs=${root}/inputs
python=/public/home/iai806/.conda/envs/icu/bin/python
abacus=/public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/app/bin/abacus
pseudo_dir=/public/home/iai806/abacus/PSEUDO/ABACUS-orbitals/Dojo-NC-FR/Pseudopotential
basis_dir=/public/home/iai806/abacus/PSEUDO/ABACUS-orbitals/Dojo-NC-FR/Orb-DZP-10au
env_script=/public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/scripts/env.sh

case "${SLURM_ARRAY_TASK_ID:-}" in
  0)
    point=q030_q000
    prior_hash=e29621ac255b70e980a290a98e313a9bfda0b046d9419729cff9d8bb7354aa05
    audit_hash=4f2f5a4a9d93a981ce5f943a16cda1d2e3b875cce3bb0d1ae33151d727bc4bff
    initial_array_index=1
    ;;
  1)
    point=q030_q030
    prior_hash=3f286ddce30961f3bd6e340a58606597ea68735c088b855522e2257a5d773351
    audit_hash=e006a07b15a952b374ab4665b7eab3061720135f88143a4adfef83280caca4bc
    initial_array_index=2
    ;;
  *) echo 'unexpected array index' >&2; exit 2 ;;
esac
workdir=${root}/run-${point}
prior=${workdir}/${point}_restricted_result.json
audit=${root}/audit-${point}-27794649_${initial_array_index}.json
preflight=${inputs}/bto_qy0_${point}_preflight_2026-09-28.json
result=${point}_restricted_tight_result.json
[[ -s "${workdir}/contract.json" && -s "${prior}" && -s "${audit}" \
   && ! -e "${workdir}/${result}" ]] || {
  echo 'requires audited prior result and unused continuation result path' >&2
  exit 2
}
echo "${prior_hash}  ${prior}" | sha256sum -c -
echo "${audit_hash}  ${audit}" | sha256sum -c -
echo "4db5823766b88c4e9e103827bee198e51a9a44954bfb532ade6ff5cdcf6c3f1a  ${code}/examples/run_bto_q1q2_conditional_abacus.py" | sha256sum -c -
"${python}" - "${prior}" "${audit}" <<'PY'
import json
from pathlib import Path
import sys

prior, audit = (json.loads(Path(path).read_text()) for path in sys.argv[1:])
if (prior.get('status') != 'orthogonal_gradient_converged_stress_target_failed_curvature_unchecked'
        or prior.get('stress_target_kbar') != 2.0
        or prior.get('stress_target_passed') is not False
        or audit.get('status') != 'verified_gradient_stationary_stress_target_failed_curvature_unchecked'
        or audit.get('n_individually_audited_DFT_points') != 7
        or audit.get('final_maximum_absolute_stress_kbar', 0) <= 2.0
        or abs(audit.get('final_third_soft_y_amplitude_sqrt_amu_A', 1)) > 1e-8):
    raise SystemExit('the source is not the audited stress-only failure')
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
  --gradient-tolerance 0.001 --orthogonal-amplitude-bound 3.0 \
  --max-iterations 32 --max-new-evaluations 60 \
  --max-absolute-stress-kbar 2.0 \
  --start-from-result "${prior}" \
  --result-name "${result}" --resume --run-dft

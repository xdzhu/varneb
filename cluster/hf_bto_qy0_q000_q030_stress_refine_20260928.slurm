#!/bin/bash
# Continue only the audited Q_y=0 point whose 2-kbar stress gate narrowly failed.
# The ABACUS Hamiltonian, fixed modes and previous complete evaluations stay unchanged.
#SBATCH -J varneb-bto-qy0-stress
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=03:00:00
#SBATCH -o /public/home/iai806/abacus/agent-runs/20260928-varneb-bto-qy0-pilot/slurm-stress-%j.out
#SBATCH -e /public/home/iai806/abacus/agent-runs/20260928-varneb-bto-qy0-pilot/slurm-stress-%j.err

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
workdir=${root}/run-q000_q030
prior=${workdir}/q000_q030_restricted_result.json
audit=${root}/audit-q000_q030-27794649_0-v2.json
preflight=${inputs}/bto_qy0_q000_q030_preflight_2026-09-28.json
python=/public/home/iai806/.conda/envs/icu/bin/python
abacus=/public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/app/bin/abacus
pseudo_dir=/public/home/iai806/abacus/PSEUDO/ABACUS-orbitals/Dojo-NC-FR/Pseudopotential
basis_dir=/public/home/iai806/abacus/PSEUDO/ABACUS-orbitals/Dojo-NC-FR/Orb-DZP-10au
env_script=/public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/scripts/env.sh

[[ -s "${workdir}/contract.json" && -s "${prior}" && -s "${audit}" \
   && ! -e "${workdir}/q000_q030_restricted_tight_result.json" ]] || {
  echo 'requires audited prior result and unused continuation result path' >&2
  exit 2
}
echo "3527b27df926b0aaacc138265679df0c9e96ffc461baf854d8cbacf6685b5f73  ${prior}" | sha256sum -c -
echo "1037844b5b2fb435a041f69a591257e23ba612e20acb28d96a89c5458293b3b6  ${audit}" | sha256sum -c -
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
  --result-name q000_q030_restricted_tight_result.json --resume --run-dft

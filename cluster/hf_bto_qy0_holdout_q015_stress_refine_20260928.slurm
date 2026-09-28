#!/bin/bash
# Same Hamiltonian and predeclared 2-kbar physical gate; only optimizer-coordinate
# stopping tolerance is tightened after a one-point raw-audited stress miss.
#SBATCH -J varneb-bto-qy0-h015r
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=05:00:00
#SBATCH -o /public/home/iai806/abacus/agent-runs/20260928-varneb-bto-qy0-pilot/slurm-h015-refine-%j.out
#SBATCH -e /public/home/iai806/abacus/agent-runs/20260928-varneb-bto-qy0-pilot/slurm-h015-refine-%j.err

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
workdir=${root}/run-holdout-q015_q015
prior=${workdir}/q015_q015_holdout_result.json
audit=${root}/audit-holdout-q015_q015-27794991_0.json
preflight=${inputs}/bto_qy0_q015_q015_holdout_preflight_2026-09-28.json
result=q015_q015_holdout_tight_result.json
python=/public/home/iai806/.conda/envs/icu/bin/python
abacus=/public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/app/bin/abacus
pseudo_dir=/public/home/iai806/abacus/PSEUDO/ABACUS-orbitals/Dojo-NC-FR/Pseudopotential
basis_dir=/public/home/iai806/abacus/PSEUDO/ABACUS-orbitals/Dojo-NC-FR/Orb-DZP-10au
env_script=/public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/scripts/env.sh

[[ -s "${workdir}/contract.json" && ! -e "${workdir}/${result}" ]] || {
  echo 'requires prior evaluator cache and unused continuation result' >&2
  exit 2
}
echo "cb7dfc09a63e852e659c2929c4266a2ef246126b9995292520a128523657797d  ${prior}" | sha256sum -c -
echo "3578c665d3aaca3cbc040ad8d088b6ea16d478b44e3b01db20cb5d8ea5f24b3a  ${audit}" | sha256sum -c -
echo "787a495b2119dfce8ccfbfe3be2fb21d969b40412139f74b9f9993f5b34ebe63  ${preflight}" | sha256sum -c -
echo "23fa2b5230bf1904c759d79f6b373be855e5db20e6cf9183bdc2cda3b8eaca4b  ${code}/examples/run_bto_q1q2_conditional_abacus.py" | sha256sum -c -
"${python}" - "${prior}" "${audit}" <<'PY'
import json
from pathlib import Path
import sys

prior, audit = (json.loads(Path(path).read_text()) for path in sys.argv[1:])
if (prior.get('status') != 'orthogonal_gradient_converged_stress_target_failed_curvature_unchecked'
        or prior.get('stress_target_kbar') != 2.0
        or prior.get('stress_target_passed') is not False
        or audit.get('status') != 'verified_gradient_stationary_stress_target_failed_curvature_unchecked'
        or audit.get('n_individually_audited_DFT_points') != 1
        or audit.get('final_maximum_absolute_stress_kbar', 0) <= 2.0
        or abs(audit.get('final_third_soft_y_amplitude_sqrt_amu_A', 1)) > 1e-8):
    raise SystemExit('source is not the audited stress-only holdout miss')
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
  --gradient-tolerance 0.001 --orthogonal-amplitude-bound 4.0 \
  --max-iterations 64 --max-new-evaluations 60 \
  --max-absolute-stress-kbar 2.0 \
  --start-from-result "${prior}" \
  --result-name "${result}" --resume --run-dft

#!/bin/bash
# Fifteen missing points of the 9x3 Qy=0 restricted BTO sheet.
# Same PBE/100-Ry/10-au-DZP/4x4x4 contract and optimizer gates as prior holdouts.
#SBATCH -J varneb-bto-qy0-grid27
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --array=0-14%8
#SBATCH --time=05:00:00
#SBATCH -o /public/home/iai806/abacus/agent-runs/20260929-varneb-bto-qy0-grid27/slurm-%A_%a.out
#SBATCH -e /public/home/iai806/abacus/agent-runs/20260929-varneb-bto-qy0-grid27/slurm-%A_%a.err

set -euo pipefail
[[ "${RUN_DFT:-0}" == 1 ]] || { echo 'set RUN_DFT=1 explicitly' >&2; exit 2; }
[[ "${SLURM_JOB_PARTITION:-}" == hfacnormal01 && "${SLURM_JOB_NUM_NODES:-0}" == 1 \
   && "${SLURM_NTASKS:-0}" == 1 && "${SLURM_CPUS_PER_TASK:-0}" -ge 32 ]] || {
  echo 'requires one 32-CPU hfacnormal01 Slurm task' >&2; exit 2;
}

root=/public/home/iai806/abacus/agent-runs/20260929-varneb-bto-qy0-grid27
prior=/public/home/iai806/abacus/agent-runs/20260928-varneb-bto-qy0-pilot
code=${prior}/code-holdouts
inputs=${prior}/inputs
plan=${root}/plan.json
python=/public/home/iai806/.conda/envs/icu/bin/python
abacus=/public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/app/bin/abacus
pseudo_dir=/public/home/iai806/abacus/PSEUDO/ABACUS-orbitals/Dojo-NC-FR/Pseudopotential
basis_dir=/public/home/iai806/abacus/PSEUDO/ABACUS-orbitals/Dojo-NC-FR/Orb-DZP-10au
env_script=/public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/scripts/env.sh

echo '4449b9be73ff2ea1e3d025d1f8bb83194e5b79a38bc23fc882b292369439e2a1  '${plan} | sha256sum -c -
echo '23fa2b5230bf1904c759d79f6b373be855e5db20e6cf9183bdc2cda3b8eaca4b  '${code}/examples/run_bto_q1q2_conditional_abacus.py | sha256sum -c -
index=${SLURM_ARRAY_TASK_ID}
mapfile -t point_info < <("${python}" - "${plan}" "${index}" <<'PY'
import hashlib
import json
from pathlib import Path
import sys

path, index = Path(sys.argv[1]), int(sys.argv[2])
plan = json.loads(path.read_text())
if (plan.get('status') != '15_missing_Qy0_inputs_preflighted_no_DFT'
        or plan.get('n_new_conditional_points') != 15
        or plan.get('gradient_tolerance_eV_per_sqrt_amu_A') != 0.003
        or plan.get('maximum_absolute_stress_kbar') != 2.0):
    raise SystemExit('restricted-sheet plan changed')
row = plan['new'][index]
preflight = path.parent / 'preflights' / (row['name'] + '.json')
digest = hashlib.sha256(preflight.read_bytes()).hexdigest()
if digest != row['preflight_sha256']:
    raise SystemExit('restricted-sheet preflight hash changed')
print(row['name'])
print(str(preflight))
PY
)
[[ "${#point_info[@]}" == 2 ]] || { echo 'invalid point mapping' >&2; exit 2; }
point=${point_info[0]}
preflight=${point_info[1]}
workdir=${root}/run-${point}
[[ ! -e "${workdir}" ]] || { echo 'refusing to overwrite workdir' >&2; exit 2; }

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
  --result-name "${point}_result.json" --run-dft

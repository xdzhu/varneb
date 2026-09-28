#!/usr/bin/env bash
# Add only 22 never-calculated nodes to the 59-point frozen BTO soft/soft cut.
#SBATCH -J varneb-bto-soft22
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=01:00:00
#SBATCH -o /public/home/iai806/abacus/agent-runs/bto_transverse_missing22_20260928/slurm-%x-%j.out
#SBATCH -e /public/home/iai806/abacus/agent-runs/bto_transverse_missing22_20260928/slurm-%x-%j.err

set -euo pipefail
[[ "${RUN_DFT:-0}" == 1 ]] || exit 2
[[ "${SLURM_JOB_PARTITION:-}" == hfacnormal01 && "${SLURM_JOB_NUM_NODES:-0}" == 1 \
   && "${SLURM_NTASKS:-0}" == 1 && "${SLURM_CPUS_PER_TASK:-0}" -ge 32 ]] || exit 2

root=/public/home/iai806/abacus/agent-runs/bto_transverse_missing22_20260928
code=${root}/code
points=${root}/bto_transverse_missing22_20260928
workdir=${root}/run-missing22-001
preflight=${root}/input-preflight-001/result_manifest.json
python=/public/home/iai806/.conda/envs/icu/bin/python
abacus=/public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/app/bin/abacus
pseudo_dir=/public/home/iai806/abacus/PSEUDO/ABACUS-orbitals/Dojo-NC-FR/Pseudopotential
basis_dir=/public/home/iai806/abacus/PSEUDO/ABACUS-orbitals/Dojo-NC-FR/Orb-DZP-10au

[[ ! -e "${workdir}" && -s "${points}/manifest.json" && -s "${preflight}" ]] || exit 3
echo "187ef8f1cffa0c574ed83ec8a5f930b8035312c19ade28f186b544a419c79fba  ${points}/manifest.json" | sha256sum -c -
echo "26acc0c0dfea66aa932cf47d6b7fac48b90488d0bd8ef158a46584668a8d02a5  ${code}/examples/run_bto_q1q2_frozen_abacus.py" | sha256sum -c -
echo "9d12c8907abbca66bcff6e5deebbdef09b2fa33d61593442002b821a5db7d6f8  ${preflight}" | sha256sum -c -
grep -q '"status": "inputs_validated_no_dft"' "${preflight}"
[[ "$(find "${points}" -mindepth 2 -maxdepth 2 -name POSCAR | wc -l)" == 22 ]] || exit 4

module purge
source /public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/scripts/env.sh
unset I_MPI_PMI_LIBRARY I_MPI_HYDRA_BOOTSTRAP_EXEC_EXTRA_ARGS I_MPI_HYDRA_BOOTSTRAP
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "${code}"
"${python}" examples/run_bto_q1q2_frozen_abacus.py \
  --points-dir "${points}" --workdir "${workdir}" \
  --pseudo-dir "${pseudo_dir}" --basis-dir "${basis_dir}" \
  --abacus-bin "${abacus}" --mpi-ranks 32 --run-dft

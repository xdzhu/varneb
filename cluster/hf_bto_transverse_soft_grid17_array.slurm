#!/usr/bin/env bash
# Eight disjoint 26-point shards; 81 old audited nodes are reused, never rerun.
# Submit explicitly: sbatch --export=ALL,RUN_DFT=1 hf_bto_transverse_soft_grid17_array.slurm
#SBATCH -J varneb-bto-soft17
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --array=0-7%8
#SBATCH --time=02:00:00
#SBATCH -o /public/home/iai806/abacus/agent-runs/20260929-varneb-bto-soft-grid17/slurm-%A_%a.out
#SBATCH -e /public/home/iai806/abacus/agent-runs/20260929-varneb-bto-soft-grid17/slurm-%A_%a.err

set -euo pipefail
[[ "${RUN_DFT:-0}" == 1 ]] || { echo 'RUN_DFT=1 is required' >&2; exit 2; }
[[ "${SLURM_JOB_PARTITION:-}" == hfacnormal01 && "${SLURM_JOB_NUM_NODES:-0}" == 1 \
   && "${SLURM_NTASKS:-0}" == 1 && "${SLURM_CPUS_PER_TASK:-0}" -ge 32 ]] || {
  echo 'requires one 32-CPU Slurm task on hfacnormal01' >&2; exit 2;
}
[[ "${SLURM_ARRAY_TASK_ID:-}" =~ ^[0-7]$ ]] || exit 2
printf -v shard '%02d' "${SLURM_ARRAY_TASK_ID}"
root=/public/home/iai806/abacus/agent-runs/20260929-varneb-bto-soft-grid17
code=${root}/code
points=${root}/qparallel_qtransverse_frozen_grid17_increment_2026-09-29/shard-${shard}
preflight=${root}/input-preflight-${shard}/result_manifest.json
workdir=${root}/run-shard-${shard}-001
python=/public/home/iai806/.conda/envs/icu/bin/python
abacus=/public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/app/bin/abacus
pseudo_dir=/public/home/iai806/abacus/PSEUDO/ABACUS-orbitals/Dojo-NC-FR/Pseudopotential
basis_dir=/public/home/iai806/abacus/PSEUDO/ABACUS-orbitals/Dojo-NC-FR/Orb-DZP-10au

case "${shard}" in
  00) expected=03b39e0b8ff2b5a812d51ff69938b45a92bc56c3be567e1bc4c030fd6e310215 ;;
  01) expected=507f989fbd924aab1b6ed0f2e0a6e0556dd7b8c144b593f120c579ee2afc9379 ;;
  02) expected=346a11f4adc5b2800fd6ed4683d1532661f1365c76c867d9962777f06e2cd22e ;;
  03) expected=692d8566940224e259ca43a00985186fdd933bb0df2b98ad081d1fc2a8239e5e ;;
  04) expected=aabf39cd49daf1ea5e6a343376bfdcc52732787d694b40923df07ddb02096cf3 ;;
  05) expected=29ab8b44c6bca62dbdfb6878be7bb5af9ad199569c89988184fa1162986c2f05 ;;
  06) expected=dea19c4b6639d5685f0eade346d21007689cf3ac199714df22d03c28a029d1cc ;;
  07) expected=87acbaebfc23e717dcd3e5d9332873d53cfa14b3b7b85b18035154ea900df5b0 ;;
esac
[[ ! -e "${workdir}" && -s "${points}/manifest.json" && -s "${preflight}" \
   && -s "${code}/examples/run_bto_q1q2_frozen_abacus.py" ]] || exit 3
printf '%s  %s\n' "${expected}" "${points}/manifest.json" | sha256sum -c -
printf '%s  %s\n' '87a2a4e52504bea1eda45e1e65729e6065a592576de5c924ce1d149e73494b05' \
  "${code}/examples/run_bto_q1q2_frozen_abacus.py" | sha256sum -c -
[[ "$(find "${points}" -mindepth 2 -maxdepth 2 -name POSCAR | wc -l)" == 26 ]] || exit 3

"${python}" - "${preflight}" "${expected}" <<'PY'
import json
import sys
from pathlib import Path

report = json.loads(Path(sys.argv[1]).read_text())
expected = sys.argv[2]
p = report.get('calculator_parameters', {})
assert report.get('status') == 'inputs_validated_no_dft'
assert report.get('source_manifest_sha256') == expected
assert report.get('dft_requested') is False
assert report.get('mpi_ranks') == 32
assert len(report.get('points', [])) == 26
assert all(row.get('status') == 'input_validated' and
           set(row.get('input_sha256', {})) == {'INPUT', 'KPT', 'STRU'}
           for row in report['points'])
assert p.get('calculation') == 'scf' and p.get('basis_type') == 'lcao'
assert p.get('dft_functional') == 'pbe' and p.get('ecutwfc') == 100.0
assert list(p.get('kpts', ())) == [4, 4, 4]
assert p.get('scf_thr') == 1e-8 and p.get('scf_nmax') == 150
assert p.get('mixing_beta') == 0.3
PY

module purge
source /public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/scripts/env.sh
unset I_MPI_PMI_LIBRARY I_MPI_HYDRA_BOOTSTRAP_EXEC_EXTRA_ARGS I_MPI_HYDRA_BOOTSTRAP
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "${code}"
"${python}" examples/run_bto_q1q2_frozen_abacus.py \
  --points-dir "${points}" --workdir "${workdir}" \
  --pseudo-dir "${pseudo_dir}" --basis-dir "${basis_dir}" \
  --abacus-bin "${abacus}" --mpi-ranks 32 --run-dft

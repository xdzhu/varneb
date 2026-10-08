#!/bin/bash
# Two reviewed historical peaks, identical INPUT/KPT/STRU/pseudo/orb files.
#SBATCH -J hfo2-channel-replica
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --array=0-1%2
#SBATCH --time=00:30:00

set -euo pipefail
[[ "${RUN_DFT:-0}" == 1 ]] || { echo "set RUN_DFT=1" >&2; exit 2; }
audit_root=${AUDIT_ROOT:?set new isolated audit root}
[[ "${SLURM_JOB_NUM_NODES:-0}" == 1 && "${SLURM_NTASKS:-0}" == 1 \
   && "${SLURM_CPUS_PER_TASK:-0}" == 32 ]] || { echo "requires one 32-CPU task" >&2; exit 2; }

source_root=/public/home/iai806/abacus/agent-runs/20260912-vcneb-hf/validation/hfo2_t_to_po
case "${SLURM_ARRAY_TASK_ID}" in
  0) source_dir=${source_root}/vcneb_n7_fire_distributed_resume27678004_v4/02; label=ordinary_peak02 ;;
  1) source_dir=${source_root}/vcneb_n7_mode_guided_unconstrained_continue300_maxstep003/02; label=guided_peak02 ;;
  *) echo "unreviewed array index" >&2; exit 2 ;;
esac
case_dir=${audit_root}/replicas/${label}
[[ ! -e "${case_dir}" ]] || { echo "refusing existing replica directory" >&2; exit 2; }
files=(INPUT KPT STRU Hf.upf O.upf Hf_gga_10au_100Ry_4s2p2d1f.orb O_gga_10au_100Ry_2s2p1d.orb)
for name in "${files[@]}"; do [[ -s "${source_dir}/${name}" ]] || exit 2; done
[[ $(sha256sum "${source_dir}/INPUT" | cut -d' ' -f1) == dc6684ffa709bf3d953c647272589b05d844c1293477bd144226d495207cb9df ]] || exit 2
[[ $(sha256sum "${source_dir}/KPT" | cut -d' ' -f1) == 92b917107d9df11da28a465cc4900f3574956a057ddf9f7ad125e75c956ec508 ]] || exit 2

module purge
source /public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/scripts/env.sh
unset I_MPI_PMI_LIBRARY I_MPI_HYDRA_BOOTSTRAP_EXEC_EXTRA_ARGS I_MPI_HYDRA_BOOTSTRAP
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export PYTHONPATH="${audit_root}/source${PYTHONPATH:+:${PYTHONPATH}}"
abacus=/public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/app/bin/abacus
python_bin=/public/home/iai806/.conda/envs/icu/bin/python
mkdir -p "${case_dir}"
for name in "${files[@]}"; do cp -p "${source_dir}/${name}" "${case_dir}/${name}"; done
cd "${case_dir}"
sha256sum "${abacus}" "${audit_root}/source/scripts/audit_hfo2_static_replica.py" > executable_sha256.txt
date --iso-8601=seconds > started_at.txt
mpirun -np 32 /bin/bash -c 'printf "rank=%s host=%s " "${PMI_RANK:-unset}" "$(hostname)"; grep Cpus_allowed_list /proc/self/status' > mpi_affinity.txt
/usr/bin/time -p mpirun -np 32 "${abacus}" > abacus.out 2> abacus.err
date --iso-8601=seconds > finished_at.txt
"${python_bin}" "${audit_root}/source/scripts/audit_hfo2_static_replica.py" \
  --source "${source_dir}" --replica "${case_dir}" --output "${case_dir}/replica_audit.json"

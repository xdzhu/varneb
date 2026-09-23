#!/bin/bash
# Isolated GaN ABACUS launcher/scaling canary. Never writes to the source image.
# Compare four and eight genuine MPI ranks on the same static structure.
#SBATCH -J varneb-aba-mpi-check
#SBATCH -p hfacnormal01
#SBATCH -N 1
# Request one Slurm task owning eight CPUs: Intel MPI's Hydra starts one
# bootstrap proxy with srun, so a one-CPU proxy would pin all ranks together.
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:45:00
#SBATCH -o /public/home/iai806/abacus/agent-runs/20260921-varneb-material/slurm-%x-%j.out
#SBATCH -e /public/home/iai806/abacus/agent-runs/20260921-varneb-material/slurm-%x-%j.err

set -euo pipefail
[[ "${RUN_DFT:-0}" == 1 ]] || { echo "set RUN_DFT=1" >&2; exit 2; }

source_dir=${SOURCE_IMAGE_DIR:?set SOURCE_IMAGE_DIR to a reviewed GaN image}
workdir=${WORKDIR:?set WORKDIR to a fresh canary directory}
abacus=${ABACUS_BIN:-/public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/app/bin/abacus}
env_script=${ABACUS_ENV:-/public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/scripts/env.sh}

[[ "${SLURM_JOB_NUM_NODES:-0}" == 1 && "${SLURM_NTASKS:-0}" == 1 \
   && "${SLURM_CPUS_PER_TASK:-0}" -ge 8 ]] || {
  echo "canary requires one Slurm task owning at least eight CPUs" >&2
  exit 2
}
[[ ! -e "${workdir}" ]] || { echo "refusing existing WORKDIR: ${workdir}" >&2; exit 2; }
for name in INPUT KPT STRU Ga.upf N.upf \
  Ga_gga_10au_100Ry_2s2p2d1f.orb N_gga_10au_100Ry_2s2p1d.orb; do
  [[ -s "${source_dir}/${name}" ]] || { echo "missing ${source_dir}/${name}" >&2; exit 2; }
done

module purge
source "${env_script}"
unset I_MPI_PMI_LIBRARY I_MPI_HYDRA_BOOTSTRAP_EXEC_EXTRA_ARGS I_MPI_HYDRA_BOOTSTRAP
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
mkdir -p "${workdir}"
printf 'source_image_dir=%s\nabacus=%s\njob=%s\n' \
  "${source_dir}" "${abacus}" "${SLURM_JOB_ID}" > "${workdir}/provenance.txt"
sha256sum "${source_dir}/INPUT" "${source_dir}/KPT" "${source_dir}/STRU" \
  >> "${workdir}/provenance.txt"
mpirun -np 4 /bin/bash -c 'printf "rank=%s host=%s " "${PMI_RANK:-unset}" "$(hostname)"; grep Cpus_allowed_list /proc/self/status' \
  > "${workdir}/mpi_affinity.txt"

for np in 4 8; do
  case_dir="${workdir}/np${np}"
  mkdir "${case_dir}"
  cp "${source_dir}/INPUT" "${source_dir}/KPT" "${source_dir}/STRU" \
    "${source_dir}/Ga.upf" "${source_dir}/N.upf" \
    "${source_dir}/Ga_gga_10au_100Ry_2s2p2d1f.orb" \
    "${source_dir}/N_gga_10au_100Ry_2s2p1d.orb" "${case_dir}/"
  (
    cd "${case_dir}"
    printf 'starting np=%s at %s\n' "${np}" "$(date --iso-8601=seconds)"
    # Match the installation's own case/abacus.slurm: mpirun -np N abacus.
    /usr/bin/time -p mpirun -np "${np}" "${abacus}" \
      > abacus.out 2> abacus.err
    if grep -q 'PMI server not found' abacus.out abacus.err; then
      echo "MPI launcher fell back to singleton mode for np=${np}" >&2
      exit 1
    fi
    printf 'finished np=%s at %s\n' "${np}" "$(date --iso-8601=seconds)"
  )
done

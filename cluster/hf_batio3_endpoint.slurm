#!/bin/bash
#SBATCH -J vcneb-bto-end
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH -o /public/home/iai806/abacus/agent-runs/20260912-vcneb-hf/slurm-%x-%j.out
#SBATCH -e /public/home/iai806/abacus/agent-runs/20260912-vcneb-hf/slurm-%x-%j.err

set -euo pipefail

endpoint=${1:?usage: sbatch hf_batio3_endpoint.slurm cubic|tetragonal}
case "${endpoint}" in
  cubic) structure="BaTiO3_cubic.vasp" ;;
  tetragonal) structure="BaTiO3_tetragonal.vasp" ;;
  *) echo "unknown endpoint: ${endpoint}" >&2; exit 2 ;;
esac

repo=/public/home/iai806/abacus/agent-runs/20260912-vcneb-hf
env_script=/public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/scripts/env.sh
python=/public/home/iai806/.conda/envs/icu/bin/python
abacus=/public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/app/bin/abacus
pseudo_dir=/public/home/iai806/abacus/PSEUDO/ABACUS-orbitals/Dojo-NC-FR/Pseudopotential
basis_dir=/public/home/iai806/abacus/PSEUDO/ABACUS-orbitals/Dojo-NC-FR/selected_Orbs

module purge
source "${env_script}"
export OMP_NUM_THREADS=1
export I_MPI_PMI_LIBRARY=/opt/gridview/slurm/lib/libpmi.so
cd "${repo}"

exec "${python}" examples/relax_abacus_endpoint.py \
  --structure "validation/batio3_cubic_to_tetragonal/source_structures/${structure}" \
  --workdir "validation/batio3_cubic_to_tetragonal/endpoint_relax_${endpoint}" \
  --command "srun --mpi=pmix_v3 -n ${SLURM_NTASKS} ${abacus}" \
  --pseudo-dir "${pseudo_dir}" \
  --basis-dir "${basis_dir}" \
  --pp Ba=Ba.upf --pp Ti=Ti.upf --pp O=O.upf \
  --basis Ba=Ba_gga_10au_100Ry_4s2p1d.orb \
  --basis Ti=Ti_gga_8au_100Ry_4s2p2d1f.orb \
  --basis O=O_gga_7au_100Ry_2s2p1d.orb \
  --ecutwfc 100 --scf-thr 1e-8 --scf-nmax 150 --mixing-beta 0.3 \
  --kpts 4 4 4 --fmax 0.01 --steps 200 --optimizer FIRE --maxstep 0.02

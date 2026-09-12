#!/bin/bash
#SBATCH -J vcneb-hfo2
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH -o /public/home/iai806/abacus/agent-runs/20260912-vcneb-hf/slurm-%x-%j.out
#SBATCH -e /public/home/iai806/abacus/agent-runs/20260912-vcneb-hf/slurm-%x-%j.err

set -eo pipefail

repo=/public/home/iai806/abacus/agent-runs/20260912-vcneb-hf
env_script=/public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/scripts/env.sh
python=/public/home/iai806/.conda/envs/icu/bin/python
abacus=/public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/app/bin/abacus
pseudo_dir=/public/home/iai806/abacus/PSEUDO/ABACUS-orbitals/Dojo-NC-FR/Pseudopotential
basis_dir=/public/home/iai806/abacus/PSEUDO/ABACUS-orbitals/Dojo-NC-FR/selected_Orbs

n_images=${N_IMAGES:-7}
steps=${STEPS:-300}
fmax=${FMAX:-0.05}
spring=${SPRING:-0.2}
optimizer=${OPTIMIZER:-FIRE}
maxstep=${MAXSTEP:-0.02}
climb_after=${CLIMB_AFTER:-100}
ecutwfc=${ECUTWFC:-60}
scf_thr=${SCF_THR:-1e-6}
scf_nmax=${SCF_NMAX:-100}
mixing_beta=${MIXING_BETA:-0.7}
kpt_x=${KPT_X:-1}
kpt_y=${KPT_Y:-1}
kpt_z=${KPT_Z:-1}
workdir=${WORKDIR:-validation/hfo2_t_to_po/vcneb_n${n_images}_${optimizer,,}_levelA}

module purge
source "${env_script}"
export OMP_NUM_THREADS=1
export I_MPI_PMI_LIBRARY=/opt/gridview/slurm/lib/libpmi.so
cd "${repo}"

climb_args=()
if [[ "${NO_CLIMB:-0}" != 0 ]]; then
  climb_args+=(--no-climb)
fi
resume_args=()
if [[ "${RESUME:-0}" != 0 ]]; then
  resume_args+=(--resume)
fi

exec "${python}" examples/run_vcneb_abacus.py \
  --initial validation/hfo2_t_to_po/relaxed_T/CONTCAR \
  --final validation/hfo2_t_to_po/relaxed_PO/CONTCAR \
  --workdir "${workdir}" --n-images "${n_images}" --fmax "${fmax}" \
  --steps "${steps}" --k "${spring}" --optimizer "${optimizer}" \
  --maxstep "${maxstep}" \
  --command "srun --mpi=pmix_v3 -n ${SLURM_NTASKS} ${abacus}" \
  --pseudo-dir "${pseudo_dir}" --basis-dir "${basis_dir}" \
  --pp Hf=Hf.upf --pp O=O.upf \
  --basis Hf=Hf_gga_7au_100Ry_4s2p2d1f.orb \
  --basis O=O_gga_7au_100Ry_2s2p1d.orb \
  --ecutwfc "${ecutwfc}" --scf-thr "${scf_thr}" --scf-nmax "${scf_nmax}" \
  --mixing-beta "${mixing_beta}" --kpts "${kpt_x}" "${kpt_y}" "${kpt_z}" \
  --mic --cell-interpolation log_strain --mapping auto --align-translation \
  --minimum-distance 1.6 --maximum-deformation 0.10 \
  --climb-after "${climb_after}" \
  "${climb_args[@]}" "${resume_args[@]}"

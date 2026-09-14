#!/bin/bash
#SBATCH -J vcneb-hfo2
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH -o /public/home/iai806/abacus/agent-runs/20260912-vcneb-hf/slurm-%x-%j.out
#SBATCH -e /public/home/iai806/abacus/agent-runs/20260912-vcneb-hf/slurm-%x-%j.err

set -eo pipefail

repo=/public/home/iai806/abacus/agent-runs/20260912-vcneb-hf
env_script=/public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/scripts/env.sh
python=/public/home/iai806/.conda/envs/icu/bin/python
abacus=/public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/app/bin/abacus
pseudo_dir=/public/home/iai806/abacus/PSEUDO/ABACUS-orbitals/Dojo-NC-FR/Pseudopotential
basis_dir=/public/home/iai806/abacus/PSEUDO/ABACUS-orbitals/Dojo-NC-FR/Orb-DZP-10au

n_images=${N_IMAGES:-7}
steps=${STEPS:-300}
fmax=${FMAX:-0.05}
spring=${SPRING:-0.2}
optimizer=${OPTIMIZER:-FIRE}
maxstep=${MAXSTEP:-0.02}
line_search_retries=${LINE_SEARCH_RETRIES:-0}
line_search_retry_factor=${LINE_SEARCH_RETRY_FACTOR:-0.5}
mode_file=${MODE_FILE:-}
mode_guided=${MODE_GUIDED:-0}
mode_amplitude=${MODE_AMPLITUDE:-0.25}
mode_envelope=${MODE_ENVELOPE:-sin}
constraint_mode=${CONSTRAINT_MODE:-none}
validate_only=${VALIDATE_ONLY:-0}
climb_after=${CLIMB_AFTER:-100}
ecutwfc=${ECUTWFC:-100}
scf_thr=${SCF_THR:-1e-8}
scf_nmax=${SCF_NMAX:-150}
mixing_beta=${MIXING_BETA:-0.3}
kpt_x=${KPT_X:-2}
kpt_y=${KPT_Y:-2}
kpt_z=${KPT_Z:-2}
minimum_distance=${MINIMUM_DISTANCE:-1.6}
maximum_deformation=${MAXIMUM_DEFORMATION:-0.25}
workdir=${WORKDIR:-validation/hfo2_t_to_po/vcneb_n${n_images}_${optimizer,,}_levelA}
initial_endpoint=${INITIAL_ENDPOINT:-validation/hfo2_t_to_po/relaxed_T/CONTCAR}
final_endpoint=${FINAL_ENDPOINT:-validation/hfo2_t_to_po/relaxed_PO/CONTCAR}

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

mode_args=()
if [[ -n "${mode_file}" ]]; then
  mode_args+=(--mode "${mode_file}")
fi
if [[ "${mode_guided}" != 0 ]]; then
  mode_args+=(--mode-guided --mode-amplitude "${mode_amplitude}" --mode-envelope "${mode_envelope}")
fi
if [[ "${constraint_mode}" != none ]]; then
  mode_args+=(--constraint-mode "${constraint_mode}")
fi
if [[ "${MODE_MASS_WEIGHTED_INPUT:-0}" != 0 ]]; then
  mode_args+=(--mode-mass-weighted-input)
fi
if [[ "${MODE_REMOVE_TRANSLATION:-0}" != 0 ]]; then
  mode_args+=(--mode-remove-translation)
fi
validate_args=()
if [[ "${validate_only}" != 0 ]]; then
  validate_args+=(--validate-only)
fi

exec "${python}" examples/run_vcneb_abacus.py \
  --initial "${initial_endpoint}" \
  --final "${final_endpoint}" \
  --workdir "${workdir}" --n-images "${n_images}" --fmax "${fmax}" \
  --steps "${steps}" --k "${spring}" --optimizer "${optimizer}" \
  --maxstep "${maxstep}" --line-search-retries "${line_search_retries}" \
  --line-search-retry-factor "${line_search_retry_factor}" \
  --command "srun --mpi=pmix_v3 -n ${SLURM_NTASKS} ${abacus}" \
  --pseudo-dir "${pseudo_dir}" --basis-dir "${basis_dir}" \
  --pp Hf=Hf.upf --pp O=O.upf \
  --basis Hf=Hf_gga_10au_100Ry_4s2p2d1f.orb \
  --basis O=O_gga_10au_100Ry_2s2p1d.orb \
  --ecutwfc "${ecutwfc}" --scf-thr "${scf_thr}" --scf-nmax "${scf_nmax}" \
  --mixing-beta "${mixing_beta}" --kpts "${kpt_x}" "${kpt_y}" "${kpt_z}" \
  --mic --cell-interpolation log_strain --mapping auto --align-translation \
  --minimum-distance "${minimum_distance}" --maximum-deformation "${maximum_deformation}" \
  --climb-after "${climb_after}" \
  "${mode_args[@]}" \
  "${validate_args[@]}" \
  "${climb_args[@]}" "${resume_args[@]}"

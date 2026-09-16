#!/bin/bash
# Gamma-point BTO phonon force constants at one stationary endpoint.
# Default RUN_DFT=0 only creates/audits the finite-displacement inputs.  This
# workflow is separate from NEB: every displacement is a static ABACUS SCF.
#SBATCH -J varneb-bto-gamma
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH -o /public/home/iai806/abacus/agent-runs/20260912-vcneb-hf/slurm-%x-%j.out
#SBATCH -e /public/home/iai806/abacus/agent-runs/20260912-vcneb-hf/slurm-%x-%j.err

set -euo pipefail

repo=${REPO:-/public/home/iai806/abacus/agent-runs/20260912-vcneb-hf}
env_script=/public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/scripts/env.sh
python=/public/home/iai806/.conda/envs/icu/bin/python
phonopy=/public/home/iai806/.conda/envs/icu/bin/phonopy
abacus=/public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/app/bin/abacus
pseudo_dir=/public/home/iai806/abacus/PSEUDO/ABACUS-orbitals/Dojo-NC-FR/Pseudopotential
basis_dir=/public/home/iai806/abacus/PSEUDO/ABACUS-orbitals/Dojo-NC-FR/Orb-DZP-10au
endpoint=${ENDPOINT:-cubic}
endpoint_tag=${ENDPOINT_TAG:-pbe100_dzp10au}
delta=${DELTA:-0.01}
workdir=${WORKDIR:-validation/batio3_cubic_to_tetragonal/gamma_phonon_${endpoint}_${endpoint_tag}}

case "${endpoint}" in
  cubic|tetragonal) ;;
  *) echo "ENDPOINT must be cubic or tetragonal" >&2; exit 2 ;;
esac
export ENDPOINT="${endpoint}"

source "${env_script}"
export OMP_NUM_THREADS=1
export I_MPI_PMI_LIBRARY=/opt/gridview/slurm/lib/libpmi.so
cd "${repo}"

source_stru="validation/batio3_cubic_to_tetragonal/endpoint_relax_${endpoint}_${endpoint_tag}/00/STRU"
if [[ ! -f "${source_stru}" ]]; then
  echo "missing audited endpoint STRU: ${source_stru}" >&2
  exit 2
fi

mkdir -p "${workdir}"
cp "${source_stru}" "${workdir}/STRU"
cd "${workdir}"

# A 1x1x1 supercell is deliberate for the first BTO Gamma-mode decomposition.
# It provides Gamma force constants, not a converged dispersion or QHA result.
"${phonopy}" --abacus -d --pm --dim "1 1 1" --amplitude "${delta}" -c STRU
mapfile -t displaced_strus < <(printf '%s\n' STRU-[0-9][0-9][0-9] | sort)
if (( ${#displaced_strus[@]} == 0 )); then
  echo "phonopy generated no displaced STRU files" >&2
  exit 2
fi

for stru in "${displaced_strus[@]}"; do
  tag=${stru#STRU-}
  dir="disp-${tag}"
  mkdir -p "${dir}"
  cp "${stru}" "${dir}/STRU"
  cat > "${dir}/INPUT" <<EOF
INPUT_PARAMETERS
calculation scf
basis_type lcao
ks_solver genelpa
dft_functional pbe
ecutwfc 100
scf_thr 1e-8
scf_nmax 150
mixing_beta 0.3
cal_force 1
cal_stress 1
out_stru 1
symmetry 0
pseudo_dir ${pseudo_dir}
orbital_dir ${basis_dir}
EOF
  cat > "${dir}/KPT" <<'EOF'
K_POINTS
0
Gamma
4 4 4 0 0 0
EOF
done

"${python}" - <<'PY'
import json
from pathlib import Path
root = Path('.')
payload = {
    'status': 'prepared_no_dft',
    'endpoint': __import__('os').environ['ENDPOINT'],
    'supercell': [1, 1, 1],
    'finite_difference_amplitude_A': float(__import__('os').environ.get('DELTA', '0.01')),
    'displacement_directories': [str(path) for path in sorted(root.glob('disp-[0-9][0-9][0-9]'))],
    'required_outputs_before_mode_claim': ['FORCE_SETS', 'phonopy.yaml', 'force-constant provenance'],
}
Path('phonon_preflight.json').write_text(json.dumps(payload, indent=2) + '\n', encoding='utf-8')
PY

if [[ "${RUN_DFT:-0}" != 1 ]]; then
  echo "[OK] BTO Gamma phonon inputs prepared; RUN_DFT=0 so no ABACUS SCF was launched"
  exit 0
fi

for dir in disp-[0-9][0-9][0-9]; do
  (
    cd "${dir}"
    srun --exclusive --nodes=1 --ntasks=32 --ntasks-per-node=32 "${abacus}" > abacus.stdout 2> abacus.stderr
  )
done

"${phonopy}" --abacus -f disp-[0-9][0-9][0-9]/OUT.ABACUS/running_scf.log
test -s FORCE_SETS
echo "[DONE] Gamma force constants assembled in ${PWD}/FORCE_SETS"

#!/usr/bin/env bash
# Add only 72 new interleaved q nodes to the previously audited 18x5 cut.
# Submit explicitly: sbatch --export=ALL,RUN_DFT=1 hf_gan_600eV_atomic_tube_q9.slurm
#SBATCH -J varneb-gan600-q9
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=32
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --array=0-71%8
#SBATCH --time=00:45:00
#SBATCH -o /public/home/iai806/abacus/agent-runs/20260929-varneb-gan600-atomic-q9/slurm-%A_%a.out
#SBATCH -e /public/home/iai806/abacus/agent-runs/20260929-varneb-gan600-atomic-q9/slurm-%A_%a.err

set -eo pipefail
source /public/home/iai806/Software/VASP/env.sh 6.3.2
set -u
[[ "${RUN_DFT:-0}" == 1 && "${SLURM_JOB_PARTITION:-}" == hfacnormal01 \
   && "${SLURM_JOB_NUM_NODES:-0}" == 1 && "${SLURM_NTASKS:-0}" == 32 \
   && "${SLURM_CPUS_PER_TASK:-0}" == 1 ]] || exit 2
root=/public/home/iai806/abacus/agent-runs/20260929-varneb-gan600-atomic-q9
work=${root}/gan_600eV_atomic_q9_20260929
python=/public/home/iai806/.conda/envs/icu/bin/python
[[ -s "${work}/manifest.json" && "${SLURM_ARRAY_TASK_ID:-}" =~ ^[0-9]+$ ]] || exit 3
printf '%s  %s\n' 'e87b630daf996dd6cb67e5abbbae9bc0c4eba51dae1d6ee7772fdfe873f86a43' \
  "${work}/manifest.json" | sha256sum -c -

case_name=$("${python}" - "${work}/manifest.json" "${SLURM_ARRAY_TASK_ID}" <<'PY'
import json
import sys
from pathlib import Path

report = json.loads(Path(sys.argv[1]).read_text())
index = int(sys.argv[2])
assert report['status'] == '72_interleaved_inputs_finalized_no_DFT'
assert report['purpose'] == 'GaN_45p7_600eV_central_18x9_frozen_atomic_transverse_cut'
assert report['pressure_GPa'] == 45.7
assert report['n_old_measured_points'] == 90 and report['n_new_statics'] == 72
assert report['prospective_error_gate_meV_per_GaN'] == 1.0
assert report['input_contract_sha256']['INCAR'] == '83ba34d4b14ff4ea641bd74dec26e462ea1cc24b575db79a07138ccdfd552772'
assert report['input_contract_sha256']['KPOINTS'] == 'b5215f3e7608c27f49d45e8129d3edca2cfaa3ba88b8c389d4b52604715712ee'
assert report['input_contract_sha256']['POTCAR'] == 'f94781ce6cf9b9454c093353320c5e80a2a0aceea6fb8353b33b01ac37b95168'
assert len(report['cases']) == 72 and 0 <= index < 72
case = report['cases'][index]
assert 5 <= case['image_index'] <= 22
assert case['q_atom_A'] in (-0.0375, -0.0125, 0.0125, 0.0375)
print(case['name'])
PY
)
case_dir=${work}/cases/${case_name}
[[ -s "${case_dir}/sha256.inputs.json" && ! -e "${case_dir}/OUTCAR" \
   && ! -e "${case_dir}/vasp.stdout" ]] || exit 3
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export I_MPI_PMI_LIBRARY=/opt/gridview/slurm/lib/libpmi2.so
cd "${case_dir}"
"${python}" - <<'PY'
import hashlib
import json
from pathlib import Path

expected = json.loads(Path('sha256.inputs.json').read_text())
assert set(expected) == {'POSCAR', 'INCAR', 'KPOINTS', 'POTCAR'}
assert all(hashlib.sha256(Path(name).read_bytes()).hexdigest() == digest
           for name, digest in expected.items())
assert expected['INCAR'] == '83ba34d4b14ff4ea641bd74dec26e462ea1cc24b575db79a07138ccdfd552772'
PY
printf 'job=%s task=%s case=%s node=%s\n' \
  "${SLURM_JOB_ID}" "${SLURM_ARRAY_TASK_ID}" "${case_name}" "${SLURM_JOB_NODELIST}"
srun --exclusive --exact --nodes=1 --ntasks=32 --ntasks-per-node=32 \
  --cpus-per-task=1 --mpi=pmi2 "${VASP_ROOT}/bin/vasp_std" > vasp.stdout 2>&1
grep -q 'General timing and accounting informations' OUTCAR

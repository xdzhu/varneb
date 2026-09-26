#!/bin/bash
# GaN four-atom local 1000 eV diagnostic static; no phonon supercell.
#SBATCH -J varneb-gan-ts-followup
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00

set -eo pipefail
source /public/home/iai806/Software/VASP/env.sh 6.3.2
set -u
[[ "${RUN_DFT:-0}" == 1 && "${SLURM_NTASKS:-0}" == 32 ]] || exit 2
: "${WORK_ROOT:?set isolated immutable staged WORK_ROOT}"
case_dir="${WORK_ROOT}/cases/newton_followup_full"
[[ -s "${WORK_ROOT}/manifest.json" && -s "${case_dir}/sha256.inputs.json" && ! -e "${case_dir}/OUTCAR" ]] || exit 3
python=/public/home/iai806/.conda/envs/icu/bin/python
"${python}" -c 'import hashlib,json,pathlib,sys; root=pathlib.Path(sys.argv[1]); m=json.loads((root/"manifest.json").read_text()); c=root/"cases/newton_followup_full"; h=json.loads((c/"sha256.inputs.json").read_text()); assert m["status"]=="inputs_finalized_no_DFT" and len(m["cases"])==1 and m["cases"][0]["input_sha256"]==h; assert all(hashlib.sha256((c/k).read_bytes()).hexdigest()==v for k,v in h.items())' "${WORK_ROOT}"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export I_MPI_PMI_LIBRARY=/opt/gridview/slurm/lib/libpmi2.so
cd "${case_dir}"
printf 'job=%s node=%s case=newton_followup_full\n' "${SLURM_JOB_ID}" "${SLURM_JOB_NODELIST}"
srun --exclusive --exact --nodes=1 --ntasks=32 --ntasks-per-node=32 --cpus-per-task=1 \
  --mpi=pmi2 "${VASP_ROOT}/bin/vasp_std" > vasp.stdout 2>&1
grep -q 'General timing and accounting informations' OUTCAR

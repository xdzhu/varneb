#!/bin/bash
# Fresh SYMPREC=1e-5 production namespace. Eighteen interior images, six
# concurrent 32-MPI workers, three batches per complete image evaluation.
#SBATCH -J hfo2-1e5-vcneb
#SBATCH -p hfacnormal01
#SBATCH -N 2
#SBATCH --ntasks=194
#SBATCH --ntasks-per-node=97
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
set -eo pipefail
source /public/home/iai806/Software/VASP/env.sh 6.3.2
set -u

[[ "${RUN_DFT:-0}" == 1 && "${SLURM_NTASKS}" == 194 && "${SLURM_JOB_NUM_NODES}" == 2 ]] || exit 2
case "${CASE_NAME:-}" in t_to_po|po_to_m) ;; *) exit 2 ;; esac

repo=/public/home/iai806/abacus/agent-runs/20260919-prl2023-hfo2-fig2a
python=/public/home/iai806/.conda/envs/icu/bin/python
probe=/public/home/iai806/abacus/agent-runs/20260918-varneb-gan-suite/diagnostics/hex_image10_basis/vasp_lattice_probe
case_root="${repo}/cases/${CASE_NAME}"
initial="${case_root}/endpoint_relax/initial"
final="${case_root}/endpoint_relax/final"
path="${case_root}/initial_corrected.traj"
workdir="${case_root}/vcneb_n20_fmax010_no_ci_symprec1e5"
export OMP_NUM_THREADS=1 PYTHONUNBUFFERED=1
export I_MPI_PMI_LIBRARY=/opt/gridview/slurm/lib/libpmi2.so
export VASP_COMMAND="srun --exclusive --exact --nodes=1 --ntasks=32 --ntasks-per-node=32 --cpus-per-task=1 --distribution=block:block --mpi=pmi2 ${VASP_ROOT}/bin/vasp_std"
export IMAGE_CACHE_NAMESPACE="hf-vasp6.3.2-hfo2-standard-symprec1e5-$(sha256sum "${VASP_ROOT}/bin/vasp_std" | cut -d' ' -f1)"
cd "${repo}"

"${python}" -c 'import json,sys; x=json.load(open(sys.argv[1])); assert x["all_passed"] and x["symprec_A"] == 1e-5 and not x["old_symprec1e4_cache_reused"]' \
  "${case_root}/symprec1e5_gate.json"
"${python}" examples/run_vcneb_vasp.py \
  --initial "${initial}" --final "${final}" --initial-trajectory "${path}" \
  --workdir "${workdir}" --n-images 20 --fmax 0.10 --steps 300 --k 0.10 \
  --pressure-gpa 0 --optimizer FIRE --image-workers 6 --image-retries 0 --no-climb \
  --initial-static-summary "${case_root}/symprec1e5_static_initial/vasp_static_summary.json" \
  --final-static-summary "${case_root}/symprec1e5_static_final/vasp_static_summary.json" \
  --image-cache-namespace "${IMAGE_CACHE_NAMESPACE}" \
  --vasp-isym -1 --vasp-symprec 1e-5 --native-lattice-probe "${probe}" \
  --candidate-step-retries 8 --candidate-step-retry-factor 0.5 \
  --minimum-distance 1.40 --maximum-deformation 0.50 \
  --ncores 32 --vasp-bin "${VASP_ROOT}/bin/vasp_std"
"${python}" scripts/audit_vcneb_result.py "${workdir}" \
  --max-min-distance 1.40 --max-deformation 0.50 --fmax-target 0.10

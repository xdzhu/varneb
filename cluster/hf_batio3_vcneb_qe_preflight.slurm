#!/bin/bash
# One-rank, no-DFT BTO QE VCNEB input/provenance gate.
# It checks endpoints, the seven-total-image geometry, static QE semantics,
# and reviewed Ba/Ti/O UPF identities without launching pw.x.
#SBATCH -J varneb-bto-qe-preflight
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -o /public/home/iai806/abacus/agent-runs/20260912-vcneb-hf/slurm-%x-%j.out
#SBATCH -e /public/home/iai806/abacus/agent-runs/20260912-vcneb-hf/slurm-%x-%j.err

set -euo pipefail

repo=${REPO:-/public/home/iai806/abacus/agent-runs/20260912-vcneb-hf}
python=${PYTHON:-/public/home/iai806/.conda/envs/icu/bin/python}
qe_bin=${QE_BIN:-pw.x}
pseudo_dir=${ESPRESSO_PSEUDO:-}
endpoint_tag=${ENDPOINT_TAG:-pbe100_dzp10au}
workdir=${WORKDIR:-validation/batio3_cubic_to_tetragonal/vcneb_tetragonal_to_cubic_n7_qe_preflight}
ecutwfc=${ECUTWFC:-100}
ecutrho=${ECUTRHO:-400}
scf_thr=${SCF_THR:-1e-8}
pp_ba=${PP_BA:-Ba.upf}
pp_ti=${PP_TI:-Ti.upf}
pp_o=${PP_O:-O.upf}
pp_manifest=${QE_PP_MANIFEST:-}

[[ -n "${pseudo_dir}" ]] || { echo "Set ESPRESSO_PSEUDO to the reviewed UPF directory" >&2; exit 2; }
[[ -n "${pp_manifest}" ]] || { echo "Set QE_PP_MANIFEST to an approved MD5-pinned UPF manifest" >&2; exit 2; }
cd "${repo}"

exec "${python}" examples/run_vcneb_qe.py \
  --initial "validation/batio3_cubic_to_tetragonal/endpoint_relax_tetragonal_${endpoint_tag}/CONTCAR" \
  --final "validation/batio3_cubic_to_tetragonal/endpoint_relax_cubic_${endpoint_tag}/CONTCAR" \
  --workdir "${workdir}" --n-images 7 --fmax 0.10 --steps 0 --k 0.2 \
  --image-workers 0 --command "${qe_bin}" --pseudo-dir "${pseudo_dir}" \
  --pp-manifest "${pp_manifest}" \
  --ecutwfc "${ecutwfc}" --ecutrho "${ecutrho}" --scf-thr "${scf_thr}" \
  --kpts 4 4 4 --no-climb --mic --cell-interpolation log_strain --mapping auto \
  --align-translation --minimum-distance 1.6 --maximum-deformation 0.10 --validate-only

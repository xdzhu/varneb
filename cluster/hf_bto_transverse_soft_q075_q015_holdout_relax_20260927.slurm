#!/usr/bin/env bash
# Pre-registered BTO (Qz,Qx)=(0.75,0.15) three-start conditional holdout.
# ABACUS PBE/100 Ry/Ba-Ti-O 10 au DZP; no physical-parameter retuning.
#SBATCH -J varneb-bto-q075-relax
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=05:00:00
#SBATCH -o /public/home/iai806/abacus/agent-runs/20260927-varneb-bto-soft-gridpilot-qx030/slurm-q075-relax-%j.out
#SBATCH -e /public/home/iai806/abacus/agent-runs/20260927-varneb-bto-soft-gridpilot-qx030/slurm-q075-relax-%j.err

set -euo pipefail
[[ "${RUN_DFT:-0}" == 1 ]] || { echo 'set RUN_DFT=1 explicitly' >&2; exit 2; }
[[ "${SLURM_JOB_PARTITION:-}" == hfacnormal01 && "${SLURM_JOB_NUM_NODES:-0}" == 1 \
   && "${SLURM_NTASKS:-0}" == 1 && "${SLURM_CPUS_PER_TASK:-0}" -ge 32 ]] || exit 2

root=/public/home/iai806/abacus/agent-runs/20260927-varneb-bto-soft-gridpilot-qx030
code=${root}/code
inputs=${root}/inputs
workdir=${root}/run-q075_q015
preflight=${inputs}/bto_transverse_soft_conditional_q075_q015_preflight_2026-09-27.json
canary=${workdir}/canary_starts_result.json
audit=${root}/audit-q075_q015-canary-27787471.json
result=conditional_q075_q015_result.json
python=/public/home/iai806/.conda/envs/icu/bin/python
abacus=/public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/app/bin/abacus
pseudo_dir=/public/home/iai806/abacus/PSEUDO/ABACUS-orbitals/Dojo-NC-FR/Pseudopotential
basis_dir=/public/home/iai806/abacus/PSEUDO/ABACUS-orbitals/Dojo-NC-FR/Orb-DZP-10au
env_script=/public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/scripts/env.sh

[[ -s "${workdir}/contract.json" && -s "${preflight}" && -s "${canary}" \
   && -s "${audit}" && ! -e "${workdir}/${result}" ]] || exit 3
echo '1eb33df77355b3c2282333994402110b705318004b5d0305c2fb745967f73ee4  '"${preflight}" | sha256sum -c -
echo '105df820654ca3bdd505fd90f056d999be4c2bef9ee6f981b47a530c263892f1  '"${canary}" | sha256sum -c -
echo '1708ed8092354c851ae19aa91e4bc99cb6546cea4c32f8963633d7f5c077652d  '"${audit}" | sha256sum -c -
echo 'de0c1c3aeac376b9e4f4ea1f9ff742e85e70f9ba08823824c46c358d547a17bf  '"${code}/examples/run_bto_q1q2_conditional_abacus.py" | sha256sum -c -
echo '91b70d48a488b5a26f1e563bbce1ba0e75581a01d9145158d57108b8dac4bbfc  '"${code}/vcneb/mode_surface.py" | sha256sum -c -

"${python}" - "${preflight}" "${canary}" "${audit}" <<'PY'
import hashlib, json, pathlib, sys
p, c, a = map(pathlib.Path, sys.argv[1:])
sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
preflight, canary, audit = (json.loads(path.read_text()) for path in (p, c, a))
if not (preflight['q_parallel_q_transverse_sqrt_amu_A'] == [0.75, 0.15]
        and preflight['phonon_supercell'] == [1, 1, 1]
        and preflight['electronic_kpoints'] == [4, 4, 4]
        and canary['status'] == 'three_start_statics_complete_optimizer_not_run'
        and canary['preflight_sha256'] == sha(p)
        and str(canary['slurm_job_id']) == '27787471'
        and audit['status'] == 'all_three_raw_DFT_points_validated'
        and audit['source_sha256']['preflight'] == sha(p)
        and audit['source_sha256']['canary_result'] == sha(c)
        and audit['q_parallel_q_transverse_sqrt_amu_A'] == [0.75, 0.15]
        and len(audit['points']) == 3):
    raise SystemExit('holdout canary/audit source contract failed')
PY

module purge
source "${env_script}"
unset I_MPI_PMI_LIBRARY I_MPI_HYDRA_BOOTSTRAP_EXEC_EXTRA_ARGS I_MPI_HYDRA_BOOTSTRAP
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
cd "${code}"
"${python}" examples/run_bto_q1q2_conditional_abacus.py \
  --preflight "${preflight}" \
  --report "${inputs}/bto_t_to_c_q1q2_strain_gradient_preflight_2026-09-25.json" \
  --reference "${inputs}/cubic_CONTCAR" \
  --force-constants "${inputs}/bto_cubic_gamma_force_constants.npz" \
  --phonopy-eigenpairs "${inputs}/bto_cubic_phonopy_gamma_eigenpairs.npz" \
  --gamma-provenance "${inputs}/bto_cubic_gamma_phonon_provenance.json" \
  --force-sets "${inputs}/bto_cubic_gamma_FORCE_SETS" \
  --eigenpairs-provenance "${inputs}/bto_cubic_phonopy_gamma_eigenpairs_provenance.json" \
  --grid-result-manifest "${inputs}/result_manifest_hf_27777454.json" \
  --workdir "${workdir}" --pseudo-dir "${pseudo_dir}" --basis-dir "${basis_dir}" \
  --abacus-bin "${abacus}" --mpi-ranks 32 \
  --gradient-tolerance 0.003 --orthogonal-amplitude-bound 3.0 \
  --max-iterations 32 --max-new-evaluations 120 \
  --max-absolute-stress-kbar 1.0 \
  --result-name "${result}" --resume --run-dft

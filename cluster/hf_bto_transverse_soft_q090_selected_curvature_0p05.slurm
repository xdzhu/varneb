#!/bin/bash
# Signed 16D fixed-Q curvature probes at two audited BTO conditional branches.
# Phonon basis is the existing 5-atom 1x1x1 Gamma calculation. ABACUS settings
# remain PBE / ecutwfc=100 Ry / 10 au DZP / 4x4x4 electronic k points.
# This is a one-step curvature screen, not a conditional-PES certificate.
#SBATCH -J varneb-bto-q090-hess
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --array=0-1%2
#SBATCH --time=02:00:00
#SBATCH -o /public/home/iai806/abacus/agent-runs/20260927-varneb-bto-soft-gridpilot-qx030/slurm-hess-%A_%a.out
#SBATCH -e /public/home/iai806/abacus/agent-runs/20260927-varneb-bto-soft-gridpilot-qx030/slurm-hess-%A_%a.err

set -euo pipefail
[[ "${RUN_DFT:-0}" == 1 ]] || { echo 'set RUN_DFT=1 explicitly' >&2; exit 2; }
[[ "${SLURM_JOB_PARTITION:-}" == hfacnormal01 && "${SLURM_JOB_NUM_NODES:-0}" == 1 \
   && "${SLURM_NTASKS:-0}" == 1 && "${SLURM_CPUS_PER_TASK:-0}" -ge 32 ]] || exit 2

root=/public/home/iai806/abacus/agent-runs/20260927-varneb-bto-soft-gridpilot-qx030
code=${root}/code
inputs=${root}/inputs
python=/public/home/iai806/.conda/envs/icu/bin/python
abacus=/public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/app/bin/abacus
pseudo_dir=/public/home/iai806/abacus/PSEUDO/ABACUS-orbitals/Dojo-NC-FR/Pseudopotential
basis_dir=/public/home/iai806/abacus/PSEUDO/ABACUS-orbitals/Dojo-NC-FR/Orb-DZP-10au
env_script=/public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/scripts/env.sh

case "${SLURM_ARRAY_TASK_ID:-}" in
  0)
    point=q090_q000
    parent_hash=f6efab62b0dcf9b6e90dc0fb553a13610c600400af10a37bf4cc8357c560aa07
    raw_hash=1fabc39e2db2bbdde37f6ac3a3dfb27df8bf988c60fef7a06d32eed026966bb1
    replay_hash=b1bab7a048dd0b3327eb24e14b539286ba9594ace58526d7e947bcf156d8d37a
    preflight_hash=2df4bdb5d4588078c3cdc89b672aa41582c4aff2ba90b643d55c409e96282f54
    ;;
  1)
    point=q090_q030
    parent_hash=7e92c953c81676ba3999e544fb05aad637315d72342034316cd53ad58a21ef87
    raw_hash=98640ba39926aa5ce6521eb18176aee0efeeedfaa2a0200609321fe978bd1fd6
    replay_hash=c152c6db670ef4029ac45e3f61759a95e862e05df4c09d0dd846c39d35314299
    preflight_hash=255806fa27b6af37d9902d366d2fbb959cfc4ffc4c2fa2265993076701f4d06e
    ;;
  *) echo 'unexpected array index' >&2; exit 2 ;;
esac

workdir=${root}/run-${point}
parent=${workdir}/conditional_${point}_result.json
raw_audit=${root}/audit-${point}-production-raw-final.json
branch_replay=${root}/replay-${point}-final.json
geometry_preflight=${inputs}/bto_transverse_soft_${point}_selected_curvature_0p05_preflight_2026-09-27.json
result=conditional_${point}_selected_curvature_0p05.json
[[ -s "${workdir}/contract.json" && -s "${parent}" && -s "${raw_audit}" \
   && -s "${branch_replay}" && -s "${geometry_preflight}" \
   && ! -e "${workdir}/${result}" ]] || exit 3
echo "${parent_hash}  ${parent}" | sha256sum -c -
echo "${raw_hash}  ${raw_audit}" | sha256sum -c -
echo "${replay_hash}  ${branch_replay}" | sha256sum -c -
echo "${preflight_hash}  ${geometry_preflight}" | sha256sum -c -
echo 'de0c1c3aeac376b9e4f4ea1f9ff742e85e70f9ba08823824c46c358d547a17bf  '"${code}/examples/run_bto_q1q2_conditional_abacus.py" | sha256sum -c -
echo '91b70d48a488b5a26f1e563bbce1ba0e75581a01d9145158d57108b8dac4bbfc  '"${code}/vcneb/mode_surface.py" | sha256sum -c -

"${python}" - "${parent}" "${raw_audit}" "${branch_replay}" "${geometry_preflight}" <<'PY'
import hashlib, json, pathlib, sys
parent, audit, replay, geometry = map(pathlib.Path, sys.argv[1:])
sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
s, a, r, p = (json.loads(path.read_text()) for path in (parent, audit, replay, geometry))
if not (s['status'] == 'orthogonal_gradient_and_stress_converged_curvature_unchecked'
        and s['stress_target_passed']
        and a['status'] == 'verified_gradient_stationary_candidate_curvature_and_branches_unchecked'
        and a['source_sha256']['summary'] == sha(parent)
        and r['status'] == 'all_three_branch_outcomes_reproduced_from_raw_audited_cache'
        and r['source_sha256']['raw_audit'] == sha(audit)
        and r['selected_start'] == s['selected_start']
        and p['kind'] == 'bto_transverse_soft_curvature_probe_geometry_preflight_no_dft'
        and p['source_sha256']['summary'] == sha(parent)
        and p['source_sha256']['audit'] == sha(audit)
        and p['source_sha256']['branch_replay'] == sha(replay)
        and p['n_signed_probes'] == 32
        and p['curvature_step_sqrt_amu_A'] == 0.05):
    raise SystemExit('curvature source and geometry contract failed')
PY

module purge
source "${env_script}"
unset I_MPI_PMI_LIBRARY I_MPI_HYDRA_BOOTSTRAP_EXEC_EXTRA_ARGS I_MPI_HYDRA_BOOTSTRAP
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
cd "${code}"
"${python}" examples/run_bto_q1q2_conditional_abacus.py \
  --preflight "${inputs}/bto_transverse_soft_conditional_${point}_preflight_2026-09-27.json" \
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
  --max-iterations 1 --max-new-evaluations 32 \
  --max-absolute-stress-kbar 1.0 --resume --run-dft \
  --start-from-result "${parent}" --curvature-only --curvature-step 0.05 \
  --result-name "${result}"

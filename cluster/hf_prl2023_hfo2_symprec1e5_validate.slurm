#!/bin/bash
# Full-chain native lattice gate, one failed-image full SCF canary, and fresh
# endpoint statics under the uniform SYMPREC=1e-5 contract.
#SBATCH -J hfo2-1e5-check
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --time=08:00:00
set -eo pipefail
source /public/home/iai806/Software/VASP/env.sh 6.3.2
set -u

[[ "${RUN_DFT:-0}" == 1 && "${SLURM_NTASKS}" == 32 ]] || exit 2
case "${CASE_NAME:-}" in
  t_to_po) canary_index=16 ;;
  po_to_m) canary_index=18 ;;
  *) exit 2 ;;
esac

repo=/public/home/iai806/abacus/agent-runs/20260919-prl2023-hfo2-fig2a
python=/public/home/iai806/.conda/envs/icu/bin/python
probe=/public/home/iai806/abacus/agent-runs/20260918-varneb-gan-suite/diagnostics/hex_image10_basis/vasp_lattice_probe
case_root="${repo}/cases/${CASE_NAME}"
initial="${case_root}/endpoint_relax/initial"
final="${case_root}/endpoint_relax/final"
path="${case_root}/initial_corrected.traj"
gate="${case_root}/symprec1e5_gate.json"
export OMP_NUM_THREADS=1 PYTHONUNBUFFERED=1
export I_MPI_PMI_LIBRARY=/opt/gridview/slurm/lib/libpmi2.so
export VASP_COMMAND="srun --exclusive --exact --nodes=1 --ntasks=32 --ntasks-per-node=32 --cpus-per-task=1 --mpi=pmi2 ${VASP_ROOT}/bin/vasp_std"
cd "${repo}"

[[ -x "${probe}" && -s "${path}" ]] || exit 3
for endpoint in initial final; do
  "${python}" -c 'import json,sys; x=json.load(open(sys.argv[1])); assert x["converged"] and x["fmax_target_eV_per_A"] == 0.001' \
    "${case_root}/endpoint_relax/${endpoint}/endpoint_relax_summary.json"
done

# This calls the licensed VASP lattice classifier on manager and exact POSCAR
# round-trip cells for all twenty images without altering their geometry.
"${python}" examples/run_vcneb_vasp.py \
  --initial "${initial}" --final "${final}" --initial-trajectory "${path}" \
  --workdir "${case_root}/symprec1e5_preflight" --n-images 20 \
  --pressure-gpa 0 --no-climb --validate-only --vasp-isym -1 --vasp-symprec 1e-5 \
  --native-lattice-probe "${probe}" --minimum-distance 1.40 --maximum-deformation 0.50 \
  --ncores 32 --vasp-bin "${VASP_ROOT}/bin/vasp_std"

# Re-run the image that failed at SYMPREC=1e-4 as a complete static SCF.
"${python}" scripts/probe_vasp_input_contract.py \
  --source-dir "${initial}" --trajectory "${path}" --n-images 20 \
  --trajectory-index "${canary_index}" \
  --workdir "${case_root}/symprec1e5_canary_image${canary_index}" \
  --symprec 1e-5 --ncores 32 --vasp-bin "${VASP_ROOT}/bin/vasp_std" \
  --timeout 1800 --full-scf --run

# Endpoint energies cached into production must use the same symmetry contract.
for endpoint in initial final; do
  "${python}" examples/run_vcneb_vasp.py \
    --initial "${initial}" --final "${final}" --initial-trajectory "${path}" \
    --workdir "${case_root}/symprec1e5_static_${endpoint}" --n-images 20 \
    --pressure-gpa 0 --no-climb --static-only --static-endpoint "${endpoint}" \
    --vasp-isym -1 --vasp-symprec 1e-5 --native-lattice-probe "${probe}" \
    --minimum-distance 1.40 --maximum-deformation 0.50 \
    --ncores 32 --vasp-bin "${VASP_ROOT}/bin/vasp_std"
done

"${python}" - "${CASE_NAME}" "${canary_index}" "${probe}" "${case_root}" "${gate}" <<'PY'
import hashlib, json, sys
from pathlib import Path
case, canary_index, probe, root, output = sys.argv[1:]
root, output, probe = Path(root), Path(output), Path(probe)
canary = json.load(open(root / f"symprec1e5_canary_image{canary_index}" / "probe_summary.json"))
endpoints = {
    side: json.load(open(root / f"symprec1e5_static_{side}" / "vasp_static_summary.json"))
    for side in ("initial", "final")
}
passed = (
    canary.get("all_passed") is True
    and all(record.get("status") == "converged_static_passed" for record in canary["records"])
    and all(item.get("status") == "completed" for item in endpoints.values())
    and all(item["calculator_parameters"]["symprec"] == 1e-5 for item in endpoints.values())
)
payload = {
    "status": "passed" if passed else "failed",
    "all_passed": passed,
    "case": case,
    "isym": -1,
    "symprec_A": 1e-5,
    "n_images_total": 20,
    "native_lattice_probe": {"path": str(probe), "sha256": hashlib.sha256(probe.read_bytes()).hexdigest()},
    "canary_image_index": int(canary_index),
    "canary_summary": str(root / f"symprec1e5_canary_image{canary_index}" / "probe_summary.json"),
    "endpoint_static_summaries": {
        side: str(root / f"symprec1e5_static_{side}" / "vasp_static_summary.json")
        for side in endpoints
    },
    "old_symprec1e4_cache_reused": False,
}
output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
if not passed:
    raise SystemExit(1)
PY

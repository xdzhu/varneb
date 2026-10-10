set -euo pipefail
source_root=/public/home/iai806/abacus/agent-runs/20261009-varneb-clamped-PO-E046-r1/source
runtime_root=/public/home/iai806/abacus/agent-runs/20261010-varneb-clamped-p0100-M-E051-r1
parameters=/public/home/iai806/abacus/agent-runs/20261009-varneb-clamped-reversing-minus-E049-r1/factory_parameters.json
[[ ! -e "$runtime_root" ]]
cd "$source_root"
export PYTHONPATH="$PWD" OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
/public/home/iai806/.conda/envs/icu/bin/python - <<'PY'
import json
from pathlib import Path
from examples.hfo2_fixed_input_factory import CONTRACT
from scripts.audit_hfo2_static_replica import sha256
from scripts.relax_clamped_ase_endpoint import load_seed
p=Path('/public/home/iai806/abacus/agent-runs/20261009-varneb-clamped-reversing-minus-E049-r1/factory_parameters.json')
assert sha256(p)=='08a9d81b964f3e7be824d0f4e96935ae801a0beab050d2b28b6b4da903ddc298'
parameters=json.loads(p.read_text())
assert set(parameters)=={'source_directory'}
for name,digest in CONTRACT.items():
    assert sha256(Path(parameters['source_directory'])/name)==digest
seed=Path('benchmarks/hfo2_channels/20261008/clamped_endpoint_seeds/strain_p0100/M/endpoint_seed.json')
assert sha256(seed)=='8fceb0980b99104ffc1ac5dad7a7da6e9856bcfcaf975bb443f4046a9144b1c5'
atoms,boundary,m=load_seed(seed)
assert m['phase_label']=='M' and m['strain']==.01 and m['pressure_gpa']==0.
assert m['physical_contract_sha256']==CONTRACT and m['allow_tilt'] is True
boundary.validate_images([atoms])
po=Path('/public/home/iai806/abacus/agent-runs/20261010-varneb-clamped-p0100-E050-r1/PO_plus/endpoint/endpoint_relax_summary.json')
s=json.loads(po.read_text())
assert s['converged'] and s['physical_contract_sha256']==CONTRACT and s['strain']==.01
assert s['max_atomic_force_eV_per_A']<.03 and s['open_traction_norm_kbar']<2.
print('E051_ZERO_DFT_PREFLIGHT_PASSED: one registered +1% M endpoint; same six physical bytes')
PY
[[ "$(sacct -X -j 28469833 --format=State -n -P | head -n 1)" == COMPLETED ]]
study_count=$(squeue -h -u iai806 --format='%j' | awk '/^hfo2-/ {n++} END {print n+0}')
[[ "$study_count" -le 1 ]]
mkdir "$runtime_root"
handle=$(sbatch --parsable --job-name=hfo2-p01-M-E051 --time=02:00:00 \
  --output="$runtime_root/endpoint-%j.out" --error="$runtime_root/endpoint-%j.err" \
  --export="ALL,RUN_DFT=1,SOURCE_ROOT=$source_root,SEED_MANIFEST=$source_root/benchmarks/hfo2_channels/20261008/clamped_endpoint_seeds/strain_p0100/M/endpoint_seed.json,ENDPOINT_ROOT=$runtime_root/endpoint,PARAMETERS=$parameters,ENDPOINT_STEPS=20" \
  "$source_root/cluster/hf_hfo2_clamped_endpoint_20261008.slurm")
printf 'E051_SUBMITTED M %s\n' "$handle"
date -Is
squeue -h -u iai806 --format='%.18i %.28j %.10T %.12M %.20R'
# End one bounded registered training endpoint; never rerun an unknown submission.

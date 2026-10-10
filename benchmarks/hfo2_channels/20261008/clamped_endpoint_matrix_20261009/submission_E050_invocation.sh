set -euo pipefail
source_root=/public/home/iai806/abacus/agent-runs/20261009-varneb-clamped-PO-E046-r1/source
runtime_root=/public/home/iai806/abacus/agent-runs/20261010-varneb-clamped-p0100-E050-r1
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
params=json.loads(p.read_text())
assert set(params)=={'source_directory'}
physical=Path(params['source_directory'])
for name,digest in CONTRACT.items():
    assert sha256(physical/name)==digest
assert sha256(Path('examples/hfo2_fixed_input_factory.py'))=='b10f73ef5d21027168ef3343efe1314a9faa688955474f3516436a14cdfee7d8'
assert sha256(Path('cluster/hf_hfo2_clamped_endpoint_20261008.slurm'))=='0ea91613ff9e343d297529e9ecc43e0ed117e7eab92a5f2f7b7345c59d685f86'
for phase,digest in [('T','c906303cc9d841b7174a4491c7c663d968d4bcedf72771520d81dc1509d12ac2'),('PO_plus','8d845848b0ee9ced2f26e558204777491e37caadb8a889ff438482fefb7566b6')]:
    seed=Path('benchmarks/hfo2_channels/20261008/clamped_endpoint_seeds/strain_p0100')/phase/'endpoint_seed.json'
    assert sha256(seed)==digest
    atoms,boundary,m=load_seed(seed)
    assert m['phase_label']==phase and m['strain']==.01 and m['pressure_gpa']==0.
    assert m['physical_contract_sha256']==CONTRACT and m['allow_tilt'] is True
    boundary.validate_images([atoms])
po=Path('/public/home/iai806/abacus/agent-runs/20261009-varneb-clamped-PO-E046-r1/endpoint/endpoint_relax_summary.json')
s=json.loads(po.read_text())
assert s['converged'] and s['physical_contract_sha256']==CONTRACT and s['strain']==0.
assert s['max_atomic_force_eV_per_A']<.03 and s['open_traction_norm_kbar']<2.
print('E050_ZERO_DFT_PREFLIGHT_PASSED: two registered +1% endpoints; no holdout or parameter change')
PY
for job in 28464144 28456312 28456313 28453655 28454923; do
    state=$(sacct -X -j "$job" --format=State -n -P | head -n 1)
    [[ "$state" == COMPLETED ]]
    [[ -z "$(squeue -h -j "$job" --format='%i')" ]]
done
[[ -z "$(squeue -h -u iai806 --format='%i %j' | awk '$2 ~ /^hfo2-/ {print $1}')" ]]
mkdir "$runtime_root"
for phase in T PO_plus; do
    mkdir "$runtime_root/$phase"
    handle=$(sbatch --parsable --job-name="hfo2-p01-${phase}-E050" --time=02:00:00 \
      --output="$runtime_root/$phase/endpoint-%j.out" --error="$runtime_root/$phase/endpoint-%j.err" \
      --export="ALL,RUN_DFT=1,SOURCE_ROOT=$source_root,SEED_MANIFEST=$source_root/benchmarks/hfo2_channels/20261008/clamped_endpoint_seeds/strain_p0100/$phase/endpoint_seed.json,ENDPOINT_ROOT=$runtime_root/$phase/endpoint,PARAMETERS=$parameters,ENDPOINT_STEPS=20" \
      "$source_root/cluster/hf_hfo2_clamped_endpoint_20261008.slurm")
    printf 'E050_SUBMITTED %s %s\n' "$phase" "$handle"
done
date -Is
squeue -h -u iai806 --format='%.18i %.28j %.10T %.12M %.20R'
# End manually bounded two registered training endpoint submissions; never rerun an unknown outcome.

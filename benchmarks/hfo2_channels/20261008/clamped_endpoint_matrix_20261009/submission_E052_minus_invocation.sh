set -euo pipefail
source_root=/public/home/iai806/abacus/agent-runs/20261009-varneb-clamped-PO-E046-r1/source
runtime_root=/public/home/iai806/abacus/agent-runs/20261010-varneb-clamped-p0100-minus-E052-r1
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
assert sha256(Path('examples/hfo2_fixed_input_factory.py'))=='b10f73ef5d21027168ef3343efe1314a9faa688955474f3516436a14cdfee7d8'
assert sha256(Path('cluster/hf_hfo2_clamped_endpoint_20261008.slurm'))=='0ea91613ff9e343d297529e9ecc43e0ed117e7eab92a5f2f7b7345c59d685f86'
p=Path('/public/home/iai806/abacus/agent-runs/20261009-varneb-clamped-reversing-minus-E049-r1/factory_parameters.json')
assert sha256(p)=='08a9d81b964f3e7be824d0f4e96935ae801a0beab050d2b28b6b4da903ddc298'
parameters=json.loads(p.read_text())
assert set(parameters)=={'source_directory'}
for name,digest in CONTRACT.items():
    assert sha256(Path(parameters['source_directory'])/name)==digest
seeds={'PO_minus_T_preserving':'f1ffcce75139eafd3cadbd5f5b82463fcba18f84060e4d724b2e8de3167f0c62',
       'PO_minus_T_reversing':'dac285b22799b7960664cdb335c102331d32f432336f168284283a5d8e858ec5'}
for phase,digest in seeds.items():
    seed=Path('benchmarks/hfo2_channels/20261008/clamped_endpoint_seeds/strain_p0100')/phase/'endpoint_seed.json'
    assert sha256(seed)==digest
    atoms,boundary,m=load_seed(seed)
    assert m['phase_label']==phase and m['strain']==.01 and m['pressure_gpa']==0.
    assert m['physical_contract_sha256']==CONTRACT and m['allow_tilt'] is True
    assert atoms.get_chemical_symbols()==['Hf']*4+['O']*8
    boundary.validate_images([atoms])
    print('E052_ZERO_DFT_SEED_PASSED',phase,digest)
po=Path('/public/home/iai806/abacus/agent-runs/20261010-varneb-clamped-p0100-E050-r1/PO_plus/endpoint/endpoint_relax_summary.json')
m=Path('/public/home/iai806/abacus/agent-runs/20261010-varneb-clamped-p0100-M-E051-r1/endpoint/endpoint_relax_summary.json')
for summary in (po,m):
    s=json.loads(summary.read_text())
    assert s['converged'] and s['physical_contract_sha256']==CONTRACT and s['strain']==.01
    assert s['max_atomic_force_eV_per_A']<.03 and s['open_traction_norm_kbar']<2.
print('E052_PREFLIGHT_PASSED: two existing training endpoints; unchanged source and physical bytes; no holdout')
PY
for id in 28469832 28469833 28472493; do
[[ "$(sacct -X -j "$id" --format=State -n -P | head -n 1)" == COMPLETED ]]
done
study_count=$(squeue -h -u iai806 --format='%j' | awk '/^hfo2-/ {n++} END {print n+0}')
[[ "$study_count" -eq 0 ]]
mkdir "$runtime_root"
for phase in PO_minus_T_preserving PO_minus_T_reversing; do
study_count=$(squeue -h -u iai806 --format='%j' | awk '/^hfo2-/ {n++} END {print n+0}')
[[ "$study_count" -le 1 ]]
mkdir "$runtime_root/$phase"
handle=$(sbatch --parsable --job-name="hfo2-p01-$phase-E052" --time=02:00:00 \
  --output="$runtime_root/$phase/endpoint-%j.out" --error="$runtime_root/$phase/endpoint-%j.err" \
  --export="ALL,RUN_DFT=1,SOURCE_ROOT=$source_root,SEED_MANIFEST=$source_root/benchmarks/hfo2_channels/20261008/clamped_endpoint_seeds/strain_p0100/$phase/endpoint_seed.json,ENDPOINT_ROOT=$runtime_root/$phase/endpoint,PARAMETERS=$parameters,ENDPOINT_STEPS=20" \
  "$source_root/cluster/hf_hfo2_clamped_endpoint_20261008.slurm")
mkdir "$runtime_root/$phase/job-$handle"
printf 'E052_SUBMITTED %s %s\n' "$phase" "$handle"
done
date --iso-8601=seconds
squeue -h -u iai806 --format='%.18i %.44j %.10T %.12M %.20R'
# End bounded manual pair; never rerun an unknown submission.

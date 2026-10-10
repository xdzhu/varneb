"""Audit one completed segment and seed its same-chain continuation; zero DFT."""
from datetime import datetime
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tarfile

ROOT = Path('/public/home/iai806/abacus/agent-runs/20261011-varneb-clamped-M-E062-r1')
PRIOR = ROOT.parent/'20261010-varneb-clamped-G2-E054-r1'
SOURCE = PRIOR/'source-fixed'
ARCHIVE_SHA = '170170ed6c3f8749c4930ca8af6a67fc17689dcdf57481e916577e40aea70d9c'
OLD_SCRIPT = SOURCE/'cluster/hf_hfo2_clamped_chain_resume_20261010.slurm'
OLD_SCRIPT_SHA = '79d9990b42a8da90fe915550332d8060761d0ec108551035e284b69f7768ef02'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit():
    if (ROOT/'audit_receipt.json').exists() or (ROOT/'observations').exists():
        raise FileExistsError('single-use audit; inspect prior artifacts, never overwrite')
    row = subprocess.run(['sacct', '-X', '-n', '-P', '-j', '28661020',
        '--format=JobIDRaw,State,ExitCode,AllocCPUS,Partition,User'],
        check=True, capture_output=True, text=True, timeout=30).stdout.strip()
    if row != '28661020|COMPLETED|0:0|32|hfacnormal01|iai806':
        raise ValueError('registered parent completion unproved')
    archive = PRIOR/'source-fixed.tar'
    if sha(archive) != ARCHIVE_SHA or sha(OLD_SCRIPT) != OLD_SCRIPT_SHA:
        raise ValueError('immutable tested runtime changed')
    code = {}
    with tarfile.open(archive) as bundle:
        for member in bundle.getmembers():
            if (member.isfile() and member.name.endswith('.py')
                    and member.name.startswith(('vcneb/', 'scripts/', 'examples/'))):
                content = bundle.extractfile(member).read()
                if content != (SOURCE/member.name).read_bytes():
                    raise ValueError('historical runtime changed: '+member.name)
                code[member.name] = hashlib.sha256(content).hexdigest()
    import numpy as np
    from ase.io import read
    import examples.hfo2_fixed_input_factory as transport
    from scripts.audit_hfo2_static_replica import audited_results
    from scripts.export_hfo2_clamped_observation import export
    from scripts.prepare_hfo2_clamped_resume import prepare
    work = PRIOR/'PO_to_M/band'
    summary = json.loads((work/'vcneb_summary.json').read_text())
    if (summary['status'] != 'max_steps_reached' or summary['converged']
            or (work/'vcneb_failure.json').exists()):
        raise ValueError('completed but unconverged segment required; no blind failed-run restart')
    real_run = subprocess.run
    def forbid(*a, **k):
        raise RuntimeError('audit/cache preflight prohibits external launches')
    subprocess.run = forbid
    try:
        observations = []
        for step in (0, 20):
            out = ROOT/'observations'/f'step_{step:04d}'
            record = export(work, step, out, '28661020', production_script=OLD_SCRIPT)
            observations.append({k:record[k] for k in ('snapshot_step', 'replayed_fmax_eV_A',
                'ordinary_residual_pass', 'sampled_forward_barrier_meV_fu',
                'sampled_reverse_barrier_meV_fu', 'evaluated_chain_sha256')})
            for image in record['raw_image_evaluations']:
                raw = Path(image['raw_source'])
                dest = out/'raw'/f"image_{image['image_index']:04d}"
                (dest/'OUT.ABACUS').mkdir(parents=True)
                for name in ('INPUT', 'KPT', 'STRU', 'OUT.ABACUS/running_scf.log'):
                    shutil.copyfile(raw/name, dest/name)
                shutil.copyfile(Path(image['audit_path']), dest/'source_audit.json')
        calls = sorted(work.glob('image_*/scf_*/call_audit.json'))
        if len(calls) != 140:
            raise ValueError('twenty complete seven-interior updates required')
        call_reports = []
        for path in calls:
            record = json.loads(path.read_text())
            if ({n:sha(path.parent/n) for n in (*transport.CONTRACT, 'STRU')} != record['input_sha256']
                    or any(record['input_sha256'][n] != h for n,h in transport.CONTRACT.items())
                    or sha(path.parent/'OUT.ABACUS/running_scf.log') != record['raw_log_sha256']):
                raise ValueError('native physical bytes or raw log changed')
            raw = audited_results(path.parent)
            if any(not np.allclose(raw[k], record['results'][k], atol=1e-12, rtol=0)
                   for k in ('energy', 'forces', 'stress')):
                raise ValueError('raw native EFS differs from pinned call audit')
            call_reports.append(dict(path=str(path), audit_sha256=sha(path),
                elapsed_seconds=record['elapsed_seconds'], raw_log_sha256=record['raw_log_sha256'],
                input_sha256=record['input_sha256']))
        seed = ROOT/'seed'
        manifest = prepare(ROOT/'observations/step_0020', seed)
        parameters = json.loads((seed/'factory_parameters.json').read_text())
        images = read(seed/'seed.traj', index=':')
        factory = transport.make_clamped_resume_cached_factory(parameters=parameters,
            command='mpirun -np 32 /public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/app/bin/abacus')
        for i, image in enumerate(images):
            old = dict(image.calc.results)
            image.calc = factory(i, image, ROOT/'cache_preflight'/f'image_{i:04d}')
            for key, value in (('energy',image.get_potential_energy()),
                               ('forces',image.get_forces()), ('stress',image.get_stress())):
                if not np.allclose(value, old[key], atol=1e-12, rtol=0) or image.calc.next_call:
                    raise ValueError('all-nine exact current-frame cache replay failed')
    finally:
        subprocess.run = real_run
    receipt = dict(status='audited_same_chain_step20_continuation_ready', check_CST=datetime.now().isoformat(),
        parent_job='28661020', scheduler_row=row, observations=observations, new_DFT_calls=0,
        fresh_parent_interior_SCFs=len(calls), parent_SCF_wall_seconds_sum=sum(r['elapsed_seconds'] for r in call_reports),
        parent_SCF_core_hours_32=sum(r['elapsed_seconds'] for r in call_reports)*32/3600,
        all_parent_calls_original_six_physical_bytes_native_DSIZE32_and_raw_EFS_checked=True,
        current_frame_exact_caches_checked=9, seed_manifest_sha256=sha(seed/'manifest.json'),
        tested_immutable_runtime_archive_sha256=ARCHIVE_SHA, runtime_code_sha256=code,
        analysis_source_sha256={n:sha(ROOT/'audit_source/scripts'/n) for n in
            ('export_hfo2_clamped_observation.py','prepare_hfo2_clamped_resume.py')},
        executed_audit_source_sha256=sha(Path(__file__)), physical_inputs_changed=False,
        new_independent_chains=0, FIRE_state_restored=False, holdout_generated_or_read=False,
        ordinary_converged=False, calls=call_reports)
    with (ROOT/'audit_receipt.json').open('x') as stream:
        json.dump(receipt, stream, indent=2); stream.write('\n')
    # Observable export excludes proprietary PP/orbital/binary files and old runtime trees.
    with tarfile.open(ROOT/'observables.tar', 'x') as bundle:
        for name in ('observations', 'seed', 'cache_preflight', 'audit_receipt.json'):
            bundle.add(ROOT/name, arcname=name)
    print(json.dumps({k:v for k,v in receipt.items() if k not in ('calls','runtime_code_sha256')}))


if __name__ == '__main__':
    audit()

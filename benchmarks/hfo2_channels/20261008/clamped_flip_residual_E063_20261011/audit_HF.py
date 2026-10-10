"""Single-use read-only E063 diagnosis of three completed live-chain frames."""
from datetime import datetime
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tarfile

ROOT = Path('/public/home/iai806/abacus/agent-runs/20261011-varneb-clamped-flip-residual-E063-r1')
PRIOR = ROOT.parent/'20261010-varneb-clamped-G2-E054-r1'
SOURCE = PRIOR/'source-fixed'
STEPS = (13, 14, 15)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit():
    if (ROOT/'observations').exists() or (ROOT/'audit_receipt.json').exists():
        raise FileExistsError('single-use observation namespace; preserve any partial evidence')
    scheduler = subprocess.run(['sacct', '-X', '-n', '-P', '-j', '28661019',
        '--format=JobIDRaw,State,ExitCode,AllocCPUS,Partition,User'],
        check=True, capture_output=True, text=True, timeout=30).stdout.strip()
    fields = scheduler.split('|')
    if (len(fields) != 6 or fields[0] != '28661019' or fields[1] not in ('RUNNING', 'COMPLETED')
            or fields[2:] != ['0:0', '32', 'hfacnormal01', 'iai806']):
        raise ValueError('known 32CPU production source state required')
    archive = PRIOR/'source-fixed.tar'
    if sha(archive) != '170170ed6c3f8749c4930ca8af6a67fc17689dcdf57481e916577e40aea70d9c':
        raise ValueError('registered runtime archive changed')
    checked = 0
    with tarfile.open(archive) as bundle:
        for member in bundle.getmembers():
            if (member.isfile() and member.name.endswith('.py')
                    and member.name.startswith(('vcneb/', 'scripts/', 'examples/'))):
                if bundle.extractfile(member).read() != (SOURCE/member.name).read_bytes():
                    raise ValueError('immutable production runtime changed: '+member.name)
                checked += 1
    from scripts.export_hfo2_clamped_observation import export
    from scripts.analyze_hfo2_clamped_residual import analyze
    script = SOURCE/'cluster/hf_hfo2_clamped_chain_resume_20261010.slurm'
    work = PRIOR/'PO_flip_T_pattern_preserving/band'
    original_run = subprocess.run
    subprocess.run = lambda *a, **k: (_ for _ in ()).throw(RuntimeError('external execution prohibited during diagnosis'))
    observations = []
    try:
        for step in STEPS:
            out = ROOT/'observations'/f'step_{step:04d}'
            report = export(work, step, out, '28661019', production_script=script)
            for image in report['raw_image_evaluations']:
                raw = Path(image['raw_source'])
                dest = out/'raw'/f"image_{image['image_index']:04d}"
                (dest/'OUT.ABACUS').mkdir(parents=True)
                for name in ('INPUT', 'KPT', 'STRU', 'OUT.ABACUS/running_scf.log'):
                    shutil.copyfile(raw/name, dest/name)
                shutil.copyfile(Path(image['audit_path']), dest/'source_audit.json')
            result = analyze(out)
            with (out/'residual.json').open('x') as stream:
                json.dump(result, stream, indent=2); stream.write('\n')
            observations.append({k:result[k] for k in ('snapshot_step','fmax_eV_A',
                'limiting_image_index','limiting_vector_block',
                'perpendicular_only_fmax_same_geometry_eV_A','spring_only_fmax_same_geometry_eV_A')})
    finally:
        subprocess.run = original_run
    receipt = dict(status='three_completed_frames_same_boundary_zero_DFT_diagnosis',
        check_CST=datetime.now().isoformat(), source_job='28661019', scheduler_row=scheduler,
        steps=list(STEPS), observations=observations, new_DFT_calls=0, new_submissions=0,
        physical_inputs_changed=False, production_source_or_jobs_changed=False,
        holdout_generated_or_read=False, G3_bottleneck_selected=False,
        all_exported_native_EFS_and_original_six_physical_bytes_checked=True,
        archived_runtime_python_files_byte_checked=checked, runtime_archive_sha256=sha(archive),
        runtime_scientific_source_sha256={name:sha(SOURCE/name) for name in
            ('vcneb/core.py','vcneb/epitaxial_boundary.py','examples/hfo2_fixed_input_factory.py',
             'scripts/audit_hfo2_static_replica.py')},
        analysis_source_sha256={name:sha(ROOT/'audit_source/scripts'/name) for name in
            ('export_hfo2_clamped_observation.py','analyze_hfo2_clamped_residual.py')},
        executed_audit_source_sha256=sha(Path(__file__)))
    with (ROOT/'audit_receipt.json').open('x') as stream:
        json.dump(receipt, stream, indent=2); stream.write('\n')
    with tarfile.open(ROOT/'observables.tar', 'x') as bundle:
        for name in ('observations','audit_receipt.json'):
            bundle.add(ROOT/name, arcname=name)
    print(json.dumps(receipt))


if __name__ == '__main__':
    audit()

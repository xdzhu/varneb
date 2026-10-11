"""Pinned audit-only adapter for the two accepted E091 children.

Preflight inspects native identity and immutable code, not a terminal result.
Terminal mode refuses live/failed jobs before reading their material outputs.
Neither mode launches DFT, submits/releases a job or changes production code.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shlex
import subprocess
import sys

RUNS = Path('/public/home/iai806/abacus/agent-runs')
PRODUCTION = RUNS/'20261011-varneb-clamped-rolling-E091-r1'
AUDIT_ROOT = RUNS/'20261011-varneb-clamped-child-audit-E093-r1'
SUPPORT = AUDIT_ROOT/'scripts'
SCRIPT = PRODUCTION/'hf_hfo2_clamped_continue_E090_20261011.slurm'
SCRIPT_SHA = 'bf8529bd7f4b182c90f5038aed13059eab157d7cdfc94a2e275abd58d4a45bb5'
SUPPORT_SHA = {
    'audit_hfo2_clamped_terminal.py': '25079147be017329a435749fe94354c6c5f1740889bf2866c21efd40e8b320c0',
    'export_hfo2_clamped_observation.py': 'a8b23b9355bebb0e1909c7d48d6d77584460c83d4d5df6a4153c23edca855ee8',
    'prepare_hfo2_clamped_resume.py': '007f6245cbcf326106cb60037ef18a6c0c9db975420c6b4314097a150cf82cab',
    'analyze_hfo2_clamped_residual.py': 'a2ed0dc5389524b10b77f3329945ab04ad188fe99e57647b574df860a9920d84',
}
REGISTRY = {
    '28727368': dict(parent='28709789', channel='PO_to_M', parent_terminal='terminal', count=338,
        source=RUNS/'20261010-varneb-nested-controls-E055-r1/source',
        archive=RUNS/'varneb-E055-clean-86d2637.tar',
        archive_sha='cbe6552a3fa4237fec2098c08c8056735a9a95df4e18691d09fed6b4f24a613c',
        parent_audit_sha='1a05735749d8893b2c6f992f0ab90f347da9cd65b674a8c1d920227dd8565178',
        seed_sha='172a060107250729ae08e6411d6b3c8d76f6a527350cd7b46acbf2e5a2cd4abd'),
    '28727382': dict(parent='28723655', channel='PO_flip_T_pattern_preserving', parent_terminal='terminal-r2', count=334,
        source=RUNS/'20261010-varneb-clamped-G2-E054-r1/source-fixed',
        archive=RUNS/'20261010-varneb-clamped-G2-E054-r1/source-fixed.tar',
        archive_sha='170170ed6c3f8749c4930ca8af6a67fc17689dcdf57481e916577e40aea70d9c',
        parent_audit_sha='61edc4548d85984e78bbafc29d8b9d366157cfef48ec8175cabf3adfc7ce7e93',
        seed_sha='06729c708ad6b81c210ad65d359a31fa5ccef9ea9cb5d3e2c844084d00c801d0'),
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def paths(registration):
    parent = PRODUCTION/registration['parent']
    terminal = parent/registration['parent_terminal']
    return dict(work=parent/registration['channel']/'band',
        seed=terminal/'seed', parent_audit=terminal/'audit_receipt.json',
        stdout=parent/(registration['channel']+'.slurm.out'),
        stderr=parent/(registration['channel']+'.slurm.err'))


def validate_native_identity(body, job, registration, *, preflight_only=False):
    """Do not transform RUNNING into a synthetic COMPLETED accounting row."""
    values = body.strip().split('|')
    if len(values) != 9:
        raise ValueError('one complete native accounting row required')
    actual_job, state, user, cpus, partition, name, stdout, stderr, submit = values
    p = paths(registration)
    allowed = ('RUNNING', 'COMPLETED') if preflight_only else ('COMPLETED',)
    if (actual_job != job or state not in allowed or user != 'iai806' or cpus != '32'
            or partition != 'hfacnormal01' or name != 'hfo2-G2-E091-'+registration['parent']
            or stdout != str(p['stdout']) or stderr != str(p['stderr'])):
        raise ValueError('registered native identity or required terminal state differs')
    argv = shlex.split(submit)
    exports = [v[len('--export='):] for v in argv if v.startswith('--export=')]
    if not argv or argv[0] != 'sbatch' or argv[-1] != str(SCRIPT) or len(exports) != 1:
        raise ValueError('registered native submitted wrapper required')
    environment = dict(v.split('=', 1) for v in exports[0].split(',') if '=' in v)
    expected = dict(RUN_DFT='1', WORKDIR=str(p['work']), SOURCE_ROOT=str(registration['source']),
        PARENT_JOB_ID=registration['parent'], PARENT_AUDIT_PATH=str(p['parent_audit']),
        PARENT_AUDIT_SHA256=registration['parent_audit_sha'], SEED_ROOT=str(p['seed']),
        SEED_MANIFEST_SHA256=registration['seed_sha'])
    if environment != expected:
        raise ValueError('registered submitted ancestry, source or seed environment differs')
    return state


def load_auditor(registration):
    if sha(SCRIPT) != SCRIPT_SHA:
        raise ValueError('original production wrapper changed')
    for name, digest in SUPPORT_SHA.items():
        if sha(SUPPORT/name) != digest:
            raise ValueError('separately pinned complete audit support changed')
    p = paths(registration)
    if (sha(p['parent_audit']) != registration['parent_audit_sha']
            or sha(p['seed']/'manifest.json') != registration['seed_sha']):
        raise ValueError('accepted parent receipt or seed changed')
    sys.path.insert(0, str(registration['source']))
    import scripts
    scripts.__path__ = [str(SUPPORT), *list(scripts.__path__)]
    from scripts import audit_hfo2_clamped_terminal as auditor
    from scripts import export_hfo2_clamped_observation as exporter
    from scripts import prepare_hfo2_clamped_resume as preparer
    from scripts import analyze_hfo2_clamped_residual as residual
    for module in (auditor, exporter, preparer, residual):
        if sha(Path(module.__file__)) != SUPPORT_SHA[Path(module.__file__).name]:
            raise ValueError('actual imported audit support differs from pinned bytes')
    import vcneb
    import examples.hfo2_fixed_input_factory as factory
    for module in (vcneb, factory):
        if not Path(module.__file__).resolve().is_relative_to(registration['source'].resolve()):
            raise ValueError('production replay must use unchanged original runtime')
    if exporter.REGISTERED_FACTORIES.get(SCRIPT_SHA) != 'examples.hfo2_fixed_input_factory:make_clamped_resume_cached_factory':
        raise ValueError('audit exporter does not recognize the accepted continuation')
    return auditor


def audit_child(job, *, preflight_only=False):
    if job not in REGISTRY:
        raise ValueError('only the two accepted E091 children are admitted')
    registration = REGISTRY[job]
    output = AUDIT_ROOT/job/'terminal'
    if not preflight_only and output.exists():
        raise FileExistsError('inspect an existing terminal audit; never overwrite or replay it')
    body = subprocess.run(['sacct', '-X', '-n', '-P', '-j', job,
        '--format=JobIDRaw,State,User,AllocCPUS,Partition,JobName%100,StdOut%1500,StdErr%1500,SubmitLine%4000'],
        check=True, capture_output=True, text=True, timeout=45).stdout.strip()
    state = validate_native_identity(body, job, registration, preflight_only=preflight_only)
    auditor = load_auditor(registration)
    work = paths(registration)['work']
    if preflight_only:
        code = auditor.verify_runtime(registration['source'], registration['archive'], registration['archive_sha'])
        if len(code) != registration['count']:
            raise ValueError('original immutable runtime file count differs')
        preflight = json.loads((work/'vcneb_preflight.json').read_text())
        if auditor.fresh_interior_call_count(40, preflight, SCRIPT_SHA) != (0, 280):
            raise ValueError('nine-cache continuation requires zero fresh initial SCFs')
        return dict(status='audit_support_ready_not_terminal_material_audit', job_id=job,
            native_state=state, native_accounting_source_raw=body,
            executed_driver_sha256=sha(Path(__file__)), support_sha256=SUPPORT_SHA,
            original_runtime_files_checked=len(code), runtime_archive_sha256=registration['archive_sha'],
            runtime_preflight_sha256=sha(work/'vcneb_preflight.json'),
            production_script_sha256=SCRIPT_SHA, expected_cap_updates=40,
            initial_fresh_SCFs=0, maximum_nominal_fresh_SCFs=280,
            all_native_SCFs_or_final_convergence_audited=False,
            output_files_created=0, new_DFT_calls=0, scheduler_writes=0,
            physical_inputs_or_running_source_changed=False, holdout_generated_or_read=False)
    summary = json.loads((work/'vcneb_summary.json').read_text())
    step = auditor.completed_step(work/'vcneb.opt.log')
    capped = summary.get('status') == 'max_steps_reached'
    if not 0 <= step <= 40 or (capped and step != 40):
        raise ValueError('unexpected finite segment history; inspect instead of silently continuing')
    receipt = auditor.audit(work, job, SCRIPT, registration['source'], registration['archive'],
        registration['archive_sha'], output, outcome='step-cap' if capped else 'converged', seed_latest=capped)
    if len(receipt['runtime_code_sha256']) != registration['count']:
        raise ValueError('audited original runtime file count differs')
    with (output/'driver_receipt.json').open('x') as stream:
        json.dump(dict(job_id=job, parent_job=registration['parent'],
            executed_driver_sha256=sha(Path(__file__)), support_sha256=SUPPORT_SHA,
            terminal_audit_receipt_sha256=sha(output/'audit_receipt.json'),
            new_DFT_calls=0, scheduler_writes=0, new_independent_chains=0,
            physical_inputs_or_running_source_changed=False), stream, indent=2)
    return receipt


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--job-id', choices=tuple(REGISTRY), required=True)
    parser.add_argument('--preflight-only', action='store_true')
    args = parser.parse_args()
    print(json.dumps(audit_child(args.job_id, preflight_only=args.preflight_only), indent=2))

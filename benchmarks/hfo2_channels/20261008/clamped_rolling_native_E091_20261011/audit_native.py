"""Audit either actual registered terminal segment; no DFT/scheduler writes.

The production runtimes and separately pinned audit support remain immutable.
Only a successful unconverged cap receives a latest-frame geometry seed.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys

RUNS = Path('/public/home/iai806/abacus/agent-runs')
ROOT = RUNS/'20261011-varneb-clamped-rolling-E091-r1'
SUPPORT = RUNS/'20261011-varneb-clamped-T-terminal-E075-source-r1/scripts'
PREPARE_SUPPORT = ROOT/'resume_support/scripts'
PREPARE_SHA = '007f6245cbcf326106cb60037ef18a6c0c9db975420c6b4314097a150cf82cab'
SUPPORT_SHA = {
    'audit_hfo2_clamped_terminal.py': '25079147be017329a435749fe94354c6c5f1740889bf2866c21efd40e8b320c0',
    'export_hfo2_clamped_observation.py': '83b41f1f1ac8982db40606e279c5fccba3fac44c786c2238a42d5841cd723cea',
    'analyze_hfo2_clamped_residual.py': 'a2ed0dc5389524b10b77f3329945ab04ad188fe99e57647b574df860a9920d84',
}
REGISTRY = {
    '28723655': dict(cap=40, initial=0, count=334, channel='PO_flip_T_pattern_preserving',audit_subdir='terminal-r2',
        work=RUNS/'20261011-varneb-clamped-flip-E069-r1/PO_flip_T_pattern_preserving/band',
        source=RUNS/'20261010-varneb-clamped-G2-E054-r1/source-fixed',
        script=RUNS/'20261011-varneb-clamped-flip-E069-r1/hf_hfo2_clamped_flip_continue_E069_20261011.slurm',
        script_sha='61199684c543a016a546f02cbfe1eb7ee659c6310db4f7ad0c0340e714647379',
        archive=RUNS/'20261010-varneb-clamped-G2-E054-r1/source-fixed.tar',
        archive_sha='170170ed6c3f8749c4930ca8af6a67fc17689dcdf57481e916577e40aea70d9c'),
    '28709789': dict(cap=10, initial=7, count=338, channel='PO_to_M',audit_subdir='terminal',
        work=RUNS/'20261010-varneb-clamped-G2-E060-r1/PO_to_M/band',
        source=RUNS/'20261010-varneb-nested-controls-E055-r1/source',
        script=RUNS/'20261010-varneb-nested-controls-E055-r1/source/cluster/hf_hfo2_clamped_chain_pilot_20261010.slurm',
        script_sha='7bf764238c4cc857ca2b7e3e94642f0517e1f3b84f4e4c5b17833379ac1101ed',
        archive=RUNS/'varneb-E055-clean-86d2637.tar',
        archive_sha='cbe6552a3fa4237fec2098c08c8056735a9a95df4e18691d09fed6b4f24a613c'),
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit_parent(parent):
    if parent not in REGISTRY:
        raise ValueError('only the two registered G2 parents are admitted')
    r = REGISTRY[parent]
    if sha(r['script']) != r['script_sha']:
        raise ValueError('original production wrapper changed')
    for name, digest in SUPPORT_SHA.items():
        if sha(SUPPORT/name) != digest:
            raise ValueError('separately pinned audit support changed')
    sys.path.insert(0, str(r['source']))
    import scripts
    paths = [str(SUPPORT), *list(scripts.__path__)]
    if parent=='28723655':
        if sha(PREPARE_SUPPORT/'prepare_hfo2_clamped_resume.py')!=PREPARE_SHA:
            raise ValueError('reviewed continuation-aware audit-only preparation changed')
        paths.insert(0,str(PREPARE_SUPPORT))
    scripts.__path__ = paths
    from scripts import audit_hfo2_clamped_terminal as auditor
    from scripts import export_hfo2_clamped_observation as exporter
    from scripts import analyze_hfo2_clamped_residual as residual
    for module in (auditor, exporter, residual):
        if sha(Path(module.__file__)) != SUPPORT_SHA[Path(module.__file__).name]:
            raise ValueError('actual imported auditor differs from pinned support')
    if parent=='28723655':
        from scripts import prepare_hfo2_clamped_resume as prepare_module
        if sha(Path(prepare_module.__file__))!=PREPARE_SHA:
            raise ValueError('actual imported continuation preparation differs')
    summary = json.loads((r['work']/'vcneb_summary.json').read_text())
    capped = summary.get('status') == 'max_steps_reached'
    if capped and auditor.completed_step(r['work']/'vcneb.opt.log') != r['cap']:
        raise ValueError('unexpected cap; inspect rather than silently continue')
    receipt = auditor.audit(r['work'], parent, r['script'], r['source'],
        r['archive'], r['archive_sha'], ROOT/parent/r['audit_subdir'],
        outcome='step-cap' if capped else 'converged', seed_latest=capped)
    if len(receipt['runtime_code_sha256']) != r['count']:
        raise ValueError('original runtime file count differs')
    with (ROOT/parent/'audit_driver_receipt.json').open('x') as stream:
        json.dump(dict(parent_job=parent,terminal_audit_directory=r['audit_subdir'],
            executed_driver_sha256=sha(Path(__file__)),
            executed_prepare_sha256=PREPARE_SHA if parent=='28723655' else None,
            original_runtime_or_physical_inputs_changed=False,new_DFT_calls=0,scheduler_writes=0),stream,indent=2)
    return receipt


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--parent', choices=tuple(REGISTRY), required=True)
    receipt = audit_parent(parser.parse_args().parent)
    print(json.dumps({k: receipt[k] for k in
        ('status','job_id','terminal_step','fresh_interior_SCFs','ordinary_converged')}))

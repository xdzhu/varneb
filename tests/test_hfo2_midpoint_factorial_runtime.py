"""Four real terminal-source decks and finite dispatch guards, not DFT results."""
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

from examples.hfo2_fixed_input_factory import CONTRACT,read_fixed_hfo2_stru
from scripts.audit_hfo2_static_replica import sha256

ROOT=Path(__file__).resolve().parents[1]
CASE=ROOT/'benchmarks/hfo2_channels/20261008/midpoint_factorial_E274_20261019'
spec=importlib.util.spec_from_file_location('E274_runtime_test',CASE/'control.py')
control=importlib.util.module_from_spec(spec);spec.loader.exec_module(control)


def test_only_four_completed_source_crossovers_with_unchanged_electronic_decks():
    prep=json.loads((CASE/'prepared/preparation.json').read_bytes())
    assert prep['upper_stage']=='ordinary_terminal' and prep['additional_SCFs_required']==4
    assert prep['new_DFT_calls']==0 and not prep['execution_admitted']
    assert prep['physical_input_sha256']==CONTRACT
    assert [(r['channel'],r['corner']) for r in prep['decks']]==control.KEYS
    channels={r['channel']:r for r in prep['channels']}
    assert channels['preserving']['sources']['B']['step']==26
    assert channels['preserving']['sources']['B']['source_job']=='28946180'
    assert channels['T']['sources']['B']['step']==2
    for row in prep['decks']:
        folder=CASE/'prepared'/row['path']
        assert {n:sha256(folder/n) for n in ('INPUT','KPT','STRU')}==row['input_sha256']
        assert all(sha256(folder/n)==CONTRACT[n] for n in ('INPUT','KPT'))
        atoms=read_fixed_hfo2_stru(folder/'STRU')
        assert atoms.get_chemical_symbols()==['Hf']*4+['O']*8
        assert np.linalg.det(atoms.cell.array)>0


def test_frozen_crossover_keeps_ordered_internal_pattern_and_only_swaps_cell(tmp_path):
    parent=CASE.with_name('midpoint_reverse_E269_20261018')
    def load(name):
        spec=importlib.util.spec_from_file_location('E274_'+name,parent/f'{name}.py')
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module
    cohort=load('analyze_completed');prepare=load('prepare_factorial')
    exported,_=cohort.restore(tmp_path)
    for channel in ('preserving','T'):
        a,*_=prepare.source_peak(ROOT/prepare.OLD_RESULTS/channel)
        b,*_=prepare.source_peak(exported/channel)
        for corner,internal,cell in (('AB',a,b),('BA',b,a)):
            actual=read_fixed_hfo2_stru(CASE/'prepared'/channel/corner/'STRU')
            # Canonical native STRU wraps coordinates. Keep the ordered analysis
            # lift separately; only integer shifts, never atom remapping, are allowed.
            delta=actual.get_scaled_positions(wrap=False)-internal.get_scaled_positions(wrap=False)
            np.testing.assert_allclose(delta-np.rint(delta),0,atol=1e-12,rtol=0)
            np.testing.assert_allclose(actual.cell.array,cell.cell.array,atol=1e-12,rtol=0)


def test_held_job_must_have_exact_one_node32CPU_one_hour_resources():
    good='JobId=100 JobState=PENDING Reason=JobHeldUser Partition=hfacnormal01 NumCPUs=32 NumTasks=1 CPUs/Task=32 TimeLimit=01:00:00 JobName=hfo2-E274-crossovers NumNodes=1'
    control.validate_held(good,'100')
    for old,new in (('NumCPUs=32','NumCPUs=64'),('TimeLimit=01:00:00','TimeLimit=06:00:00'),
                    ('NumTasks=1','NumTasks=32'),('JobHeldUser','Dependency'),('JobId=100','JobId=101')):
        with pytest.raises(ValueError):control.validate_held(good.replace(old,new),'100')


def test_native_preflight_and_immutable_calculator_writer_match_each_other(tmp_path):
    from examples.hfo2_fixed_input_factory import same_ordered_geometry
    from scripts.prepare_hfo2_channel_work_probes import write_structure
    prep=json.loads((CASE/'prepared/preparation.json').read_bytes())
    for row in prep['decks']:
        atoms=read_fixed_hfo2_stru(CASE/'prepared'/row['path']/'STRU')
        preflight=tmp_path/'native_preflight_STRU';formal=tmp_path/'formal_call_STRU'
        write_structure(preflight,atoms);write_structure(formal,atoms.copy())
        assert sha256(preflight)==sha256(formal)
        assert same_ordered_geometry(atoms,read_fixed_hfo2_stru(formal))


def test_login_or_wrong_rank_thread_environment_cannot_execute_DFT():
    env=dict(RUN_DFT='1',SLURM_JOB_ID='100',SLURM_JOB_NAME='hfo2-E274-crossovers',
        SLURM_JOB_PARTITION='hfacnormal01',SLURM_JOB_NUM_NODES='1',SLURM_NTASKS='1',
        SLURM_CPUS_PER_TASK='32',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
    control.require_resources(env,'100')
    with pytest.raises(ValueError):control.require_resources({},'100')
    for key in env:
        wrong=dict(env);wrong[key]='wrong'
        with pytest.raises(ValueError):control.require_resources(wrong,'100')


def test_entry_is_bounded_to_four_frozen_calls_and_no_CI_or_relaxation():
    shell=(CASE/'run.slurm').read_text(encoding='utf-8')
    assert '#SBATCH --time=01:00:00' in shell and '#SBATCH --cpus-per-task=32' in shell
    assert 'I_MPI_HYDRA_BOOTSTRAP=fork' in shell and 'source /public/home/iai806/apprepo/abacus/' in shell
    assert b'\r' not in (CASE/'control.py').read_bytes() and b'\r' not in (CASE/'run.slurm').read_bytes()
    import ast
    import inspect
    calls=[node.func.attr for node in ast.walk(ast.parse(inspect.getsource(control.run)))
           if isinstance(node,ast.Call) and isinstance(node.func,ast.Attribute)]
    assert calls.count('get_potential_energy')==calls.count('get_forces')==calls.count('get_stress')==1

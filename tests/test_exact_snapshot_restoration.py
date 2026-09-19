import json
import numpy as np
import pytest
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import read,write
from vcneb.executor import ImageEvaluation,ThreadedCalculatorExecutor
from scripts.restore_exact_cached_snapshot import restore,physics_namespace,file_sha


def case(tmp_path):
    p={"n_images":3,"command_line":["test","--image-cache-namespace","binary-fixed"],
       "calculator_parameters":{"encut":455},"licensed_input_fingerprints":{"POTCAR":{"sha256":"frozen"}},"vca":None}
    pre=tmp_path/'pre.json'; pre.write_text(json.dumps(p))
    cache=tmp_path/'cache'
    executor=ThreadedCalculatorExecutor(1,cache_dir=cache,cache_namespace=physics_namespace(p))
    chain=[]
    for i in range(3):
        a=Atoms('CdSe',cell=np.eye(3)*6,positions=[[0,0,0],[2+i*.01,2,2]],pbc=True)
        a.calc=SinglePointCalculator(a,energy=float(i),forces=np.zeros((2,3)),stress=np.zeros(6))
        chain.append(a)
    executor._store_cached(1,chain[1],ImageEvaluation(1.,np.zeros((2,3)),np.zeros((3,3))))
    reference=tmp_path/'ref.traj'; write(reference,chain)
    missing=[a.copy() for a in chain]
    for i in (0,2):
        missing[i].calc=chain[i].calc
    source=tmp_path/'missing.traj'; write(source,missing)
    return source,reference,pre,pre,cache,tmp_path/'archival.traj'


def test_restore_exact_results_without_geometry_or_source_mutation(tmp_path):
    args=case(tmp_path); before=file_sha(args[0])
    r=restore(*args)
    assert r['DFT_calls']==0 and not r['geometry_modified']
    assert file_sha(args[0])==before
    assert r['images'][1]['filled_missing_results']
    a,b=read(args[0],index=':'),read(args[-1],index=':')
    assert all(np.array_equal(x.positions,y.positions) for x,y in zip(a,b))
    assert b[1].get_potential_energy()==1.


def test_wrong_namespace_rejected(tmp_path):
    args=case(tmp_path)
    (args[4]/'cache_metadata.json').write_text(json.dumps({'format_version':1,'namespace':'other'}))
    with pytest.raises(ValueError,match='namespace'):
        restore(*args)
    assert not args[-1].exists()


def test_nearby_geometry_is_not_an_exact_state(tmp_path):
    args=case(tmp_path)
    a=read(args[0],index=':'); a[1].positions[1,0]+=1e-12; write(args[0],a)
    with pytest.raises(ValueError,match='byte-exact'):
        restore(*args)


def test_corrupt_cache_does_not_enter_archival_copy(tmp_path):
    args=case(tmp_path)
    path=next(args[4].glob('*.npz'))
    with path.open('wb') as f:
        np.savez(f,energy=999.,forces=np.zeros((2,3)),stress=np.zeros((3,3)))
    with pytest.raises(ValueError,match='disagree'):
        restore(*args)
    assert not args[-1].exists()


def test_missing_exact_cache_rejected(tmp_path):
    args=case(tmp_path)
    next(args[4].glob('*.npz')).unlink()
    with pytest.raises(FileNotFoundError):
        restore(*args)


def test_output_cannot_replace_source_or_previous_copy(tmp_path):
    args=case(tmp_path)
    with pytest.raises(FileExistsError):
        restore(*args[:-1],args[0])
    restore(*args)
    with pytest.raises(FileExistsError):
        restore(*args)

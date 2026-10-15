"""Replay provenance spelling only after proving exact SCF path identity.

Synthetic E/F/stress below test the auditor, not DFT or material stability.
"""
from pathlib import Path
from types import SimpleNamespace
import hashlib,importlib.util,json,zipfile

import numpy as np
import pytest
from ase.calculators.singlepoint import SinglePointCalculator
from scripts import hfo2_G3_conditional_faces as faces

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'benchmarks/hfo2_channels/20261008'


def module(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m


old=module('original_E183_alias_test',BASE/'G3_released_audit_E183_20261015/audit.py')
new=module('revised_E190_alias_test',BASE/'G3_released_audit_E190_20261015/audit.py')


def test_only_the_provenance_validator_parent_changes_not_scientific_checks():
    previous=(BASE/'G3_released_audit_E183_20261015/audit.py').read_text(encoding='utf-8')
    assert hashlib.sha256(previous.encode()).hexdigest()=='d477fc059b8a38f820254608f093dc0691857dfe7c3b1c0f1afe4d7af71c1393'
    expected='# E190: preserve E183; fix lexical provenance replay after exact resolved-path checks.\n'+previous.replace(
        '    checked=w.study.result_validator(0,atoms,path.parent);e,g=chart.energy_and_gradient(x,atoms)\n',
        '    # Resolved identity was checked above; replay the stored lexical source.\n'
        '    # Do not mistake bundle/../ provenance spelling for a different SCF.\n'
        '    checked=w.study.result_validator(0,atoms,raw.parent);e,g=chart.energy_and_gradient(x,atoms)\n')
    assert (BASE/'G3_released_audit_E190_20261015/audit.py').read_text(encoding='utf-8')==expected


@pytest.fixture
def record(tmp_path):
    d=faces.load_inputs('T');x=np.asarray(d['origin_record']['coordinates'])
    producer=tmp_path/'producer';(producer/'bundle').mkdir(parents=True)
    namespace=producer/'released_curvature/T/evaluations'
    path=namespace/'eval-000000-test/result.json';raw=path.parent/'scf_000000';raw.mkdir(parents=True)
    alias=producer/'bundle/../released_curvature/T/evaluations/eval-000000-test/scf_000000'
    inputs={'INPUT':'synthetic-not-a-DFT-input'}
    (raw/'input_sha256.json').write_text(json.dumps(inputs))
    native=dict(energy=1.,forces=np.zeros((12,3)),stress=np.zeros(6))
    atoms=d['chart'].to_atoms(x);atoms.calc=SinglePointCalculator(atoms,**native)
    e,g=d['chart'].energy_and_gradient(x,atoms)
    distances=atoms.get_all_distances(mic=True);np.fill_diagonal(distances,np.inf)
    def validate(i,atoms,directory):
        return dict(raw_directory=str(directory/'scf_000000'),input_sha256=inputs,genuine_MPI_ranks=32)
    w=SimpleNamespace(study=SimpleNamespace(audited_results=lambda p:native,result_validator=validate))
    row=dict(status='complete',coordinate_sha256=hashlib.sha256(x.tobytes()).hexdigest(),
        coordinates=x.tolist(),contract_sha256=new.physical_contract(d['chart'])['contract_sha256'],
        enthalpy_eV=e,gradient=g.tolist(),stress_eV_per_A3=atoms.get_stress(voigt=False).tolist(),
        volume_A3=atoms.get_volume(),minimum_distance_A=float(distances.min()),maximum_atomic_force_eV_per_A=0.,
        validation=validate(0,atoms,alias.parent))
    path.write_text(json.dumps(row))
    return w,d,x,path,namespace,row,alias


def test_original_bug_reproduces_and_new_reader_preserves_record_bytes(record):
    w,d,x,path,namespace,row,alias=record
    assert str(alias)!=str(path.parent/'scf_000000')
    assert alias.resolve(strict=True)==(path.parent/'scf_000000').resolve(strict=True)
    before=path.read_bytes()
    with pytest.raises(ValueError,match='native energy/gradient/stress/volume/force/distance'):
        old.native_record(w,d,x,path,namespace)
    read,source=new.native_record(w,d,x,path,namespace)
    assert read==row and source['validation']==row['validation']
    assert path.read_bytes()==before


@pytest.mark.parametrize('fault',['energy','gradient','stress','volume','distance','force','MPI','coordinate'])
def test_real_field_or_contract_changes_are_still_rejected(record,fault):
    w,d,x,path,namespace,row,alias=record
    if fault=='energy':row['enthalpy_eV']+=1e-5
    elif fault=='gradient':row['gradient'][0]+=1e-5
    elif fault=='stress':row['stress_eV_per_A3'][0][0]+=1e-5
    elif fault=='volume':row['volume_A3']+=1e-5
    elif fault=='distance':row['minimum_distance_A']+=1e-5
    elif fault=='force':row['maximum_atomic_force_eV_per_A']+=1e-5
    elif fault=='MPI':row['validation']['genuine_MPI_ranks']=1
    else:row['coordinates'][0]=np.nextafter(row['coordinates'][0],1.)
    path.write_text(json.dumps(row))
    with pytest.raises(ValueError):new.native_record(w,d,x,path,namespace)


def test_alias_to_another_real_directory_is_not_accepted(record,tmp_path):
    w,d,x,path,namespace,row,alias=record
    wrong=tmp_path/'other/scf_000000';wrong.mkdir(parents=True)
    row['validation']['raw_directory']=str(wrong);path.write_text(json.dumps(row))
    with pytest.raises(ValueError,match='raw directory differs'):
        new.native_record(w,d,x,path,namespace)


def test_namespace_membership_remains_strict(record,tmp_path):
    w,d,x,path,namespace,row,alias=record
    other=tmp_path/'other_namespace';other.mkdir()
    with pytest.raises(ValueError,match='namespace identity'):
        new.native_record(w,d,x,path,other)


def test_real_first_probe_archive_reproduces_native_energy_force_stress_and_gradient(tmp_path):
    source=BASE/'G3_released_audit_E190_20261015/first_probe_native.zip'
    assert hashlib.sha256(source.read_bytes()).hexdigest()=='d7bd98e90ca97969dd208729739703eba2a00ceed4dc7bb5f902e155e0640e32'
    with zipfile.ZipFile(source) as archive:
        manifest=json.loads(archive.read('source_manifest.json'))
        expected={'result.json','scf_000000/INPUT','scf_000000/KPT','scf_000000/STRU',
                  'scf_000000/input_sha256.json','scf_000000/call_audit.json',
                  'scf_000000/OUT.ABACUS/running_scf.log'}
        assert set(manifest['files'])==expected and set(archive.namelist())==expected|{'source_manifest.json'}
        for name,digest in manifest['files'].items():
            content=archive.read(name);assert hashlib.sha256(content).hexdigest()==digest
            target=tmp_path/name;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(content)
    row=json.loads((tmp_path/'result.json').read_bytes())
    assert hashlib.sha256((tmp_path/'result.json').read_bytes()).hexdigest()=='c35035e2ff2af599e806b0ed76e2dcdb71bc96b9a3df51fa60169e3c044ecf49'
    native=faces.study.audited_results(tmp_path/'scf_000000')
    d=faces.load_inputs('T');atoms=d['chart'].to_atoms(row['coordinates'])
    atoms.calc=SinglePointCalculator(atoms,**native)
    energy,gradient=d['chart'].energy_and_gradient(row['coordinates'],atoms)
    assert energy==row['enthalpy_eV']
    np.testing.assert_allclose(gradient,row['gradient'],atol=2e-12,rtol=0)
    np.testing.assert_array_equal(atoms.get_stress(voigt=False),row['stress_eV_per_A3'])
    files=json.loads((tmp_path/'scf_000000/input_sha256.json').read_bytes())
    for name in ('INPUT','KPT','STRU'):
        assert hashlib.sha256((tmp_path/'scf_000000'/name).read_bytes()).hexdigest()==files[name]
    assert all(files[n]==h for n,h in faces.study.CONTRACT.items())
    assert 'bundle/../released_curvature' in row['validation']['raw_directory']

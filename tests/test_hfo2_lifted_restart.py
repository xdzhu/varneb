import json
from pathlib import Path

import numpy as np
import pytest
from ase.io import read,write

import scripts.prepare_hfo2_lifted_restart as repair
from scripts.prepare_hfo2_switching_chains import ordered_seed
from scripts.audit_hfo2_static_replica import sha256
from vcneb import validate_periodic_path_lift


def test_complete_snapshot_repair_reuses_only_exact_raw_geometry(tmp_path,monkeypatch):
    repository=Path(__file__).resolve().parents[1]/'benchmarks/hfo2_channels/20261008/reference_variants'
    chain=ordered_seed(read(repository/'PO.vasp'),read(repository/'PO_minus_T_preserving_inversion.vasp'))
    work=tmp_path/'stopped'
    snapshot=work/'snapshots/step_0007'
    snapshot.mkdir(parents=True)
    input_file=tmp_path/'original_INPUT'
    input_file.write_text('fixed test Hamiltonian')
    monkeypatch.setattr(repair,'CONTRACT',{'INPUT':sha256(input_file)})
    real_read=repair.read
    monkeypatch.setattr(repair,'read',lambda path,format=None: real_read(path,format='vasp' if format=='abacus' else format))
    monkeypatch.setattr(repair,'audited_results',lambda source: {'energy':-10.,'forces':np.zeros((12,3)),'stress':np.zeros(6)})
    for i,a in enumerate(chain):
        a.set_scaled_positions(a.get_scaled_positions(wrap=True))
        write(snapshot/f'POSCAR_{i:02d}',a,format='vasp',direct=True)
        source=work/f'image_{i:04d}/scf_000007'
        source.mkdir(parents=True)
        write(source/'STRU',a,format='vasp',direct=True)
        (source/'INPUT').write_text(input_file.read_text())
        (source/'OUT.ABACUS').mkdir()
        (source/'OUT.ABACUS/running_scf.log').write_text('complete test output')
        (source/'call_audit.json').write_text(json.dumps({'input_sha256':{n:sha256(source/n) for n in ('INPUT','STRU')},'raw_log_sha256':sha256(source/'OUT.ABACUS/running_scf.log')}))
    output=tmp_path/'repaired'
    report=repair.prepare(work,snapshot,output,'synthetic-no-DFT')
    assert report['n_total_images']==9 and len(report['seed_cache_audits'])==9
    assert not report['DFT_executed'] and not report['old_optimizer_state_reused']
    validate_periodic_path_lift(read(output/'seed.traj',index=':'))
    assert all(a['snapshot_to_raw_ordered_periodic_match'] for a in report['seed_cache_audits'])
    # A stale structure is not accepted merely because it has a complete log.
    stale=chain[3].copy()
    stale.positions[0,0]+=.01
    with pytest.raises(ValueError,match='no complete exact ordered SCF'):
        repair.exact_cached_source(stale,work/'image_0003')
    # Changes to actual evidence fail instead of silently replacing the record.
    (work/'image_0003/scf_000007/OUT.ABACUS/running_scf.log').write_text('changed')
    with pytest.raises(ValueError,match='raw evidence changed'):
        repair.exact_cached_source(chain[3],work/'image_0003')
    with pytest.raises(FileExistsError):
        repair.prepare(work,snapshot,output,'synthetic-no-DFT')


def test_partial_or_unrelated_snapshot_is_rejected_before_output(tmp_path):
    work=tmp_path/'work'
    snapshot=work/'snapshots/step_0003'
    snapshot.mkdir(parents=True)
    with pytest.raises(ValueError,match='complete nine/ten'):
        repair.prepare(work,snapshot,tmp_path/'out','synthetic')
    with pytest.raises(ValueError,match='belong'):
        repair.prepare(work,tmp_path,tmp_path/'out','synthetic')
    assert not (tmp_path/'out').exists()


@pytest.mark.parametrize('label', ['gap_lifted_seed', 'PO_M_lifted_seed'])
def test_actual_hf_force_attribution_replays_from_public_E_F_stress(label):
    # Replays genuine cached hf data, not a new DFT or accelerator benchmark.
    from ase.calculators.singlepoint import SinglePointCalculator
    from vcneb import VCNEB
    root=Path(__file__).resolve().parents[1]/'benchmarks/hfo2_channels/20261008/lift_recovery'/label
    manifest=json.loads((root/'manifest.json').read_text())
    evidence=json.loads((root/'raw_evaluations.json').read_text())
    preflight=json.loads((root/'preflight.json').read_text())
    assert evidence['all_raw_evidence_freshly_verified'] and not evidence['DFT_executed']
    assert sha256(root/'seed.traj')==evidence['seed_traj_sha256']
    assert sha256(root/'manifest.json')==evidence['manifest_sha256']
    lifted=read(root/'seed.traj',index=':')
    original=[a.copy() for a in lifted]
    for a,shift in zip(original,manifest['periodic_lift_audit']['integer_lattice_shifts_by_image_atom']):
        a.set_scaled_positions(a.get_scaled_positions(wrap=False)-np.asarray(shift))
    with pytest.raises(ValueError,match='continuous periodic lift'):
        validate_periodic_path_lift(original)
    validate_periodic_path_lift(lifted)
    for name,images in (('original_broken_lift',original),('explicit_continuous_lift',lifted)):
        for i,(a,raw) in enumerate(zip(images,evidence['points'])):
            assert raw['image_index']==i
            a.calc=SinglePointCalculator(a,energy=raw['energy_eV_cell'],
                forces=np.array(raw['forces_eV_A']),stress=np.array(raw['stress_eV_A3_voigt']))
        fmax=np.linalg.norm(VCNEB(images,k=.2,climb=False).get_forces(),axis=1).max()
        assert fmax==pytest.approx(preflight['force_attribution'][name]['fmax_eV_A'],abs=1e-10)
    assert preflight['not_an_acceleration_claim']

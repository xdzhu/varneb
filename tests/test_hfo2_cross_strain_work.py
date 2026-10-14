"""Native training contrasts plus independent Cartesian work identities."""
import numpy as np
import pytest
from ase.build import bulk
from ase.calculators.emt import EMT
from ase.calculators.singlepoint import SinglePointCalculator

import scripts.analyze_hfo2_cross_strain_work as study
from vcneb import BiaxialClampedModeCoordinates,clamped_plane_vcneb_boundary


def fixture(rotated,tilt):
    atoms=bulk('Cu','fcc',a=3.6,cubic=True)
    cell=atoms.cell.array+[[0.,.1,.15],[.05,0.,.08],[.11,-.07,0.]]
    if rotated:
        u,_=np.linalg.qr(np.random.default_rng(63).normal(size=(3,3)));cell=cell@u.T
    atoms.set_cell(cell,scale_atoms=True)
    boundary=clamped_plane_vcneb_boundary(4,cell,allow_tilt=tilt)
    chart=BiaxialClampedModeCoordinates(atoms,boundary=boundary,cell_scale_A=3.6,
        anchor_strain=0.,strain_scale_A=3.6)
    pair=[]
    for epsilon in (0.,.01):
        y=np.random.default_rng(int(epsilon*100)+91).normal(size=chart.size)*.01
        y[-1]=epsilon*chart.strain_scale_A
        a=chart.to_atoms(y);a.calc=EMT();pair.append(a)
    return chart,pair


@pytest.mark.parametrize('rotated',[False,True])
@pytest.mark.parametrize('tilt',[False,True])
def test_dual_and_independent_physical_work_agree_without_intermediate_sampling(rotated,tilt):
    chart,pair=fixture(rotated,tilt)
    before=[a.positions.copy() for a in pair]
    row=study.cross_pair(chart,*pair)
    assert row['maximum_identity_defect_eV_cell']<1e-12
    assert sum(row['parts_eV_cell'].values())==pytest.approx(row['endpoint_trapezoid_work_eV_cell'],abs=1e-12)
    assert row['raw_energy_change_eV_cell']==pair[1].get_potential_energy()-pair[0].get_potential_energy()
    assert row['energy_minus_endpoint_work_eV_cell']==pytest.approx(row['raw_energy_change_eV_cell']-row['endpoint_trapezoid_work_eV_cell'])
    for a,b in zip(pair,before):np.testing.assert_array_equal(a.positions,b)
    np.testing.assert_allclose(row['control_parameter_values'],[0.,.01],atol=1e-12)


def test_explicit_translation_work_keeps_raw_nonzero_force_sum():
    chart,pair=fixture(False,True)
    for i,a in enumerate(pair):
        a.positions+=i*np.array([.03,.02,-.01])
        a.calc=SinglePointCalculator(a,energy=float(i),forces=np.ones((4,3))*.02,stress=np.zeros(6))
    before=[a.get_forces().copy() for a in pair]
    row=study.cross_pair(chart,*pair)
    assert abs(row['parts_eV_cell']['rigid_translation'])>1e-5
    for a,b in zip(pair,before):np.testing.assert_array_equal(a.get_forces(),b)


def test_long_periodic_lift_is_reported_not_minimum_imaged():
    chart,pair=fixture(False,True)
    a=pair[1];energy=a.get_potential_energy();forces=a.get_forces();stress=a.get_stress()
    a.positions[0]+=a.cell[0]
    a.calc=SinglePointCalculator(a,energy=energy,forces=forces,stress=stress)
    row=study.cross_pair(chart,*pair)
    assert row['cross_condition_lift_change_above_half_cell'] is True
    assert row['geometry_rematched_wrapped_or_projected'] is False


@pytest.fixture(scope='module')
def material():return study.analyze()


def test_all_actual_bands_and_normalizations_are_retained(material):
    assert material['matched_image_pairs']==36 and material['source_frame_records_used']==72
    assert material['distinct_native_log_STRU_identities']==66
    assert material['new_DFT_calls']==material['scheduler_mutations']==0
    assert material['physical_holdout_generated_or_read'] is False
    assert material['intermediate_geometry_generated'] is False
    assert material['barrier_error_bound_stationarity_or_TS_certified'] is False
    assert material['maximum_identity_defect_eV_cell']<1e-12
    assert material['maximum_initial_ordered_lift_serialization_difference_A']<1e-12
    assert material['maximum_common_initial_work_roundoff_eV_cell']<1e-12


@pytest.mark.parametrize('name,change,atomic,control,residual',[
    ('M',-1.2304034162298194,2.7206868091840772,-3.83524744300778,.058484810590644226),
    ('T',12.347026513452874,-.228947111475388,12.27509497966245,.08081112517457217),
    ('preserving',28.5962481552815,36.89528782282829,-12.858887855964417,4.638764068266934),
    ('reversing',-44.02433096583991,-32.33208278039087,.07281392818832971,-11.49506100624212),
])
def test_native_peak_work_contrasts_and_residual_are_not_removed(material,name,change,atomic,control,residual):
    row=next(r for r in material['candidates'] if r['candidate']==name);peak=row['matched_peak_pair']
    assert row['observed_sampled_barrier_change_meV_fu']==pytest.approx(change,abs=1e-9)
    assert peak['relative_PO_parts_meV_fu']['internal_atomic']==pytest.approx(atomic,abs=1e-9)
    assert peak['relative_PO_parts_meV_fu']['prescribed_substrate']==pytest.approx(control,abs=1e-9)
    assert peak['relative_PO_energy_minus_work_meV_fu']==pytest.approx(residual,abs=1e-9)
    assert sum(peak['relative_PO_parts_meV_fu'].values())+residual==pytest.approx(change,abs=1e-9)


def test_missing_training_source_is_refused(monkeypatch):
    monkeypatch.setattr(study,'CASES',study.CASES[1:])
    with pytest.raises(ValueError,match='eight unique'):study.analyze()


def test_common_cached_physics_does_not_authorize_initial_lift_remapping(monkeypatch):
    reader=study.terminal_images
    def altered(*args):
        images,obs,residual,receipt=reader(*args)
        if receipt['job_id']=='28722320':
            a=images[0];e=a.get_potential_energy();f=a.get_forces();s=a.get_stress()
            a.positions[0]+=a.cell[0]
            a.calc=SinglePointCalculator(a,energy=e,forces=f,stress=s)
        return images,obs,residual,receipt
    monkeypatch.setattr(study,'terminal_images',altered)
    with pytest.raises(ValueError,match='different unwrapped lift'):study.analyze()

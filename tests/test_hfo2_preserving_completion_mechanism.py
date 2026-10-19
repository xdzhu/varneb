"""Real zero-strain rows, incomplete-target refusal and explicit figure scope."""
from pathlib import Path
import importlib.util
import json
import pytest

BASE=Path(__file__).resolve().parents[1]/'benchmarks/hfo2_channels/20261008'
CASE=BASE/'preserving_local_completion_E276_20261019'


def module(name):
    s=importlib.util.spec_from_file_location('E276_mechanism_test_'+name,CASE/(name+'.py'))
    m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m


def test_actual_zero_triplet_yields_native_local_energies_and_not_mode_energy_fractions():
    old=json.loads((BASE/'preserving_pattern_branches_E272_20261018/complete_analysis/analysis.json').read_bytes())
    row=next(r for r in old['rows'] if r['index']==0 and r['sign']=='minus');a,j=row['arms']
    e0=a['energy_eV_cell']-a['dE_minus_F_meV_fu']/250
    r=module('mechanism_data').measured_row(0.,'minus',e0,a['energy_eV_cell'],j['energy_eV_cell'],
        a['actual_geometry'],j['actual_geometry'],a['released_gradient_eV_A'],j['released_gradient_eV_A'],.5,'28948425')
    assert r['atomic_lowering_meV_fu']==pytest.approx(13.5992583905)
    assert r['total_lowering_meV_fu']==pytest.approx(13.6170819223)
    assert r['atomic_norm_F_A_A']==pytest.approx(.476962220640566)
    assert r['transverse_proxy_change_F_A_A']==pytest.approx(.12185693858624)
    assert r['oxygen_squared_displacement_fraction_A']==pytest.approx(.803640885869185)
    assert r['open_y_change_F_A_percentage_points']==0.


def test_missing_original_terminal_cannot_generate_a_completed_comparison(tmp_path,monkeypatch):
    m=module('mechanism_data');monkeypatch.setattr(m,'CASE',tmp_path/'missing_native_terminal')
    with pytest.raises(FileNotFoundError):m.build()


def test_missing_triplets_refuse_plot_before_any_output_or_backend_import(tmp_path):
    p=module('plot_mechanism');source=tmp_path/'incomplete.json'
    source.write_text(json.dumps(dict(rows=[],all_four_triplets_A_J_qualified=False)))
    with pytest.raises(ValueError,match='four actual qualified triplets'):p.draw(source,tmp_path/'figure')
    assert not (tmp_path/'figure').exists()


def test_native_atomic_cell_and_four_observable_gates_are_not_relaxed():
    old=json.loads((BASE/'preserving_pattern_branches_E272_20261018/complete_analysis/analysis.json').read_bytes())
    row=old['rows'][0];a,j=row['arms'];e0=a['energy_eV_cell']-a['dE_minus_F_meV_fu']/250
    m=module('mechanism_data');args=[0.,'minus',e0,a['energy_eV_cell'],j['energy_eV_cell'],
        a['actual_geometry'],j['actual_geometry'],.031,j['released_gradient_eV_A'],.5,'28948425']
    with pytest.raises(ValueError,match='raw released-gradient'):m.measured_row(*args)
    args[7]=a['released_gradient_eV_A'];bad=dict(a['actual_geometry'],open_y_stretch_change_percentage_points=.01);args[5]=bad
    with pytest.raises(ValueError,match='entire original F cell'):m.measured_row(*args)


@pytest.mark.parametrize('seed',[BASE/'preserving_pattern_branches_E272_20261018/prepared/seeds/0',CASE/'prepared/seed'])
def test_transverse_mean_proxy_is_not_the_fixed_Qy_shuffle(seed):
    data=module('mechanism_data').transverse_independence(seed)
    assert data['retained_rank']==4 and data['with_transverse_proxy_rank']==5
    assert .60<data['released_transverse_covector_norm']<.62
    assert data['geometric_Qy_not_mean_polar_proxy_y']
    assert data['not_electronic_polarization_or_a_new_eigenmode']


def test_offline_plot_and_mechanism_entry_points_work_outside_repo_cwd(tmp_path):
    import subprocess,sys
    for name in ('mechanism_data.py','plot_mechanism.py'):
        result=subprocess.run([sys.executable,str(CASE/name),'--help'],cwd=tmp_path,capture_output=True,text=True)
        assert result.returncode==0,result.stderr
        assert '--output' in result.stdout


def test_actual_completed_comparison_links_native_local_and_free_path_without_substitution():
    m=module('mechanism_data');d=m.build();rows=d['rows']
    assert d['new_native_calls']==27 and d['new_DFT_calls']==0
    assert d['original_plus1_failure_preserved'] and d['old_and_new_numerical_radii_differ']
    assert d['J_point_to_free_peak_difference_is_not_an_additive_causal_mode_energy']
    assert [(r['strain_percent'],r['sign']) for r in rows]==[(0.,'minus'),(0.,'plus'),(1.,'minus'),(1.,'plus')]
    for r in rows[2:]:
        assert 17.17<r['atomic_lowering_meV_fu']<17.18
        assert .48<r['additional_joint_lowering_meV_fu']<.49
        assert r['open_y_change_F_A_percentage_points']==0.
        assert -.33<r['open_y_change_F_J_percentage_points']<-.31
        assert 16.<r['J_point_minus_connected_free_peak_meV_fu']<16.1
        assert .229<r['Qx_original_F_to_free_peak_change_A']<.230
    assert not d['statistical_replicates_or_inferred_error_bars']


def test_manuscript_contains_native_completion_and_its_nonadditive_mechanism():
    text=(BASE.parents[2]/'paper/VARNEB_JCTC/MANUSCRIPT_DRAFT.md').read_text(encoding='utf-8')
    d=json.loads((CASE/'mechanism_analysis/analysis.json').read_bytes())
    for r in d['rows'][2:]:
        assert f"{r['atomic_lowering_meV_fu']:.3f}" in text
        assert f"{r['additional_joint_lowering_meV_fu']:.3f}" in text
    assert 'Figure 10. Matched preserving F/A/J controls' in text
    assert 'held observables had a displacement scale of 0.810' in text
    assert 'Both original start signs are shown, not averaged or selected by energy' in text
    assert 'not additive\nenergy components' in text
    assert 'original failures as successes' in text
    assert 'remains incomplete at +1%' not in text

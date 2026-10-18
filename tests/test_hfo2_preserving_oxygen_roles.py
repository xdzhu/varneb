"""Preserved reference patterns are not preserved oxygen coordination roles."""
import importlib.util
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
PATH=ROOT/'benchmarks/hfo2_channels/20261008/preserving_pattern_branches_E272_20261018/inspect_oxygen_roles.py'
spec=importlib.util.spec_from_file_location('preserving_oxygen_roles_test',PATH)
reader=importlib.util.module_from_spec(spec);spec.loader.exec_module(reader)


def test_all_three_preserving_endpoint_pairs_exchange_the_same_ordered_roles():
    data=reader.analyze()
    assert [r['substrate_strain'] for r in data['rows']]==[0.,.005,.01]
    for row in data['rows']:
        assert row['all_eight_ordered_oxygen_3_4_roles_exchanged']
        initial,final=row['endpoints']
        assert [r['cutoff_robust_coordination'] for r in initial['oxygen']]==[4,4,3,3,4,4,3,3]
        assert [r['cutoff_robust_coordination'] for r in final['oxygen']]==[3,3,4,4,3,3,4,4]
        for end in (initial,final):
            for atom in end['oxygen']:
                assert len(set(atom['coordination_by_cutoff']))==1
                d3,d4,d5=atom['third_fourth_fifth_Hf_distances_A']
                assert d3<d4<d5
                assert (d4-d3 if atom['cutoff_robust_coordination']==3 else d5-d4)>.9


def test_geometry_only_check_does_not_assign_irrep_electronic_polarity_or_full_path_family():
    data=reader.analyze()
    assert data['geometric_Hf_O_cutoffs_A']==[2.4,2.5,2.6,2.7]
    assert not data['literature_C_N_S_T_family_assigned']
    assert not data['periodic_Hf_plane_crossing_classified']
    assert not data['phonon_irrep_or_electronic_polarization_inferred']
    assert not data['native_EFS_reparsed_in_this_geometry_only_analysis']
    assert not data['atom_permutations_or_DFT_parameter_changes'] and data['new_DFT_calls']==0

from dataclasses import replace

import pytest

from vcneb import ResponseEvaluationCost,recorded_response_dataset_cost


def row(name,**changes):
    args=dict(evaluation_id=name,raw_audit_sha256='a'*64,outcome='completed',
              elapsed_seconds=120.,allocated_cpu_cores=32)
    args.update(changes)
    return ResponseEvaluationCost(**args)


def test_all_shared_unused_and_failed_preparation_counted_once():
    shared=row('reference')
    records=[shared,row('initial'),shared,row('saddle'),row('unused'),row('failed',outcome='failed',elapsed_seconds=60)]
    report=recorded_response_dataset_cost(records,required_evaluation_ids=('initial','saddle','reference'))
    assert report['unique_source_DFT_evaluations']==5 and report['failed_evaluations']==1
    assert report['recorded_wall_seconds_sum']==540
    assert report['recorded_transport_cpu_core_hours_sum']==pytest.approx(4.8)
    assert report['provided_evaluations_not_in_required_subset']==2
    assert report['cost_complete'] and report['new_DFT_calls_by_this_analysis']==0


def test_unknown_time_or_resources_are_not_free_calls():
    records=[row('known'),row('missing-wall',elapsed_seconds=None),row('missing-CPU',allocated_cpu_cores=None)]
    report=recorded_response_dataset_cost(records,required_evaluation_ids=['known'])
    assert report['recorded_wall_seconds_sum'] is None
    assert report['recorded_transport_cpu_core_hours_sum'] is None
    assert report['known_wall_seconds_partial_sum']==240
    assert report['known_transport_core_hours_partial_sum']==pytest.approx(32/30)
    assert not report['cost_complete']


@pytest.mark.parametrize('change',[{'elapsed_seconds':121},{'allocated_cpu_cores':16},
    {'raw_audit_sha256':'b'*64},{'outcome':'failed'}])
def test_conflicting_shared_identity_rejected(change):
    with pytest.raises(ValueError,match='conflicting'):
        recorded_response_dataset_cost([row('same'),row('same',**change)])


def test_missing_required_calls_rejected():
    with pytest.raises(ValueError,match='unrecorded'):
        recorded_response_dataset_cost([row('visible')],required_evaluation_ids=['missing'])


@pytest.mark.parametrize('change',[{'elapsed_seconds':True},{'elapsed_seconds':-1},
    {'elapsed_seconds':float('nan')},{'allocated_cpu_cores':True},{'allocated_cpu_cores':0},
    {'allocated_cpu_cores':1.2},{'raw_audit_sha256':'bad'},{'outcome':'cached'},
    {'evaluation_id':''}])
def test_invalid_records_fail_closed(change):
    with pytest.raises(ValueError): row('original',**change)


def test_empty_analytic_dataset_is_distinct_from_missing_required_DFT():
    report=recorded_response_dataset_cost([])
    assert report['unique_source_DFT_evaluations']==0 and report['cost_complete']
    with pytest.raises(ValueError): recorded_response_dataset_cost([],required_evaluation_ids=['DFT-needed'])


def test_owned_frozen_record():
    r=row('shared')
    assert replace(r,elapsed_seconds=120)==r

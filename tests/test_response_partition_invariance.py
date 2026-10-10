"""A complete harmonic response cannot gain accuracy by re-partitioning Q/R."""
import subprocess
import sys
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from scripts.check_response_partition_invariance import (
    chart_cost_bounds, check_partitions, independent_full_stationary, partition,
)
from scripts.check_nested_response import synthetic_centre
from vcneb import nested_stationary_controls


def test_independent_full_space_partition_and_paired_null():
    r = check_partitions()
    assert r['stationary_points'] == 114 and r['paired_points'] == 54
    assert max(r['maximum_errors'].values()) < 1e-9
    assert len(r['correct_joint_release_refusals']) == 2
    assert r['frozen_vs_joint_max_energy_difference_eV_cell'] > 1e-5


@pytest.mark.parametrize('rotated', [False, True])
def test_promotion_restores_same_full_stationary_model_not_extra_accuracy(rotated):
    c, h, g = synthetic_centre(1, rotated)
    c = partition(c, 'omitted_negative')
    results = {r.control:r for r in nested_stationary_controls(c)}
    assert results['B4_joint_release'].branch is None
    promoted = results['B5_training_instability_promotion']
    assert promoted.promoted_dimension == 1
    assert promoted.branch is not None
    y, e, d, k = independent_full_stationary(c, h, g, .04)
    point = promoted.branch.evaluate(.04)
    np.testing.assert_allclose(point.full_displacement, y, rtol=0., atol=1e-12)
    np.testing.assert_allclose([point.energy_change, point.control_gradient, promoted.branch.control_curvature],
                               [e, d, k], rtol=0., atol=1e-12)


def test_actual_registered_geometry_dimensions_not_a_DFT_hessian():
    rows = chart_cost_bounds()
    assert len(rows) == 10
    for row in rows:
        assert (row['atomic_dimension'],row['open_cell_dimension']) == (33,3)
        assert row['joint_chart_dimension'] == 40
        assert row['full_internal_plus_control_measured_columns'] == 37
        assert row['full_central_pair_one_step_uncached_evaluations'] == 75
        assert row['full_central_pair_two_steps_uncached_evaluations'] == 149
        assert row['actual_new_DFT_calls'] == 0 and not row['actual_bottleneck_or_TS']


def test_cli_refuses_existing_evidence_without_mutation(tmp_path):
    target = tmp_path/'receipt.json'
    target.write_text('previous evidence')
    result = subprocess.run([sys.executable,'-m','scripts.check_response_partition_invariance',
                             '--output',str(target)],capture_output=True,text=True)
    assert result.returncode != 0 and 'fresh evidence namespace' in result.stderr
    assert target.read_text() == 'previous evidence'


def test_historical_HF_material_chart_and_analytic_replay_do_not_claim_predictions():
    case = Path(__file__).resolve().parents[1]/'benchmarks/hfo2_channels/20261008/partition_null_E057_20261010'
    local = json.loads((case/'analytic_clean.json').read_text())
    remote = json.loads((case/'analytic_HF.json').read_text())
    receipt = json.loads((case/'verification_receipt.json').read_text())
    assert local['source_SHA256'] == remote['source_SHA256']
    assert local['registered_geometry_chart_design_counts'] == remote['registered_geometry_chart_design_counts']
    assert remote['analytic']['stationary_points'] == 114 and remote['analytic']['paired_points'] == 54
    assert max(remote['analytic']['maximum_errors'].values()) < 1e-9
    assert not remote['full_Hessian_budget_registered'] and not remote['material_prediction_advantage_proven']
    assert remote['DFT_calls'] == 0 and remote['job_submissions'] == 0 and not remote['held_out_condition_accessed']
    assert receipt['external_executable_prohibited'] and not receipt['HF_pytest_claimed']
    assert not receipt['environment_installed_or_upgraded'] and not receipt['production_sources_or_jobs_mutated']
    for name, expected in receipt['source_SHA256'].items():
        assert hashlib.sha256((case/'executed_source'/name).read_bytes()).hexdigest() == expected


def test_frozen_and_joint_contrast_is_not_zero_or_a_partition_accuracy_claim():
    case = Path(__file__).resolve().parents[1]/'benchmarks/hfo2_channels/20261008/partition_null_E057_20261010'
    remote = json.loads((case/'analytic_HF.json').read_text())
    assert remote['analytic']['frozen_vs_joint_max_energy_difference_eV_cell'] > 1e-5
    refused = remote['analytic']['correct_joint_release_refusals']
    assert len(refused) == 2 and all(r['control'] == 'B4_joint_release' for r in refused)
    assert 'cost evidence' in remote['next_action']

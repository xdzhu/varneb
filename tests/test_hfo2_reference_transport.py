"""Scalar material-adapter evidence; no private author arrays needed in tests."""
import json
from pathlib import Path

import pytest

import scripts.verify_hfo2_reference_transport as study

CASE=Path(__file__).resolve().parents[1]/'benchmarks/hfo2_channels/20261008/reference_transport_E154_20261014'


@pytest.mark.parametrize('channel',['T','preserving'])
def test_actual_registered_reference_bank_keeps_all_frames_ranks_and_measured_coverage(channel):
    report=json.loads((CASE/'scalar_report.json').read_text())
    rows=[r for r in report['rows'] if r['channel']==channel]
    expected={'T_geometric_triplet':3,'T_two_complete_optical_doublets':4,'T_all_optical':33,
              **{f'Cmma_frame_{n}':4 for n in (103,106,114,127)}}
    assert len(rows)==7 and {r['reference']:r['rank'] for r in rows}==expected
    for row in rows:
        assert row['rank']==row['expected_rank']
        assert row['all_cell_and_control_components_exact_zero']
        assert row['reconstruction_max_defect']<1e-13
        assert row['independent_source_local_atomic_work_max_defect']<1e-13
        assert row['internal_plus_translation_work_max_defect']<1e-13
        assert row['measured_span_projection_defect']<1e-12
        assert row['local_translation_overlap_norm']<1e-12
        assert row['excluded_translation_column_Frobenius_norm']>=0
        assert all(not isinstance(value,(dict,list)) for value in row.values())
    assert report['curvature_manifest_sha256']==study.MANIFEST_SHA
    assert report['new_DFT_or_scheduler_writes']==0
    for flag in ('new_raw_SCF_reaudit','private_raw_or_transformed_author_arrays_exported',
        'source_modes_or_forces_ASR_corrected','source_frequencies_used_as_PBE_curvature',
        'conditional_basis_or_forecast_selected','physical_reserved_geometry_generated_or_read'):
        assert report[flag] is False


@pytest.mark.parametrize('existing',[True,False])
def test_scalar_export_never_overwrites_or_writes_into_author_source(tmp_path,monkeypatch,existing):
    source=tmp_path/'private';source.mkdir()
    output=source/'audit.json'
    if existing:output.write_text('original')
    calls=[];monkeypatch.setattr(study,'verify',lambda *a:calls.append(True))
    monkeypatch.setattr('sys.argv',['verify','--source-root',str(source),'--output',str(output)])
    with pytest.raises((FileExistsError,ValueError)):study.main()
    assert calls==[]
    if existing:assert output.read_text()=='original'
    else:assert not output.exists()

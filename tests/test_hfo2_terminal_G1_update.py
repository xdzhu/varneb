from pathlib import Path

import pytest

from scripts.analyze_hfo2_channel_network import read_evaluated_observation
from examples.hfo2_fixed_input_factory import same_ordered_geometry

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "benchmarks/hfo2_channels/20261008"


@pytest.mark.parametrize("name,job,step,n,fmax,peak", [
    ("reversing_step45","28319571",45,9,.09402284566626806,392.8229051971357),
    ("M_refined_step01","28392675",1,12,.09696401269978178,82.67995730602706),
])
def test_actual_terminal_force_and_energy_replay(name,job,step,n,fmax,peak):
    images,r,_ = read_evaluated_observation(BASE / "terminal_G1_update_20261009_17" / name)
    assert len(images) == n and r["source_job_id"] == job and r["snapshot_step"] == step
    assert r["replayed_fmax_eV_A"] == pytest.approx(fmax) and fmax <= .10
    assert max((a.get_potential_energy()+9783.249675811956)*250 for a in images) == pytest.approx(peak)
    assert r["new_DFT_calls"] == 0 and not r["climb"] and r["formula_units"] == 4
    assert all(not p["new_SCF_for_export"] for p in r["raw_image_evaluations"])
    assert any("stationary" in s for s in r["limitations"])


def test_all_four_candidates_keep_common_initial_well_and_pass_ordinary_only():
    folders = [
        ("converged_gap/gap_converged_step06",-1),
        ("switching_converged_update_20261009_1255/preserving_step69",0),
        ("terminal_G1_update_20261009_17/reversing_step45",0),
        ("terminal_G1_update_20261009_17/M_refined_step01",0),
    ]
    well = None
    for folder,index in folders:
        images,r,_ = read_evaluated_observation(BASE / folder)
        assert r["replayed_fmax_eV_A"] <= .10 and r["pressure_GPa"] == 0
        a=images[index]
        if well is not None:
            assert same_ordered_geometry(well,a)
        well=a
        assert a.get_potential_energy() == pytest.approx(-9783.249675811956,abs=1e-10)

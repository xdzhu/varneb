from pathlib import Path
import json
import shutil

import pytest
from ase.io import read

from scripts import analyze_hfo2_G1_peak_sampling as audit
from scripts import prepare_hfo2_M_sampling_refinement as refinement

ROOT = Path(__file__).resolve().parents[1]
CASE = ROOT / "benchmarks/hfo2_channels/20261008/G1_peak_sampling_20261009/completed_HF"


def test_actual_five_raw_statics_distinguish_sampling_and_force_residual(tmp_path):
    r = audit.analyze(ROOT, CASE, tmp_path / "analysis.json")
    t, m = r["channels"]["T_to_PO"], r["channels"]["PO_to_M"]
    assert r["new_DFT_calls_for_analysis"] == 0 and r["unique_MPI_ranks"] == 32
    assert r["new_SCF_calls"] == 5 and r["SCF_core_hours"] == pytest.approx(4.282321521639824)
    assert t["highest_union_sample_meV_fu"] == pytest.approx(115.22152409497721)
    assert t["ordinary_residual_passed"]
    assert m["highest_union_sample_meV_fu"] == pytest.approx(82.74707088730793)
    assert m["new_minus_old_meV_fu"] == pytest.approx(11.164863719841378)
    assert m["replayed_ordinary_fmax_eV_A"] == pytest.approx(.1234081213779406)
    assert not m["ordinary_residual_passed"]
    for c in (t, m):
        assert c["forward_sampled_barrier_meV_fu"]-c["reverse_sampled_barrier_meV_fu"] == pytest.approx(c["endpoint_difference_meV_fu"])
        assert c["optimizer_steps"] == c["additional_SCF_calls"] == 0


@pytest.mark.parametrize("target", ["raw_log", "mpi_probe", "maximum"])
def test_actual_raw_or_summary_changes_fail_closed(tmp_path, target):
    copied = tmp_path / "case"
    shutil.copytree(CASE, copied)
    if target == "raw_log":
        p = copied / "calculations/00/scf_000000/OUT.ABACUS/running_scf.log"
        p.write_bytes(p.read_bytes()+b"\nchanged\n")
    elif target == "mpi_probe":
        (copied / "mpi_affinity.txt").write_text("rank=0 host=one\n"*32)
    else:
        p = copied / "summary.json"
        m = json.loads(p.read_text())
        m["channels"]["PO_to_M"]["highest_new_sample_meV_fu"] += 1
        p.write_text(json.dumps(m))
    with pytest.raises(ValueError):
        audit.analyze(ROOT, copied, tmp_path / "analysis.json")


def completed_job(job):
    # Offline test fixture, not evidence of a live scheduler query.
    return {"job_id":job,"state":"COMPLETED","exit_code":"0:0","end":"2026-10-09T14:56:17"}


def test_actual_twelve_image_seed_is_cache_only_same_lift_and_bounded(tmp_path):
    output = tmp_path / "seed"
    m = refinement.prepare(ROOT, CASE, output, job_reader=completed_job)
    images = read(output / "seed.traj", index=":")
    assert len(images) == m["n_total_images"] == 12 and m["n_active_images"] == 10
    assert all(a.calc is None for a in images)
    p = json.loads((output / "factory_parameters.json").read_text())
    assert len(p["seed_static_directories"]) == 12
    assert [s.split("/calculations/")[1] for s in p["seed_static_directories"][4:7]] == ["02/scf_000000","03/scf_000000","04/scf_000000"]
    assert m["step_cap"] == 20 and m["wall_cap_hours"] == 8
    assert m["maximum_moving_geometry_SCFs"] == 200 and m["new_DFT_calls_for_preparation"] == 0
    assert not m["climb"] and not m["physical_inputs_changed"] and not m["new_independent_channel"]


def test_live_or_failed_sources_and_existing_outputs_refused(tmp_path):
    with pytest.raises(ValueError, match="confirmed completed"):
        refinement.prepare(ROOT, CASE, tmp_path / "seed", job_reader=lambda j:{"job_id":j,"state":"RUNNING","exit_code":"0:0"})
    existing = tmp_path / "existing"
    existing.mkdir()
    with pytest.raises(FileExistsError):
        refinement.prepare(Path("missing"), CASE, existing)
    p = tmp_path / "analysis.json"
    p.write_text("retain")
    with pytest.raises(FileExistsError):
        audit.analyze(Path("missing"), Path("missing"), p)
    assert p.read_text() == "retain"

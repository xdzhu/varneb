"""Contract for the PBS manager plus five 28-rank VASP image workers."""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_235_vasp_template_has_one_manager_and_five_fixed_28_rank_workers() -> None:
    text = (ROOT / "cluster" / "235_batio3_vasp_distributed.pbs").read_text(encoding="utf-8")
    for fragment in (
        "#PBS -q gold5120",
        "#PBS -l nodes=6:ppn=28",
        "PBS_NP:-0}\" == 168",
        "VASP_INITIAL_STATIC",
        "VASP_FINAL_STATIC",
        "--image-workers 5",
        "--initial-static-summary",
        "--final-static-summary",
        "for image in 01 02 03 04 05",
        "scripts/run_235_vasp_image_worker.sh",
        "scripts/validate_vasp_static_baseline_gate.py",
        "--no-climb",
    ):
        assert fragment in text


def test_235_vasp_worker_wrapper_only_accepts_interior_images_and_one_node() -> None:
    text = (ROOT / "scripts" / "run_235_vasp_image_worker.sh").read_text(encoding="utf-8")
    assert "01|02|03|04|05" in text
    assert "mpirun -np 28 -hosts" in text

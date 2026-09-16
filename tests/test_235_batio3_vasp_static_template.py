"""Contract for the PBS gold5120 VASP static baseline."""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_235_bto_vasp_static_template_requests_one_28_core_gold5120_node() -> None:
    text = (ROOT / "cluster" / "235_batio3_vasp_static_baseline.pbs").read_text(encoding="utf-8")
    for fragment in (
        "#PBS -q gold5120",
        "#PBS -l nodes=1:ppn=28",
        "PBS_NP:-0}\" == 28",
        "RUN_DFT:-0",
        "VASP_BIN",
        "VCNEB_GIT_REVISION",
        "mpirun -np 28",
        "--static-only",
        "--n-images 7",
        "scripts/validate_endpoint_identity_gate.py",
    ):
        assert fragment in text

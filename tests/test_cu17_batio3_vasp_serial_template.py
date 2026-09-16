"""Contract for direct serial BTO VASP VCNEB runs on the owned cu17 node."""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_cu17_runner_uses_all_cores_per_serial_interior_image_and_cached_endpoints() -> None:
    text = (ROOT / "cluster" / "cu17_batio3_vasp_serial.sh").read_text(encoding="utf-8")
    for fragment in (
        '"$(hostname -s)" == cu17',
        '"$(nproc)" == 40',
        "mpirun -np 40",
        "--n-images 7",
        "--image-workers 0",
        "--initial-static-summary",
        "--final-static-summary",
        "--no-climb",
        "scripts/validate_vasp_static_baseline_gate.py",
        "scripts/audit_vcneb_result.py",
    ):
        assert fragment in text

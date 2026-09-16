from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_bst50_vca_vcneb_runner_uses_cached_endpoints_and_physical_adapter() -> None:
    text = (ROOT / "cluster" / "cu17_bst50_vca_vcneb.sh").read_text(encoding="utf-8")
    for fragment in (
        '"${RUN_DFT:-0}" == 1',
        '"$(hostname -s)" == cu17',
        "--initial-static-summary",
        "--final-static-summary",
        "--vca-virtual-symbol Ba --vca-components Ba Sr",
        "--no-climb",
        "--cell-interpolation log_strain",
        "scripts/audit_vcneb_result.py",
    ):
        assert fragment in text

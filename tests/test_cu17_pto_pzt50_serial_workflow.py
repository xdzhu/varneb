from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_serial_workflow_orders_pto_before_pzt_and_uses_b_site_vca() -> None:
    text = (ROOT / "cluster" / "cu17_pto_pzt50_serial_workflow.sh").read_text(encoding="utf-8")
    for fragment in (
        '"${RUN_DFT:-0}" == 1', '"$(hostname -s)" == cu17', "wait_or_resume_pto_relax",
        "run_vcneb \"${pto_relax}\" PTO", "--template-a Pb --template-b Ti", "--vca-b-symbol Zr",
        "--virtual-symbol Ti --components Ti Zr", "run_vcneb \"${pzt_relax}\" PZT50 Ti \"Ti Zr\"",
    ):
        assert fragment in text

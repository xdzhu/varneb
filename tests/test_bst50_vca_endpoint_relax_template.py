from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_endpoint_relax_runner_is_guarded_and_symmetry_preserving() -> None:
    runner = (ROOT / "cluster" / "cu17_bst50_vca_endpoint_relax.sh").read_text(encoding="utf-8")
    script = (ROOT / "scripts" / "relax_vasp_vca_endpoint.py").read_text(encoding="utf-8")
    for fragment in ('"${RUN_DFT:-0}" == 1', '"$(hostname -s)" == cu17', "for endpoint in tetragonal cubic", "--fmax 0.03"):
        assert fragment in runner
    for fragment in ("FixSymmetry", "FrechetCellFilter", "CONTCAR.current", "endpoint_relax_summary.json"):
        assert fragment in script

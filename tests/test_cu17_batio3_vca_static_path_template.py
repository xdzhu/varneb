from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_vca_path_runner_is_guarded_and_serial() -> None:
    text = (ROOT / "cluster" / "cu17_batio3_vca_static_path.sh").read_text(encoding="utf-8")
    for fragment in (
        '"${RUN_DFT:-0}" == 1',
        '"$(hostname -s)" == cu17',
        '"$(nproc)" == 40',
        "for index in 00 01 02 03 04 05 06",
        'refusing to overwrite ${image_dir}/OUTCAR',
        'mpirun -np 40 "${vasp_bin}"',
    ):
        assert fragment in text

from pathlib import Path
import importlib.util


MODULE_PATH = Path(__file__).parents[1] / "benchmarks" / "convergence" / "summarize_optimizer_log.py"
SPEC = importlib.util.spec_from_file_location("summarize_optimizer_log", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_parse_and_midnight_elapsed(tmp_path):
    log = tmp_path / "opt.log"
    log.write_text(
        " Step Time Energy fmax\n"
        "FIRE: 0 23:59:50 -1.0 0.20\n"
        "BlockFIRE: 1 00:00:10 -1.1 0.09\n",
        encoding="utf-8",
    )
    rows = MODULE.parse(log)
    assert [row["step"] for row in rows] == [0, 1]
    assert MODULE.elapsed_seconds(rows[0]["time"], rows[1]["time"]) == 20.0


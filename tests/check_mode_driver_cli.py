"""Validate that production drivers expose the documented mode controls."""

from __future__ import annotations

from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    required = (
        "--mode",
        "--mode-guided",
        "--constraint-mode",
        "--mode-amplitude",
        "--mode-envelope",
        "--mode-cell-scale",
    )
    for name in ("run_vcneb_abacus.py", "run_vcneb_vasp.py"):
        result = subprocess.run(
            [sys.executable, str(ROOT / "examples" / name), "--help"],
            check=False,
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            raise SystemExit(f"{name} --help failed: {result.stderr}")
        missing = [option for option in required if option not in result.stdout]
        if missing:
            raise SystemExit(f"{name} is missing mode options: {missing}")

    for name in (
        "hf_batio3_vcneb.slurm",
        "hf_batio3_vcneb_parallel.slurm",
        "hf_hfo2_vcneb.slurm",
        "hf_hfo2_vcneb_distributed.slurm",
    ):
        text = (ROOT / "cluster" / name).read_text(encoding="utf-8")
        for marker in ("MODE_FILE", "MODE_GUIDED", "CONSTRAINT_MODE", "VALIDATE_ONLY", "mode_args"):
            if marker not in text:
                raise SystemExit(f"{name} is missing mode wiring: {marker}")

    print("mode_driver_cli_regression=ok")


if __name__ == "__main__":
    main()

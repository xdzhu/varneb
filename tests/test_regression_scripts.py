"""Pytest entry points for the repository's script-style regression checks.

The checks remain executable as standalone ``python tests/check_*.py``
commands for cluster smoke workflows.  This small adapter makes the same
checks visible to normal ``pytest`` discovery without duplicating their
assertions or launching any DFT calculator.
"""

from __future__ import annotations

from pathlib import Path
import runpy

import pytest


ROOT = Path(__file__).resolve().parents[1]
CHECK_SCRIPTS = (
    "check_abacus_parallel_sort.py",
    "check_hfo2_endpoint_promotion.py",
    "check_image_comparison.py",
    "check_line_search_retry.py",
    "check_mode_template.py",
    "check_native_abacus_driver.py",
    "check_native_endpoint_audit.py",
    "check_release_metadata.py",
    "check_vcneb_audit.py",
    "check_vcneb_forces.py",
    "check_vcneb_metrics.py",
    "check_vcneb_plot.py",
)


@pytest.mark.parametrize("script_name", CHECK_SCRIPTS)
def test_script_style_regression(script_name: str) -> None:
    namespace = runpy.run_path(
        str(ROOT / "tests" / script_name),
        run_name="varneb_pytest_regression",
    )
    namespace["main"]()

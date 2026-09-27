"""Public VCNEB entry points must not silently start with climbing image."""

from __future__ import annotations

import runpy
import sys
from pathlib import Path

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator

from vcneb import VCNEB, run_vcneb


ROOT = Path(__file__).resolve().parents[1]
CLI_DRIVERS = (
    ("run_vcneb_ase.py", ("--initial", "is", "--final", "fs", "--workdir", "work")),
    ("run_vcneb_abacus.py", ("--initial", "is", "--final", "fs")),
    ("run_vcneb_qe.py", ("--initial", "is", "--final", "fs")),
    ("run_vcneb_vasp.py", ()),
)


def _barrier_images() -> list[Atoms]:
    images = []
    for x, energy in ((0.0, 0.0), (1.0, 1.0), (2.0, 0.0)):
        image = Atoms("H", positions=[[x, 0.0, 0.0]], cell=[6.0, 6.0, 6.0], pbc=True)
        image.calc = SinglePointCalculator(
            image, energy=energy, forces=np.zeros((1, 3)), stress=np.zeros(6)
        )
        images.append(image)
    return images


def test_core_defaults_to_ordinary_neb_even_with_an_interior_peak() -> None:
    chain = VCNEB(_barrier_images())
    assert chain.climb is False
    assert VCNEB(_barrier_images(), climb=True).climb is True


def test_run_helper_defaults_to_ordinary_neb() -> None:
    chain, _ = run_vcneb(
        _barrier_images(), steps=0, logfile=None, trajectory=None, snapshot_dir=None
    )
    assert chain.climb is False
    climbed, _ = run_vcneb(
        _barrier_images(), climb=True, steps=0, logfile=None,
        trajectory=None, snapshot_dir=None,
    )
    assert climbed.climb is True


@pytest.mark.parametrize("script,required", CLI_DRIVERS)
def test_material_drivers_require_explicit_ci_opt_in(script, required, monkeypatch) -> None:
    namespace = runpy.run_path(str(ROOT / "examples" / script))

    def enabled(*flags: str) -> bool:
        monkeypatch.setattr(sys, "argv", [script, *required, *flags])
        return namespace["_climb_enabled"](namespace["parse_args"]())

    assert enabled() is False
    assert enabled("--no-climb") is False
    assert enabled("--climb") is True
    if script != "run_vcneb_vasp.py":
        assert enabled("--climb-after", "5") is True
        assert enabled("--no-climb", "--climb-after", "5") is False


def test_generic_slurm_template_keeps_explicit_ci_switch() -> None:
    template = (ROOT / "cluster" / "hf_material_vcneb_ase.slurm").read_text(encoding="utf-8")
    assert 'args+=(--no-climb)' in template
    assert 'args+=(--climb)' in template

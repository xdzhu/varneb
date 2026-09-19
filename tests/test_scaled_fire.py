import numpy as np
from ase import Atoms

from vcneb.core import VCNEB
from vcneb.scaled_fire import ImageScaledFIRE, StagedFIRE


def images(n=9):
    return [Atoms("Ga", positions=[[0, 0, 0]], cell=np.eye(3) * 4, pbc=True) for _ in range(n)]


def uniform_force(chain, magnitude):
    force = np.zeros(chain.ndofs())
    for start in range(0, chain.ndofs(), chain.image_ndofs):
        force[start] = magnitude
    return force.reshape(-1, 3)


def test_image_scaled_fire_removes_sqrt_image_count_cap():
    chain = VCNEB(images(9), climb=False)
    opt = ImageScaledFIRE(chain, logfile=None, maxstep=0.02)
    assert np.isclose(opt.band_scale, np.sqrt(7))
    assert np.isclose(opt.maxstep, 0.02 * np.sqrt(7))


def test_staged_fire_switches_without_second_force_evaluation(monkeypatch):
    chain = VCNEB(images(5), climb=False)
    opt = StagedFIRE(
        chain,
        logfile=None,
        maxstep=0.02,
        switch_fmax=0.12,
        refine_maxstep=0.02,
    )
    calls = []

    def unexpected_force_call():
        calls.append(1)
        raise AssertionError("stage check requested another force evaluation")

    monkeypatch.setattr(chain, "get_forces", unexpected_force_call)
    opt.step(uniform_force(chain, 0.11))
    assert not calls
    assert opt.switched_to_refine
    assert np.isclose(opt.maxstep, 0.02)
    assert opt.switch_history[0]["additional_calculator_evaluations"] == 0


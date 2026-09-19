import importlib.util
from pathlib import Path

from ase import Atoms
import numpy as np
import pytest

spec=importlib.util.spec_from_file_location("gan_failed_step",Path(__file__).resolve().parents[1]/"scripts/prepare_gan_failed_step_canary.py")
module=importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_fractional_step_preserves_unwrapped_motion_and_inputs():
    a=Atoms("GaN",scaled_positions=[[0,0,0],[0.9,0.5,0.5]],cell=np.eye(3)*4,pbc=True)
    b=Atoms("GaN",scaled_positions=[[0,0,0],[1.1,0.5,0.5]],cell=np.eye(3)*4.2,pbc=True)
    before=a.positions.copy()
    candidate=module.fractional_step(a,b,0.5)
    assert np.allclose(candidate.get_scaled_positions(wrap=False)[1],[1.0,0.5,0.5])
    assert np.allclose(candidate.cell.array,np.eye(3)*4.1)
    assert np.array_equal(a.positions,before)


@pytest.mark.parametrize("fraction",[0,1,-0.5,float("nan")])
def test_bad_fraction_rejected(fraction):
    a=Atoms("GaN",scaled_positions=[[0,0,0],[0.5,0.5,0.5]],cell=np.eye(3)*4,pbc=True)
    with pytest.raises(ValueError): module.fractional_step(a,a,fraction)

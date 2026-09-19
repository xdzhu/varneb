"""Result snapshots must preserve executor outputs even on persistent hits."""
import numpy as np
import pytest
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes
from ase.io.trajectory import Trajectory

from vcneb.core import VCNEB, run_vcneb
from vcneb.executor import ThreadedCalculatorExecutor


class ToyCalculator(Calculator):
    implemented_properties=["energy","forces","stress"]

    def __init__(self):
        super().__init__()
        self.calls=0

    def calculate(self,atoms=None,properties=None,system_changes=all_changes):
        super().calculate(atoms,properties,system_changes)
        self.calls+=1
        self.results={"energy":7.0,"forces":np.tile([0,0.3,0.4],(len(atoms),1)),
                      "stress":np.zeros((3,3))}


def images():
    frames=[Atoms("Ga",positions=[[index*0.1,0,0]],cell=np.eye(3)*6,pbc=True) for index in range(4)]
    for frame in frames:frame.calc=ToyCalculator()
    return frames


def test_resumed_cache_hits_are_saved_with_all_results_without_recalculation(tmp_path):
    executor=ThreadedCalculatorExecutor(1,cache_dir=tmp_path/"cache",cache_namespace="toy-only")
    executor.evaluate(images())
    frames=images()
    original_calculators=[frame.calc for frame in frames]
    chain,_=run_vcneb(frames,image_executor=executor,climb=False,steps=0,logfile=None,
                     trajectory=tmp_path/"path.traj",snapshot_dir=tmp_path/"snapshots")
    assert [calculator.calls for calculator in original_calculators]==[1,0,0,1]
    assert [frame.calc for frame in chain.images]==original_calculators
    for path in [tmp_path/"path.traj",tmp_path/"snapshots/chain_step_0000.traj"]:
        with Trajectory(path) as trajectory:
            assert len(trajectory)==4
            for frame in trajectory:
                assert frame.calc is not None
                assert set(frame.calc.results)>={"energy","forces","stress"}
                assert frame.get_potential_energy()==7.0
                assert np.array_equal(frame.get_forces(),[[0,0.3,0.4]])
                assert np.array_equal(frame.get_stress(voigt=False),np.zeros((3,3)))


def test_snapshot_does_not_compute_missing_results(tmp_path):
    frames=images()
    chain=VCNEB(frames,image_executor=ThreadedCalculatorExecutor(1),climb=False)
    with pytest.raises(RuntimeError,match="complete executor evaluation"):
        chain.evaluated_snapshot_images()
    assert all(frame.calc.calls==0 for frame in frames)


def test_external_geometry_change_cannot_serialize_stale_results():
    frames=images()
    chain=VCNEB(frames,image_executor=ThreadedCalculatorExecutor(1),climb=False)
    chain.get_forces()
    calls=[frame.calc.calls for frame in frames]
    frames[1].positions[0,0]+=0.01
    with pytest.raises(RuntimeError,match="geometry changed"):
        chain.evaluated_snapshot_images()
    assert [frame.calc.calls for frame in frames]==calls

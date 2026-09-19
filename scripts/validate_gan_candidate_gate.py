"""Replay the real failed candidate without DFT in the deployed ASE runtime.

The FIRE input force is synthetic delta/dt^2, solely to replay the saved trial;
it is NOT a NEB force or a new scientific path result.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

import ase
from ase.io import read
from ase.io.trajectory import Trajectory
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from vcneb.core import VCNEB
from vcneb.step_control import CandidateStepRejected, CheckedFIRE
from vcneb.vasp_lattice import NativeVaspCandidateValidator, NativeVaspLatticeProbe


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evaluated-trajectory",required=True,type=Path)
    parser.add_argument("--proposed-directory",required=True,type=Path)
    parser.add_argument("--canary-structure",required=True,type=Path)
    parser.add_argument("--native-checker",required=True)
    parser.add_argument("--output",required=True,type=Path)
    args=parser.parse_args()
    if args.output.exists():raise FileExistsError("refuse to overwrite replay evidence")
    with Trajectory(args.evaluated_trajectory) as trajectory:
        if len(trajectory)<29 or len(trajectory)%29:raise ValueError("complete 29-image chains required")
        previous=[trajectory[i].copy() for i in range(len(trajectory)-29,len(trajectory))]
    proposed=[previous[0].copy()]+[read(args.proposed_directory/f"{i:02d}"/"POSCAR") for i in range(1,28)]+[previous[-1].copy()]
    probe=NativeVaspLatticeProbe(args.native_checker)
    validator=NativeVaspCandidateValidator(probe,minimum_distance=1.4,maximum_deformation=0.5)
    chain=VCNEB(previous,climb=False,candidate_validator=validator)
    original=chain.get_x().copy()
    target=VCNEB(proposed,climb=False).get_x()
    cached=object()
    chain._last_evaluations=cached
    try:chain.set_x(target)
    except CandidateStepRejected as error:rejected_details=error.details
    else:raise AssertionError("actual failed candidate was not rejected")
    if not np.array_equal(chain.get_x(),original) or chain._last_evaluations is not cached:
        raise AssertionError("failed candidate modified the evaluated chain/cache")
    delta=target-original
    opt=CheckedFIRE(chain,logfile=None,dt=0.1,maxstep=float(np.linalg.norm(delta))+1,
                    max_candidate_retries=8,candidate_retry_factor=0.5)
    def prohibit_electronic_evaluation():raise AssertionError("replay must not evaluate a physical force")
    chain._compute_forces=prohibit_electronic_evaluation
    opt.step((delta/(0.1**2)).reshape(-1,3))
    expected=original+0.5*delta
    difference=float(np.max(np.abs(chain.get_x()-expected)))
    canary=read(args.canary_structure)
    cell_difference=float(np.max(np.abs(chain.images[10].cell.array-canary.cell.array)))
    position_difference=float(np.max(np.abs(chain.images[10].positions-canary.positions)))
    if difference>1e-12 or cell_difference>1e-12 or position_difference>1e-12:
        raise AssertionError("accepted half-step differs from the full-SCF canary geometry")
    if any(image.calc is not None for image in chain.images):raise AssertionError("DFT calculator attached in replay")
    velocity=opt.vel if hasattr(opt,"vel") else opt.v
    if np.count_nonzero(velocity) or opt.dt!=0.05:raise AssertionError("short-step FIRE momentum/dt reset failed")
    report={"mode":"diagnostic_replay_no_DFT_no_physical_force", "ase_version":ase.__version__,
            "numpy_version":np.__version__,"probe":probe.descriptor(),"all_passed":True,
            "atomic_rejection_preserved_chain_and_cache":True,
            "rejected_details":rejected_details,"accepted_generalized_coordinate_difference":difference,
            "canary_cell_difference_A":cell_difference,"canary_position_difference_A":position_difference,
            "accepted_fraction":0.5,"momentum_reset":True,"dt_after_step":opt.dt,
            "candidate_history":opt.candidate_step_history,
            "parent_sha256":{"evaluated_trajectory":sha256(args.evaluated_trajectory.read_bytes()).hexdigest(),
                             "canary":sha256(args.canary_structure.read_bytes()).hexdigest()},
            "source_sha256":{name:sha256((ROOT/name).read_bytes()).hexdigest() for name in
                             ["vcneb/core.py","vcneb/step_control.py","vcneb/vasp_lattice.py",
                              "scripts/validate_gan_candidate_gate.py"]}}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open("x",encoding="utf-8") as handle:json.dump(report,handle,indent=2)
    print(json.dumps({key:report[key] for key in ["all_passed","ase_version","accepted_fraction",
                     "accepted_generalized_coordinate_difference","canary_cell_difference_A",
                     "canary_position_difference_A"]},indent=2))


if __name__=="__main__":main()

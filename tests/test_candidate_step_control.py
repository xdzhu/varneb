import numpy as np
import pytest
from ase import Atoms
from ase.io import read

from vcneb.core import VCNEB, run_vcneb, validate_candidate_cell_step
from vcneb.step_control import CandidateStepRejected, CheckedFIRE
from vcneb.vasp_lattice import NativeVaspCandidateValidator, LatticeClassification


def images():
    return [Atoms("Ga", positions=[[0, 0, 0]], cell=np.eye(3)*4, pbc=True) for _ in range(4)]


def bounded_validator(candidates):
    if candidates[2].positions[0, 0] > 0.0075:
        raise CandidateStepRejected("late image outside domain", details=[{"image_index": 2}])


def moving_force(chain):
    forces=np.zeros(chain.ndofs())
    forces[0]=1
    forces[chain.image_ndofs]=1
    return forces.reshape(-1, 3)


def test_atomic_rejection_preserves_chain_calculators_endpoints_and_cached_forces():
    chain=VCNEB(images(), candidate_validator=bounded_validator, climb=False)
    cached=object()
    chain._last_evaluations=cached
    original=chain.get_x().copy()
    candidate=original.copy()
    candidate[0]=0.01
    candidate[chain.image_ndofs]=0.01
    with pytest.raises(CandidateStepRejected) as exc:
        chain.set_x(candidate)
    assert np.array_equal(exc.value.candidate_coordinates, candidate)
    assert np.array_equal(chain.get_x(),original)
    assert chain._last_evaluations is cached
    assert all(np.array_equal(image.positions,np.zeros((1,3))) for image in chain.images)


def test_checked_fire_accepts_exact_half_step_without_electronic_call(tmp_path):
    chain=VCNEB(images(),candidate_validator=bounded_validator,climb=False)
    opt=CheckedFIRE(chain,logfile=None,dt=0.1,max_candidate_retries=3,
                    candidate_manifest=tmp_path/"steps.jsonl")
    opt.step(moving_force(chain))
    assert np.isclose(chain.images[1].positions[0,0],0.005)
    assert np.isclose(chain.images[2].positions[0,0],0.005)
    assert opt.nsteps==0  # Only ASE.run increments accepted-step count.
    velocity=opt.vel if hasattr(opt,"vel") else opt.v
    assert np.count_nonzero(velocity)==0
    assert opt.dt==0.05
    assert [entry['event'] for entry in opt.candidate_step_history]==['rejected','accepted_backtracked']
    assert opt.candidate_step_history[-1]['proposal_fraction']==0.5
    assert len((tmp_path/"steps.jsonl").read_text().splitlines())==2
    artifact=opt.candidate_step_history[0]["candidate_trajectory"]
    assert len(read(artifact,index=":"))==len(chain.images)
    assert all(image.calc is None for image in chain.images)


def test_exhaustion_restores_original_optimizer_and_chain():
    def reject(candidates):
        if candidates[2].positions[0,0]>0:
            raise CandidateStepRejected("no valid shortened step")
    chain=VCNEB(images(),candidate_validator=reject,climb=False)
    original=chain.get_x().copy()
    opt=CheckedFIRE(chain,logfile=None,dt=0.1,max_candidate_retries=2)
    with pytest.raises(CandidateStepRejected):opt.step(moving_force(chain))
    assert np.array_equal(chain.get_x(),original)
    assert opt.dt==0.1
    velocity=opt.vel if hasattr(opt,"vel") else opt.v
    assert velocity is None
    assert opt.candidate_step_history[-1]['event']=='exhausted'


def test_unexpected_program_error_is_not_retried():
    calls=[]
    def validator(candidates):
        calls.append(1)
        if candidates[2].positions[0,0]>0:raise OSError("broken native checker")
    chain=VCNEB(images(),candidate_validator=validator,climb=False)
    opt=CheckedFIRE(chain,logfile=None,max_candidate_retries=8)
    with pytest.raises(OSError):opt.step(moving_force(chain))
    assert len(calls)==2 # constructor baseline and one proposal
    assert not opt.candidate_step_history


def test_electronic_failure_is_not_retried(monkeypatch):
    chain=VCNEB(images(),candidate_validator=bounded_validator,climb=False)
    calls=[]
    def electronic_failure():
        calls.append(1)
        raise RuntimeError("VASP electronic calculation failed")
    monkeypatch.setattr(chain,"_compute_forces",electronic_failure)
    opt=CheckedFIRE(chain,logfile=None,max_candidate_retries=8)
    with pytest.raises(RuntimeError,match="electronic calculation failed"):
        opt.step()
    assert len(calls)==1
    assert not opt.candidate_step_history


def test_validator_mutation_cannot_modify_chain():
    chain=VCNEB(images(),climb=False)
    original=chain.get_x().copy()
    def mutate(candidates):candidates[1].positions[0,0]=1
    chain.candidate_validator=mutate
    with pytest.raises(RuntimeError,match="must not modify"):
        chain.set_x(original)
    assert np.array_equal(chain.get_x(),original)


def test_late_invalid_cell_also_cannot_partially_change_chain():
    chain=VCNEB(images(),climb=False)
    original=chain.get_x().copy()
    proposed=original.copy()
    proposed[0]=0.1
    # Deformation coordinates store (F-I)*cell_scale. Set image2's Fxx=-1.
    proposed[chain.image_ndofs+3]=-2*chain.cell_scale
    with pytest.raises(ValueError):chain.set_x(proposed)
    assert np.array_equal(chain.get_x(),original)


def test_native_candidate_checks_all_images_and_serialized_cells():
    class Probe:
        def __init__(self):self.calls=[]
        def classify(self,cell):
            self.calls.append(cell.copy())
            return LatticeClassification(4,11 if cell[0,0]>4.01 else 4,4,"fixture")
    probe=Probe()
    validator=NativeVaspCandidateValidator(probe,maximum_deformation=0.5)
    candidates=images()
    candidates[2].set_cell(np.eye(3)*4.02)
    with pytest.raises(CandidateStepRejected) as exc:validator(candidates)
    assert len(probe.calls)==8
    assert {item['representation'] for item in exc.value.details}=={'manager','POSCAR_preview'}
    assert {item['image_index'] for item in exc.value.details}=={2}


def test_candidate_cell_step_rejects_single_image_volume_blowup():
    previous = images()
    candidates = [image.copy() for image in previous]
    candidates[2].set_cell(np.eye(3) * 4.12, scale_atoms=False)
    with pytest.raises(CandidateStepRejected, match="cell step") as exc:
        validate_candidate_cell_step(previous, candidates, maximum_cell_step=0.05)
    assert exc.value.details[0]["image_index"] == 2
    assert exc.value.details[0]["volume_ratio"] > 1.0


def test_candidate_cell_step_accepts_small_uniform_change():
    previous = images()
    candidates = [image.copy() for image in previous]
    candidates[1].set_cell(np.eye(3) * 4.01, scale_atoms=False)
    validate_candidate_cell_step(previous, candidates, maximum_cell_step=0.05)


@pytest.mark.parametrize("optimizer,validator",[("BFGS",bounded_validator),("FIRE",None)])
def test_backtracking_configuration_rejected_before_any_calculation(optimizer,validator):
    with pytest.raises(ValueError,match="requires"):
        run_vcneb(images(),optimizer=optimizer,candidate_validator=validator,
                  candidate_step_retries=1,validate_calculators=False,trajectory=None)

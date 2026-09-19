import numpy as np
import pytest
from ase import Atoms

from vcneb.block_fire import BlockFIRE
from vcneb.core import VCNEB
from vcneb.step_control import CandidateStepRejected


def images(n=5):
    return [Atoms("Ga", positions=[[0, 0, 0]], cell=np.eye(3) * 4, pbc=True) for _ in range(n)]


def uniform_force(chain):
    force = np.zeros(chain.ndofs())
    for start in range(0, chain.ndofs(), chain.image_ndofs):
        force[start] = 1.0
    return force.reshape(-1, 3)


def test_image_blocks_move_independently_of_band_length():
    chain = VCNEB(images(), climb=False)
    opt = BlockFIRE(chain, logfile=None, dt=1.0, dtmax=1.0, maxstep=0.2)
    opt.step(uniform_force(chain))
    for image in chain.images[1:-1]:
        assert np.isclose(image.positions[0, 0], 0.2)


def test_candidate_backtracking_only_penalizes_rejected_image(tmp_path):
    def validator(candidates):
        if candidates[2].positions[0, 0] > 0.11:
            raise CandidateStepRejected("image 2 outside domain", details=[{"image_index": 2}])

    chain = VCNEB(images(), candidate_validator=validator, climb=False)
    opt = BlockFIRE(
        chain,
        logfile=None,
        dt=1.0,
        dtmax=1.0,
        maxstep=0.2,
        max_candidate_retries=2,
        candidate_manifest=tmp_path / "candidate.jsonl",
    )
    opt.step(uniform_force(chain))
    assert np.isclose(chain.images[1].positions[0, 0], 0.2)
    assert np.isclose(chain.images[2].positions[0, 0], 0.1)
    assert np.isclose(chain.images[3].positions[0, 0], 0.2)
    assert np.isclose(opt.block_dt[0], 1.0)
    assert np.isclose(opt.block_dt[1], 0.5)
    assert np.isclose(opt.block_dt[2], 1.0)
    assert [row["event"] for row in opt.candidate_step_history] == [
        "rejected",
        "accepted_backtracked",
    ]


def test_split_mode_has_atomic_and_cell_blocks_per_image():
    chain = VCNEB(images(4), climb=False)
    opt = BlockFIRE(chain, logfile=None, block_mode="atomic-cell", cell_maxstep=0.03)
    assert len(opt._slices) == 4
    assert np.allclose(opt._caps, [0.2, 0.03, 0.2, 0.03])


@pytest.mark.parametrize("kwargs", [{"maxstep": 0.0}, {"cell_maxstep": -0.1}])
def test_invalid_trust_radius_is_rejected(kwargs):
    chain = VCNEB(images(), climb=False)
    with pytest.raises(ValueError, match="must be positive"):
        BlockFIRE(chain, logfile=None, **kwargs)

import numpy as np
import pytest
from ase import Atoms

from scripts.audit_hfo2_cmma_reference import assign_declared_frame, audit, reconstruct_table_S2


def sample():
    return Atoms("NaCl", cell=np.diag([3., 4., 5.]), scaled_positions=[[.1, .2, .3], [.6, .7, .8]], pbc=True)


def test_proper_frame_and_fixed_permutation_without_mutating_input():
    source = sample()
    before = source.positions.copy()
    frame = np.array([[-1, 0, 0], [0, 0, 1], [0, 1, 0]])
    target = source.copy()
    target.set_cell(frame @ source.cell.array @ frame.T, scale_atoms=False)
    target.set_scaled_positions(source.get_scaled_positions() @ frame.T + [.2, 0, 0])
    result = assign_declared_frame(source, target, frame, frame, [.2, 0, 0], tolerance_A=1e-10)
    assert result["source_to_target_indices"] == [0, 1]
    assert result["maximum_site_error_A"] < 1e-12
    np.testing.assert_array_equal(source.positions, before)


@pytest.mark.parametrize("bad", [np.diag([-1., 1., 1.]), np.eye(3)*1.01, np.eye(2)])
def test_improper_or_strained_frame_rejected(bad):
    with pytest.raises(ValueError, match="proper"):
        assign_declared_frame(sample(), sample(), np.eye(3), bad, [0, 0, 0], tolerance_A=1e-4)


def test_inconsistent_cell_and_species_are_not_fit():
    target = sample()
    target.cell[0, 0] += .01
    with pytest.raises(ValueError, match="cells disagree"):
        assign_declared_frame(sample(), target, np.eye(3), np.eye(3), [0, 0, 0], tolerance_A=1e-4)
    with pytest.raises(ValueError, match="composition"):
        assign_declared_frame(sample(), Atoms("NaNa", cell=np.eye(3), pbc=True), np.eye(3), np.eye(3), [0, 0, 0], tolerance_A=1e-4)


def test_ambiguous_assignment_rejected():
    source = Atoms("NaNa", cell=np.eye(3), scaled_positions=[[0, 0, 0], [.5, 0, 0]], pbc=True)
    target = Atoms("NaNa", cell=np.eye(3), scaled_positions=[[.25, 0, 0], [.75, 0, 0]], pbc=True)
    with pytest.raises(ValueError, match="ambiguous"):
        assign_declared_frame(source, target, np.eye(3), np.eye(3), [0, 0, 0], tolerance_A=.3)


def test_primitive_reconstruction_is_same_reference_not_new_production_cell():
    atoms, record = reconstruct_table_S2()
    assert len(atoms) == 12 and record["primitive_atoms"] == 6
    assert record["determinant"] == 2 and record["primitive_input_unchanged"]
    assert sorted(atoms.get_chemical_symbols()) == ["Hf"]*4+["O"]*8
    assert abs(np.linalg.det(record["Cartesian_rotation"]) - 1) < 1e-12
    assert np.linalg.det(atoms.cell.array) > 0


def test_changed_public_source_is_rejected_before_parsing(tmp_path):
    (tmp_path / "bto.in").write_text("not the pinned author data")
    with pytest.raises(ValueError, match="Git blob changed"):
        audit(tmp_path)

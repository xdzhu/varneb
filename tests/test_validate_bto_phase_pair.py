from ase import Atoms

from scripts.validate_bto_phase_pair import audit_bto_phase_pair


def _bto(*, tetragonal: bool) -> Atoms:
    cell = [4.0, 4.0, 4.1 if tetragonal else 4.0]
    ti_z = 0.53 if tetragonal else 0.5
    return Atoms(
        symbols=["Ba", "Ti", "O", "O", "O"],
        scaled_positions=[
            [0.0, 0.0, 0.0],
            [0.5, 0.5, ti_z],
            [0.5, 0.5, 0.0],
            [0.5, 0.0, 0.5],
            [0.0, 0.5, 0.5],
        ],
        cell=cell,
        pbc=True,
    )


def test_distinct_tetragonal_and_cubic_bto_pass() -> None:
    result = audit_bto_phase_pair(_bto(tetragonal=True), _bto(tetragonal=False))
    assert result["status"] == "passed"
    assert result["tetragonal"]["c_over_mean_ab"] > 1.005
    assert result["tetragonal"]["ti_equatorial_o_z_offset_A"] > 0.01


def test_cubic_input_mislabelled_as_tetragonal_is_rejected() -> None:
    result = audit_bto_phase_pair(_bto(tetragonal=False), _bto(tetragonal=False))
    assert result["status"] == "failed"
    assert any("tetragonal c/a" in issue for issue in result["issues"])
    assert any("Ti-O polar" in issue for issue in result["issues"])

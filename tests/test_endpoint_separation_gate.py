"""A converged zero-length path is not evidence for a phase transition."""

from ase import Atoms
import pytest

from vcneb import path_geometry_diagnostics, validate_path_geometry


def _pair(displacement: float) -> list[Atoms]:
    initial = Atoms("H", positions=[[0.0, 0.0, 0.0]], cell=[4.0, 4.0, 4.0], pbc=True)
    final = initial.copy()
    final.positions[0, 0] = displacement
    return [initial, final]


def test_nearly_identical_endpoints_fail_declared_separation_gate() -> None:
    report = path_geometry_diagnostics(_pair(0.0002), minimum_endpoint_separation=0.05)
    assert report["endpoint_separation_A"] == pytest.approx(0.0002)
    assert report["valid"] is False
    assert "same structure" in report["issues"][0]
    with pytest.raises(ValueError, match="same structure"):
        validate_path_geometry(_pair(0.0002), minimum_endpoint_separation=0.05)


def test_distinct_endpoints_pass_declared_separation_gate() -> None:
    report = validate_path_geometry(_pair(0.2), minimum_endpoint_separation=0.05)
    assert report["endpoint_separation_A"] == pytest.approx(0.2)
    assert report["valid"] is True


def test_endpoint_separation_threshold_must_be_positive() -> None:
    with pytest.raises(ValueError, match="minimum_endpoint_separation"):
        path_geometry_diagnostics(_pair(0.2), minimum_endpoint_separation=0.0)

from pathlib import Path
import hashlib

import numpy as np
import pytest

from vcneb.qe_modes import read_qe_gamma_modes


def fixture_text(masses=(2., 7.)):
    # Analytic normal modes, deliberately not a material result.
    rng = np.random.default_rng(721)
    eigenvectors, _ = np.linalg.qr(rng.normal(size=(6, 6)))
    displacement = eigenvectors / np.repeat(np.sqrt(masses), 3)[:, None]
    displacement /= np.linalg.norm(displacement, axis=0)
    lines = ["diagonalizing the dynamical matrix ...", "q = 0 0 0", "********"]
    for mode in range(6):
        frequency = mode - 2.
        lines.append(f"freq ({mode+1}) = {frequency:.9f} [THz] = {frequency*33.35640951981521:.9f} [cm-1]")
        phase = np.exp(1j * .17 * mode)
        for row in (displacement[:, mode] * phase).reshape(2, 3):
            lines.append("( " + " ".join(f"{z.real:.9f} {z.imag:.9f}" for z in row) + " )")
    return "\n".join(lines) + "\n", eigenvectors


def write_fixture(tmp_path, text=None):
    source = tmp_path / "matdyn.modes"
    source.write_text(fixture_text()[0] if text is None else text, encoding="ascii")
    return source


def test_physical_displacements_are_not_eigenvectors_and_complex_phase_is_retained(tmp_path):
    text, expected = fixture_text()
    source = write_fixture(tmp_path, text)
    result = read_qe_gamma_modes(source, [2., 7.])
    np.testing.assert_allclose(np.abs(expected.T @ result.mass_weighted_eigenvectors), np.eye(6), atol=2e-9)
    assert np.max(np.abs(result.cartesian_displacements.conj().T @ result.cartesian_displacements - np.eye(6))) > .1
    assert np.max(np.abs(result.mass_weighted_eigenvectors.imag)) > .1
    assert result.mass_weighted_gram_max_defect < 2e-9
    assert result.source_sha256 == hashlib.sha256(source.read_bytes()).hexdigest()
    assert not result.mass_weighted_eigenvectors.flags.writeable
    assert not result.masses_amu.flags.writeable


def test_later_nonzero_q_is_not_mislabelled_gamma_or_parsed(tmp_path):
    text, _ = fixture_text()
    result = read_qe_gamma_modes(write_fixture(tmp_path, text + "diagonalizing the dynamical matrix ...\nq = 0 .5 0\nnot parsed\n"), [2., 7.])
    assert len(result.frequencies_THz) == 6


@pytest.mark.parametrize("mutate", [
    lambda t: t.replace("q = 0 0 0", "q = .5 0 0"),
    lambda t: t.replace("q = 0 0 0", "q = nan 0 0"),
    lambda t: t.replace("freq (2)", "freq (1)"),
    lambda t: t.replace("[cm-1]", "[Ry]", 1),
    lambda t: t.replace("-66.712819040", "-65.712819040"),
    lambda t: t.rsplit("\n", 3)[0],
    lambda t: t.replace("( ", "( nan ", 1),
    lambda t: t.replace("( ", "( 1e999 ", 1),
    lambda t: t.replace("q = 0 0 0", "q = 0 0 0\nq = 0 0 0"),
    lambda t: t.replace("-2.000000000 [THz] = -66.712819040", "1e999 [THz] = 1e999"),
    lambda t: t + "freq (7) = 1 [THz] = 33.356409520 [cm-1]\n(1 0 0 0 0 0)\n(1 0 0 0 0 0)\n",
])
def test_malformed_or_incomplete_mode_block_is_rejected(tmp_path, mutate):
    with pytest.raises(ValueError):
        read_qe_gamma_modes(write_fixture(tmp_path, mutate(fixture_text()[0])), [2., 7.])


@pytest.mark.parametrize("masses", [[1., 1.], [2.], [0., 7.], [np.nan, 7.], [[2., 7.]], []])
def test_wrong_order_count_or_mass_convention_is_rejected(tmp_path, masses):
    with pytest.raises(ValueError):
        read_qe_gamma_modes(write_fixture(tmp_path), masses)


@pytest.mark.parametrize("tolerance", [0, -1, np.nan, 1])
def test_invalid_tolerance_is_rejected(tmp_path, tolerance):
    with pytest.raises(ValueError):
        read_qe_gamma_modes(write_fixture(tmp_path), [2., 7.], gram_tolerance=tolerance)


def test_non_normalized_displacements_and_non_ascii_are_rejected(tmp_path):
    text, _ = fixture_text()
    source = write_fixture(tmp_path, text)
    rows = text.splitlines()
    rows[4] = "( 0 0 0 0 0 0 )"
    with pytest.raises(ValueError, match="normalized"):
        read_qe_gamma_modes(write_fixture(tmp_path, "\n".join(rows)), [2., 7.])
    source.write_bytes(b"\xff")
    with pytest.raises(UnicodeDecodeError):
        read_qe_gamma_modes(source, [2., 7.])


def test_file_size_bound(tmp_path, monkeypatch):
    source = write_fixture(tmp_path)
    monkeypatch.setattr(Path, "stat", lambda self: type("Stat", (), {"st_size": 13*1024**2})())
    with pytest.raises(ValueError, match="12MiB"):
        read_qe_gamma_modes(source, [2., 7.])


def test_common_mass_rescaling_preserves_eigenvectors(tmp_path):
    source = write_fixture(tmp_path)
    normal = read_qe_gamma_modes(source, [2., 7.])
    large = read_qe_gamma_modes(source, [2e300, 7e300])
    np.testing.assert_allclose(normal.mass_weighted_eigenvectors, large.mass_weighted_eigenvectors, atol=1e-14)

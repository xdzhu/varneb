import numpy as np
import pytest
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator

from vcneb import VCNEB, minimum_image_path_lift, validate_periodic_path_lift


def chain(xs, pbc=True):
    return [Atoms('HHe', scaled_positions=[[x, .1, .2], [.4, .6, .7]],
                  cell=np.eye(3)*5, pbc=pbc) for x in xs]


def test_explicit_lift_preserves_geometry_order_and_input():
    originals = chain([.95, .02, .08])
    positions = [a.positions.copy() for a in originals]
    with pytest.raises(ValueError, match='continuous periodic lift'):
        validate_periodic_path_lift(originals)
    images, report = minimum_image_path_lift(originals)
    np.testing.assert_allclose([a.get_scaled_positions(wrap=False)[0, 0] for a in images], [.95, 1.02, 1.08])
    assert not report['physical_geometry_changed'] and not report['atomic_permutation_applied']
    assert report['integer_lattice_shifts_by_image_atom'][1][0] == [1, 0, 0]
    for original, image, before in zip(originals, images, positions):
        np.testing.assert_array_equal(original.positions, before)
        assert original.get_chemical_symbols() == image.get_chemical_symbols()
        delta = image.get_scaled_positions(wrap=False)-original.get_scaled_positions(wrap=False)
        np.testing.assert_allclose(delta, np.rint(delta), atol=1e-12)
        assert image.calc is None
    validate_periodic_path_lift(images)


def test_adjacent_lift_keeps_a_resolved_winding_instead_of_endpoint_shortcut():
    images, _ = minimum_image_path_lift(chain([.1, .4, .7, 0.]))
    np.testing.assert_allclose([a.get_scaled_positions(wrap=False)[0, 0] for a in images], [.1, .4, .7, 1.])
    validate_periodic_path_lift(images)


def test_nonperiodic_coordinates_are_never_shifted():
    images, report = minimum_image_path_lift(chain([.1, 1.1], pbc=[False, True, True]))
    np.testing.assert_allclose(images[-1].positions, chain([1.1], pbc=[False, True, True])[-1].positions)
    assert not np.any(report['integer_lattice_shifts_by_image_atom'])


def test_ambiguous_half_step_and_atom_pbc_changes_are_rejected():
    with pytest.raises(ValueError, match='half-cell'):
        minimum_image_path_lift(chain([.1, .6]))
    images=chain([.1, .2])
    images[-1].pbc[2]=False
    with pytest.raises(ValueError, match='ordered atoms and PBC'):
        minimum_image_path_lift(images)
    images=chain([.1, .2])
    images[-1]=images[-1][[1, 0]]
    with pytest.raises(ValueError, match='ordered atoms and PBC'):
        minimum_image_path_lift(images)


def test_neb_forces_are_invariant_after_explicit_integer_gauge_recovery():
    canonical=chain([.95, 1.02, 1.08, 1.14])
    scrambled=[a.copy() for a in canonical]
    for i,a in enumerate(scrambled[1:], start=1):
        a.set_scaled_positions(a.get_scaled_positions(wrap=False)+np.array([[i, -i, 2*i], [-i, 2*i, i]]))
    recovered,_=minimum_image_path_lift(scrambled)
    for family in (canonical,recovered):
        for i,a in enumerate(family):
            a.calc=SinglePointCalculator(a, energy=float((i-1)**2),
                forces=np.array([[.1,.2,.3],[-.1,-.2,-.3]]),stress=np.zeros(6))
    np.testing.assert_allclose(VCNEB(canonical,k=.2).get_forces(), VCNEB(recovered,k=.2).get_forces(), atol=1e-12)


def test_cli_guard_rejects_wrapped_supplied_chain_before_calculator_loading(tmp_path):
    from ase.io import write
    from vcneb.material_runner import main
    images=chain([.95,.02,.08])
    write(tmp_path/'initial.vasp', images[0])
    write(tmp_path/'final.vasp', images[-1])
    write(tmp_path/'seed.traj', images)
    def forbidden_loader(_):
        pytest.fail('calculator must not be loaded for an invalid periodic lift')
    with pytest.raises(ValueError,match='continuous periodic lift'):
        main(['--initial',str(tmp_path/'initial.vasp'),'--final',str(tmp_path/'final.vasp'),
            '--initial-chain',str(tmp_path/'seed.traj'),'--n-images','3','--workdir',str(tmp_path/'run'),
            '--factory','invalid:forbidden','--require-continuous-periodic-lift'],symbol_loader=forbidden_loader)
    assert not (tmp_path/'run/vcneb_preflight.json').exists()

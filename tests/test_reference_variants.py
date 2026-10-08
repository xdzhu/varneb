import numpy as np
import pytest
from ase import Atoms

from vcneb.reference_variants import apply_parent_operation


def scaffold():
    parent = Atoms("Si2", scaled_positions=[[.25, .25, .25], [.75, .75, .75]],
                   cell=[5, 5, 6], pbc=True)
    product = parent.copy()
    product.positions += [[.2, .1, .3], [-.1, .2, .3]]
    return parent, product


def test_inversion_has_parent_not_product_mapping_and_is_involution():
    parent, product = scaffold()
    first = apply_parent_operation(parent, product, -np.eye(3), [0, 0, 0])
    assert first.source_to_target == (1, 0)
    delta = product.positions - parent.positions
    actual = first.atoms.positions - parent.positions
    assert np.allclose(actual, -delta[::-1])
    second = apply_parent_operation(parent, first.atoms, -np.eye(3), [0, 0, 0])
    assert np.allclose(second.atoms.positions, product.positions)
    assert np.allclose(product.positions - parent.positions, delta)


def test_identity_preserves_order_without_nearest_product_relabeling():
    parent, product = scaffold()
    product.positions[0] += [2.4, 0, 0]
    result = apply_parent_operation(parent, product, np.eye(3), [0, 0, 0])
    assert result.source_to_target == (0, 1)
    assert np.allclose(result.atoms.get_scaled_positions(), product.get_scaled_positions())


@pytest.mark.parametrize("rotation,translation", [
    (np.eye(3), [.1, 0, 0]), (np.eye(3) * .5, [0, 0, 0]),
    (np.eye(2), [0, 0, 0]), (np.eye(3), [0, np.nan, 0]),
])
def test_invalid_parent_operation_rejected(rotation, translation):
    with pytest.raises(ValueError):
        apply_parent_operation(*scaffold(), rotation, translation)


def test_species_and_cell_orientation_are_not_silently_changed():
    parent, product = scaffold()
    parent[1].symbol = "O"
    product[1].symbol = "O"
    with pytest.raises(ValueError, match="species-preserving"):
        apply_parent_operation(parent, product, -np.eye(3), [0, 0, 0])
    parent, product = scaffold()
    with pytest.raises(ValueError, match="metric"):
        apply_parent_operation(parent, product, [[0, 0, 1], [0, 1, 0], [1, 0, 0]], [0, 0, 0])

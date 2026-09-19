import numpy as np
import pytest
from ase import Atoms
from scripts.prepare_cdse_sheppard_suite import atomic_endpoint, mapped_chain, symmetry


def endpoints():
    rs = Atoms("Cd4Se4", cell=[[11.46,0,0],[2.865,2.865,0],[0,0,5.73]],
               scaled_positions=[[0,0,0],[.5,0,0],[.25,0,.5],[.75,0,.5],[.25,0,0],[.75,0,0],[0,0,.5],[.5,0,.5]], pbc=True)
    wz = Atoms("Cd4Se4", cell=[[8.76,0,0],[2.19,np.sqrt(3)*2.19,0],[0,0,6.96]],
               scaled_positions=[[0,0,0],[.5,0,0],[1/6,1/3,.5],[2/3,1/3,.5],[1/6,1/3,.117],[2/3,1/3,.117],[0,0,.617],[.5,0,.617]], pbc=True)
    return rs, wz


def test_seed_phase_identity():
    rs, wz = endpoints()
    assert all(r["number"] == 225 for r in symmetry(rs))
    assert all(r["number"] == 186 for r in symmetry(wz))


def test_equivalent_basis_preserves_periodic_geometry():
    rs, _ = endpoints()
    changed = atomic_endpoint(rs)
    assert np.isclose(changed.get_volume(), rs.get_volume())
    assert np.allclose(np.sort(changed.get_all_distances(mic=True).ravel()), np.sort(rs.get_all_distances(mic=True).ravel()))
    assert all(r["number"] == 225 for r in symmetry(changed))
    assert np.isclose(changed.cell.lengths()[0], np.sqrt(5)*rs.cell.lengths()[1])


@pytest.mark.parametrize("route", ["cell_mapping", "atomic_mapping"])
def test_explicit_chain_geometry_and_counts(route):
    images = mapped_chain(*endpoints(), route)
    assert len(images) == 17
    assert all(a.get_chemical_symbols() == ["Cd"]*4+["Se"]*4 for a in images)
    assert all(np.linalg.det(a.cell.array) > 0 for a in images)
    assert all(r["number"] == 186 for r in symmetry(images[-1]))


def test_two_mappings_are_distinct_not_duplicate_jobs():
    rs, wz = endpoints()
    a = mapped_chain(rs,wz,"atomic_mapping")
    c = mapped_chain(rs,wz,"cell_mapping")
    assert not np.allclose(a[0].cell.array, c[0].cell.array)
    assert not np.allclose(a[8].positions, c[8].positions)


def test_bad_species_and_route_rejected():
    rs,wz = endpoints()
    with pytest.raises(ValueError):
        mapped_chain(rs,wz,"unknown")
    with pytest.raises(ValueError):
        mapped_chain(rs[:4],wz,"cell_mapping")

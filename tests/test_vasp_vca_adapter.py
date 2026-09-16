from __future__ import annotations

from pathlib import Path

import numpy as np
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes
from ase.calculators.vasp import Vasp

from vcneb.vasp import (
    ExplicitPotcarVasp,
    VirtualCrystalCalculator,
    expand_virtual_site,
    potcar_dataset_labels,
    potcar_setups,
)


class ExpandedFakeCalculator(Calculator):
    implemented_properties = ["energy", "forces", "stress", "free_energy"]

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        assert atoms.get_chemical_symbols() == ["Ba", "Sr", "Ti", "O", "O", "O"]
        self.results = {
            "energy": -3.0,
            "free_energy": -3.1,
            "forces": np.arange(18, dtype=float).reshape(6, 3),
            "stress": np.arange(6, dtype=float),
        }


def test_potcar_labels_and_ase_setups_preserve_sv_choices(tmp_path: Path) -> None:
    potcar = tmp_path / "POTCAR"
    potcar.write_text(
        " TITEL = PAW_PBE Ba_sv 06Sep2000\n TITEL = PAW_PBE Sr_sv 07Sep2000\n"
        " TITEL = PAW_PBE Ti_sv 26Sep2005\n TITEL = PAW_PBE O 08Apr2002\n",
        encoding="latin-1",
    )
    assert potcar_dataset_labels(potcar) == ["Ba_sv", "Sr_sv", "Ti_sv", "O"]
    assert potcar_setups(potcar) == {"Ba": "_sv", "Sr": "_sv", "Ti": "_sv", "O": ""}


def test_explicit_potcar_is_restored_after_ase_writes_inputs(tmp_path: Path, monkeypatch) -> None:
    source = tmp_path / "source.POTCAR"
    source.write_bytes(b"approved licensed bytes")
    directory = tmp_path / "image"
    directory.mkdir()

    def fake_write_input(self, atoms, properties=None, system_changes=None):
        (Path(self.directory) / "POTCAR").write_bytes(b"ASE replacement")

    monkeypatch.setattr(Vasp, "write_input", fake_write_input)
    calculator = ExplicitPotcarVasp(source_potcar=source, directory=str(directory))
    calculator.write_input(Atoms("H"))
    assert (directory / "POTCAR").read_bytes() == source.read_bytes()


def test_vca_adapter_sums_component_forces_back_to_one_physical_site(tmp_path: Path) -> None:
    atoms = Atoms("BaTiO3", scaled_positions=[[0, 0, 0]] * 5, cell=[4, 4, 4], pbc=True)
    expanded, mapping = expand_virtual_site(atoms, virtual_symbol="Ba", components=("Ba", "Sr"))
    assert expanded.get_chemical_symbols() == ["Ba", "Sr", "Ti", "O", "O", "O"]
    assert mapping == [[0, 1], [2], [3], [4], [5]]

    base = ExpandedFakeCalculator()
    base.directory = str(tmp_path / "00")
    atoms.calc = VirtualCrystalCalculator(base, virtual_symbol="Ba", components=("Ba", "Sr"))
    expected = np.arange(18, dtype=float).reshape(6, 3)
    forces = atoms.get_forces()
    np.testing.assert_allclose(forces[0], expected[0] + expected[1])
    np.testing.assert_allclose(forces[1:], expected[2:])
    assert atoms.get_potential_energy() == -3.0
    np.testing.assert_allclose(atoms.get_stress(), np.arange(6, dtype=float))

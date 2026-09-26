"""The VASP diagnostic must alter only ENCUT, not GaN production inputs."""

import json

from ase import Atoms
from ase.io import write

from scripts.prepare_gan_strain_cutoff_diagnostic import stage


def test_isolated_seven_static_cutoff_probe(tmp_path) -> None:
    source = tmp_path / "source"
    template = source / "cases/axis12_plus"
    template.mkdir(parents=True)
    center = Atoms("Ga2N2", cell=[3.2, 3.2, 4.8],
                   scaled_positions=[[0, 0, 0], [0.5, 0.5, 0.5],
                                     [0.25, 0.25, 0.25], [0.75, 0.75, 0.75]], pbc=True)
    write(source / "center_POSCAR", center, format="vasp", direct=True, sort=False)
    (template / "INCAR").write_text(
        " ENCUT = 600.000000\n EDIFF = 1e-7\n ISYM = -1\n", encoding="utf-8",
    )
    (template / "KPOINTS").write_text("Gamma 8x8x6", encoding="utf-8")
    (template / "POTCAR").write_text("mock-for-test-only", encoding="utf-8")
    names = [f"axis{i:02d}_{side}" for i in range(18) for side in ("plus", "minus")]
    for name in names:
        directory = source / "cases" / name
        directory.mkdir(exist_ok=True)
        write(directory / "POSCAR", center, format="vasp", direct=True, sort=False)
    (source / "manifest.json").write_text(json.dumps({
        "status": "inputs_finalized_no_DFT", "step_A": 0.02,
        "coordinate_count": 18, "pressure_GPa": 45.7, "cell_scale_A": 3.4,
        "cases": [{"name": name} for name in names],
        "source_sha256": {"trajectory": "trajectory", "summary": "summary"},
        "source_static_sha256": {"OUTCAR": "outcar"},
    }), encoding="utf-8")
    work = tmp_path / "cutoff800"
    manifest = stage(source, work, 800)
    assert len(manifest["cases"]) == 7
    assert manifest["encut_eV"] == 800 and manifest["baseline_encut_eV"] == 600
    assert (work / "cases/center/INCAR").read_text(encoding="utf-8") == (
        " ENCUT = 800.000000\n EDIFF = 1e-7\n ISYM = -1\n"
    )
    assert (template / "INCAR").read_text(encoding="utf-8").startswith(" ENCUT = 600")
    assert all(record["name"] in ("center", "axis12_plus", "axis12_minus",
                                  "axis13_plus", "axis13_minus",
                                  "axis14_plus", "axis14_minus")
               for record in manifest["cases"])

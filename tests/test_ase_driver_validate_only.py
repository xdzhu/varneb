import json
import sys

from ase import Atoms
from ase.io import write

from examples import run_vcneb_ase


def test_validate_only_does_not_instantiate_external_calculators(tmp_path, monkeypatch) -> None:
    initial = Atoms("H", scaled_positions=[[0.0, 0.0, 0.0]], cell=[4.0, 4.0, 4.0], pbc=True)
    final = initial.copy()
    final.cell = [4.1, 4.0, 4.0]
    initial_path = tmp_path / "initial.vasp"
    final_path = tmp_path / "final.vasp"
    write(initial_path, initial, format="vasp")
    write(final_path, final, format="vasp")

    def builder(**kwargs):
        def forbidden_factory(*args, **factory_kwargs):
            raise AssertionError("validate-only instantiated an external calculator")

        return forbidden_factory

    monkeypatch.setattr(run_vcneb_ase, "_load_symbol", lambda spec: builder)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "run_vcneb_ase.py",
            "--initial",
            str(initial_path),
            "--final",
            str(final_path),
            "--workdir",
            str(tmp_path / "run"),
            "--factory",
            "test:builder",
            "--n-images",
            "3",
            "--cell-interpolation",
            "linear",
            "--no-align-cells",
            "--mapping",
            "identity",
            "--validate-only",
        ],
    )
    run_vcneb_ase.main()
    payload = json.loads((tmp_path / "run" / "vcneb_preflight.json").read_text())
    assert payload["calculator_validation"] == "factory_configuration_only"
    assert payload["calculator_reports"] == []

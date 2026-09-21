from pathlib import Path


TEMPLATE = Path(__file__).parents[1] / "cluster" / "hf_material_vcneb_ase.slurm"


def test_cp2k_shell_default_is_single_rank() -> None:
    text = TEMPLATE.read_text(encoding="utf-8")
    assert 'ntasks=1 --ntasks-per-node=1 cp2k_shell.psmp' in text
    assert 'ntasks=${image_mpi} --ntasks-per-node=${image_mpi} cp2k_shell.psmp' not in text


def test_production_requires_endpoint_static_gate() -> None:
    text = TEMPLATE.read_text(encoding="utf-8")
    assert 'ENDPOINT_STATIC_SUMMARY is required before production VC-NEB' in text
    assert 'scripts/validate_ase_static_gate.py' in text
    assert 'STATIC_ONLY:-0' in text

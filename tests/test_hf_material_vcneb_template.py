from pathlib import Path
import json


TEMPLATE = Path(__file__).parents[1] / "cluster" / "hf_material_vcneb_ase.slurm"
ENDPOINT_TEMPLATE = Path(__file__).parents[1] / "cluster" / "hf_ase_endpoint_relax.slurm"


def test_cp2k_shell_default_is_single_rank() -> None:
    text = TEMPLATE.read_text(encoding="utf-8")
    assert 'ntasks=1 --ntasks-per-node=1 cp2k_shell.psmp' in text
    assert 'ntasks=${image_mpi} --ntasks-per-node=${image_mpi} cp2k_shell.psmp' not in text


def test_production_requires_endpoint_static_gate() -> None:
    text = TEMPLATE.read_text(encoding="utf-8")
    assert 'ENDPOINT_STATIC_SUMMARY is required before production VC-NEB' in text
    assert 'scripts/validate_ase_static_gate.py' in text
    assert 'STATIC_ONLY:-0' in text


def test_endpoint_template_uses_single_rank_cp2k_shell() -> None:
    text = ENDPOINT_TEMPLATE.read_text(encoding="utf-8")
    assert "scripts/relax_ase_endpoint.py" in text
    assert 'ntasks=1 --ntasks-per-node=1 cp2k_shell.psmp' in text
    assert 'RUN_DFT:-0' in text


def test_abinit_templates_use_matching_intel_hydra_launcher() -> None:
    production = TEMPLATE.read_text(encoding="utf-8")
    endpoint = ENDPOINT_TEMPLATE.read_text(encoding="utf-8")
    assert 'module load compiler/intel/2017.5.239' in production
    assert 'module load mpi/intelmpi/2017.4.239' in production
    assert 'mpiexec.hydra -bootstrap slurm -n ${image_mpi} abinit' in production
    assert 'unset I_MPI_PMI_LIBRARY I_MPI_HYDRA_BOOTSTRAP_EXEC_EXTRA_ARGS I_MPI_HYDRA_BOOTSTRAP' in production
    assert 'srun --exclusive --nodes=1 --ntasks=${image_mpi} --ntasks-per-node=${image_mpi} abinit' not in production
    assert 'module load compiler/intel/2017.5.239' in endpoint
    assert 'module load mpi/intelmpi/2017.4.239' in endpoint
    assert 'mpiexec.hydra -bootstrap slurm -n 32 abinit' in endpoint
    assert 'unset I_MPI_PMI_LIBRARY I_MPI_HYDRA_BOOTSTRAP_EXEC_EXTRA_ARGS I_MPI_HYDRA_BOOTSTRAP' in endpoint
    assert 'srun --exclusive --nodes=1 --ntasks=32 --ntasks-per-node=32 abinit' not in endpoint


def test_gan_cp2k_stable_profile_uses_conservative_mixing() -> None:
    profile = json.loads(
        (Path(__file__).parents[1] / "examples" / "material_profiles" / "cp2k_gan_pbe_dzvp_stable.json").read_text()
    )
    assert profile["max_scf"] == 1000
    assert "DIRECT_P_MIXING" in profile["inp"]
    assert "ALPHA 0.05" in profile["inp"]

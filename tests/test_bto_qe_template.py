"""Guard the BTO QE production and lightweight-preflight templates."""

from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_bto_qe_production_template_preserves_the_five_by_32_worker_contract() -> None:
    text = (ROOT / "cluster" / "hf_batio3_vcneb_qe_distributed.slurm").read_text(encoding="utf-8")
    for fragment in (
        "#SBATCH -N 2",
        "#SBATCH --ntasks=160",
        "n_interior=$((n_images - 2))",
        "image_workers=${IMAGE_WORKERS:-${n_interior}}",
        "image_mpi=${IMAGE_MPI:-32}",
        "five 32-MPI interior workers (5 x 32 = 160)",
        "--n-images \"${n_images}\"",
        "--no-climb",
        "if [[ \"${RUN_DFT:-0}\" != 1 ]]",
        "Refusing a 160-rank allocation without RUN_DFT=1",
        "${QE_ENV_SCRIPT:?Set QE_ENV_SCRIPT",
        "${ENDPOINT_IDENTITY_GATE:?Set ENDPOINT_IDENTITY_GATE",
        "scripts/validate_endpoint_identity_gate.py",
        "QE_PP_MANIFEST",
        "--pp-manifest \"${pp_manifest}\"",
        "--cell-interpolation log_strain --mapping auto --align-translation",
    ):
        assert fragment in text


def test_bto_qe_preflight_template_is_one_rank_and_never_launches_pw_x() -> None:
    text = (ROOT / "cluster" / "hf_batio3_vcneb_qe_preflight.slurm").read_text(encoding="utf-8")
    for fragment in (
        "#SBATCH --ntasks=1",
        "--n-images 7",
        "--image-workers 0",
        "--steps 0",
        "--validate-only",
        "QE_PP_MANIFEST",
        "--pp-manifest \"${pp_manifest}\"",
        "--no-climb",
        "--cell-interpolation log_strain --mapping auto",
    ):
        assert fragment in text
    assert "srun" not in text


def test_bto_qe_static_baseline_is_one_fixed_endpoint_32_mpi_scf() -> None:
    text = (ROOT / "cluster" / "hf_batio3_qe_static_baseline.slurm").read_text(encoding="utf-8")
    for fragment in (
        "#SBATCH --ntasks=32",
        "if [[ \"${RUN_DFT:-0}\" != 1 ]]",
        "QE_PP_MANIFEST",
        "--pp-manifest \"${pp_manifest}\"",
        "--ecutwfc \"${ecutwfc}\" --ecutrho \"${ecutrho}\"",
        "--static-only",
        "--no-climb",
        "ecutrho=${ECUTRHO:-600}",
        "srun --exclusive --nodes=1 --ntasks=32 pw.x",
    ):
        assert fragment in text

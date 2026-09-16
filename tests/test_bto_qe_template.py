"""Guard the no-DFT QE BTO distributed-validation template."""

from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_bto_qe_template_preserves_the_five_by_32_worker_contract() -> None:
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
        "--validate-only",
        "if [[ \"${RUN_DFT:-0}\" == 1 ]]",
        "--cell-interpolation log_strain --mapping auto --align-translation",
    ):
        assert fragment in text

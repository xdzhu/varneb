"""Regression guards for the licensed VASP BTO distributed submission gate."""

from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_bto_vasp_template_requires_explicit_licensed_inputs_and_worker_layout() -> None:
    text = (ROOT / "cluster" / "hf_batio3_vcneb_vasp_distributed.slurm").read_text(encoding="utf-8")
    for fragment in (
        "#SBATCH --ntasks=160",
        "RUN_DFT:-0",
        "VASP_BIN",
        "VASP_INITIAL_DIR",
        "VASP_FINAL_DIR",
        "ENDPOINT_IDENTITY_GATE",
        "scripts/validate_endpoint_identity_gate.py",
        "POTCAR",
        "five 32-MPI interior workers",
        "N_IMAGES:-7",
        "FMAX:-0.10",
        "VCNEB_GIT_REVISION",
        "VASP_COMMAND=\"srun --exclusive --nodes=1 --ntasks=${image_mpi}",
        "--image-workers \"${image_workers}\"",
        "--no-climb --mic --cell-interpolation log_strain --mapping auto --align-translation",
    ):
        assert fragment in text

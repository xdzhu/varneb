"""Static contract for the non-executable QE BTO pseudopotential candidate."""

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_sssp_bto_candidate_preserves_review_and_cutoff_boundaries() -> None:
    path = ROOT / "examples" / "bto_qe_sssp_1_3_pbe_precision_candidate.json"
    payload = json.loads(path.read_text(encoding="utf-8"))
    assert payload["approval_status"] == "candidate_requires_user_approval"
    assert payload["proposed_bto_run_settings"] == {
        "ecutwfc_Ry": 100,
        "ecutrho_Ry": 600,
        "rationale": "100 Ry preserves the project BTO wavefunction-cutoff target; 600 Ry is the largest SSSP precision density-cutoff recommendation among Ba, Ti and O.",
    }
    assert payload["species"]["Ba"]["md5"] == "ecceda5fc736cf81896add41b7758c6c"
    assert payload["species"]["Ti"]["md5"] == "88a00a6731bd790ddea75d31a80cb452"
    assert payload["species"]["O"]["md5"] == "0234752ac141de4415c5fc33072bef88"

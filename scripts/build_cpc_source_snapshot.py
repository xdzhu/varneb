"""Build an auditable CPC source-snapshot candidate from one Git commit.

This is deliberately not a publisher submission or a substitute for a manual
licence review. The archive contains only committed files, never worktree
changes, and retains the paper's examples, tests, figure sources, and audits.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile
import zipfile


ROOT = Path(__file__).resolve().parents[1]
FIGURE_PREFIXES = (
    "bto_qy0_restricted_sheet_27_20260929",
    "bto_frozen_soft_mode_289_20260929",
    "gan_600eV_local_joint_dft289_20260930_v2",
    "gan_600eV_atomic_transverse_162_20260929_v2",
)
REQUIRED_PATHS = frozenset({
    "LICENSE", "README.md", "MANIFEST.in", "pyproject.toml",
    "docs/USER_MANUAL.md", "vcneb/__init__.py",
    "examples/run_toy_vcneb.py", "scripts/rebuild_cpc_landscapes.py",
    "scripts/build_cpc_source_snapshot.py", "tests/test_cpc_numeric_claims.py",
    "tests/test_cpc_source_snapshot.py",
    "paper/VARNEB_CPC/varneb_CPC.tex", "paper/VARNEB_CPC/varneb.bib",
    "paper/VARNEB_CPC/MANUSCRIPT_EVIDENCE.md",
    "paper/VARNEB_CPC/SUBMISSION_READINESS.md",
    "paper/VARNEB_CPC/CPC_SOURCE_SNAPSHOT.md",
    "paper/VARNEB_CPC/evidence/gan_45p7_multibackend_vcneb_20260924.json",
    "paper/VARNEB_CPC/evidence/gan_45p7_final_chains_20260927/gan_vasp_45p7_final_chain.traj",
    "outputs/batio3_t_to_c_pbe100_dzp10au/bto_validation_provenance.json",
    "benchmarks/numerical_integrity/bto_qy0_27_node_sheet_20260929.json",
    "benchmarks/numerical_integrity/bto_transverse_soft_frozen289_20260929.json",
    "benchmarks/numerical_integrity/gan_600eV_ts_2d_17x17_refinement_20260930.json",
    "benchmarks/numerical_integrity/gan_600eV_atomic_tube_q9_20260929.json",
} | {
    f"paper/VARNEB_CPC/figures/{prefix}{suffix}"
    for prefix in FIGURE_PREFIXES
    for suffix in (".pdf", ".png", "_source_data.csv", "_qa.json")
})
RESTRICTED_NAMES = frozenset({
    ".env", "id_rsa", "id_ed25519", "wavecar", "chgcar",
    "aeccar0", "aeccar1", "aeccar2",
})
RESTRICTED_SUFFIXES = frozenset({
    ".upf", ".psp", ".psp8", ".vps", ".psml", ".orb",
    ".pem", ".key", ".p12",
})


def _git(repo: Path, *args: str) -> bytes:
    return subprocess.run(
        ["git", "-C", str(repo), *args], check=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    ).stdout


def committed_paths(repo: Path, commit: str) -> set[str]:
    names = _git(repo, "ls-tree", "-r", "-z", "--name-only", commit)
    return {name.decode("utf-8") for name in names.split(b"\0") if name}


def check_paths(paths: set[str], required: frozenset[str] = REQUIRED_PATHS) -> None:
    missing = sorted(required - paths)
    if missing:
        raise ValueError("CPC source snapshot is missing required committed files: "
                         + ", ".join(missing))
    restricted = []
    for path in paths:
        name = path.rsplit("/", 1)[-1].lower()
        bare_name = name
        for compression in (".gz", ".xz", ".bz2", ".zst"):
            if bare_name.endswith(compression):
                bare_name = bare_name[: -len(compression)]
                break
        if (bare_name == "potcar" or bare_name.startswith("potcar.")
                or bare_name in RESTRICTED_NAMES
                or Path(bare_name).suffix in RESTRICTED_SUFFIXES):
            restricted.append(path)
    if restricted:
        raise ValueError("potential, orbital, or secret-like files require removal/review: "
                         + ", ".join(sorted(restricted)))


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_snapshot(repo: Path, revision: str, output: Path,
                   required: frozenset[str] = REQUIRED_PATHS) -> dict:
    """Archive the pinned commit and verify ZIP entries against its Git tree."""
    if output.suffix.lower() != ".zip":
        raise ValueError("snapshot output must have a .zip extension")
    manifest = output.with_name(output.stem + ".manifest.json")
    if output.exists() or manifest.exists():
        raise FileExistsError("snapshot or manifest exists; choose a fresh output path")
    commit = _git(repo, "rev-parse", "--verify", "--end-of-options",
                  f"{revision}^{{commit}}").decode("ascii").strip()
    paths = committed_paths(repo, commit)
    check_paths(paths, required)
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(prefix=".cpc-source-", suffix=".zip",
                                     dir=output.parent, delete=False) as handle:
        temporary = Path(handle.name)
    try:
        _git(repo, "archive", "--format=zip", "--output", str(temporary), commit)
        with zipfile.ZipFile(temporary) as archive:
            entries = {item.filename for item in archive.infolist() if not item.is_dir()}
            if entries != paths:
                raise ValueError("Git ZIP contents differ from the pinned commit tree")
            if archive.testzip() is not None:
                raise ValueError("Git ZIP has a corrupt member")
        report = {
            "status": "candidate_for_manual_CPC_source_review_not_submitted",
            "commit": commit,
            "file_count": len(paths),
            "path_inventory_sha256": hashlib.sha256(
                ("\n".join(sorted(paths)) + "\n").encode("utf-8")
            ).hexdigest(),
            "archive_sha256": _sha256(temporary),
            "archive_bytes": temporary.stat().st_size,
            "required_files_verified": len(required),
            "restricted_filename_gate": "passed; manual licence/secret review still required",
            "source_scope": "all files in the pinned Git commit; uncommitted files excluded",
        }
        temporary.rename(output)
        manifest.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        return report
    finally:
        temporary.unlink(missing_ok=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=ROOT)
    parser.add_argument("--revision", default="HEAD")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(build_snapshot(args.repo, args.revision, args.output), indent=2))


if __name__ == "__main__":
    main()

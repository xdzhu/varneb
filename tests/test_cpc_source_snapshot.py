"""The CPC candidate must be pinned, complete, and free of known restricted files."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import zipfile

import pytest

from scripts.build_cpc_source_snapshot import build_snapshot, check_paths


def _git(repo: Path, *args: str) -> str:
    return subprocess.run(["git", "-C", str(repo), *args], check=True,
                          capture_output=True, text=True).stdout.strip()


@pytest.fixture
def repository(tmp_path: Path) -> Path:
    if not shutil.which("git"):
        pytest.skip("Git is required only to build CPC source snapshots")
    repo = tmp_path / "repo"
    repo.mkdir()
    _git(repo, "init", "-q")
    (repo / "README.md").write_text("committed source\n", encoding="utf-8")
    _git(repo, "add", "README.md")
    _git(repo, "-c", "user.name=VARNEB test", "-c",
         "user.email=test@example.invalid", "commit", "-qm", "source")
    return repo


def test_archive_uses_commit_not_dirty_or_untracked_files(
    repository: Path, tmp_path: Path,
) -> None:
    (repository / "README.md").write_text("uncommitted edit\n", encoding="utf-8")
    (repository / "untracked.txt").write_text("not for archive\n", encoding="utf-8")
    output = tmp_path / "varneb-cpc.zip"
    report = build_snapshot(repository, "HEAD", output,
                            required=frozenset({"README.md"}))
    committed_readme = subprocess.run(
        ["git", "-C", str(repository), "show", "HEAD:README.md"],
        check=True, capture_output=True,
    ).stdout
    with zipfile.ZipFile(output) as archive:
        assert archive.namelist() == ["README.md"]
        # Git archive may apply platform text-export filters; the dirty text
        # must still be absent and the committed logical content preserved.
        assert archive.read("README.md").replace(b"\r\n", b"\n") == committed_readme
    assert report["commit"] == _git(repository, "rev-parse", "HEAD")
    assert report["file_count"] == report["required_files_verified"] == 1
    assert report["archive_sha256"] == hashlib.sha256(output.read_bytes()).hexdigest()
    manifest = json.loads((tmp_path / "varneb-cpc.manifest.json").read_text())
    assert manifest == report
    before = output.read_bytes()
    with pytest.raises(FileExistsError, match="fresh output"):
        build_snapshot(repository, "HEAD", output, required=frozenset({"README.md"}))
    assert output.read_bytes() == before


@pytest.mark.parametrize("restricted", [
    "POTCAR", "Ga.UPF", "Ga.UPF.gz", "Ga.psp8", "Ga.orb", ".env", "secret.pem",
])
def test_known_restricted_files_are_rejected(restricted: str) -> None:
    with pytest.raises(ValueError, match="require removal/review"):
        check_paths({"README.md", restricted}, frozenset({"README.md"}))


def test_missing_material_provenance_fails_before_writing(
    repository: Path, tmp_path: Path,
) -> None:
    output = tmp_path / "incomplete.zip"
    with pytest.raises(ValueError, match="missing required committed files"):
        build_snapshot(repository, "HEAD", output,
                       required=frozenset({"README.md", "paper/missing.json"}))
    assert not output.exists()
    assert not (tmp_path / "incomplete.manifest.json").exists()


def test_tracked_potential_fails_before_writing(repository: Path, tmp_path: Path) -> None:
    (repository / "POTCAR").write_text("licensed potential\n", encoding="utf-8")
    _git(repository, "add", "POTCAR")
    _git(repository, "-c", "user.name=VARNEB test", "-c",
         "user.email=test@example.invalid", "commit", "-qm", "restricted")
    output = tmp_path / "restricted.zip"
    with pytest.raises(ValueError, match="require removal/review"):
        build_snapshot(repository, "HEAD", output,
                       required=frozenset({"README.md"}))
    assert not output.exists()

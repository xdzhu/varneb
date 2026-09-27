"""Historical curvature audits must not mutate an append-only DFT cache."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.create_audited_cache_view import create_view


def _audit(path: Path, name: str) -> None:
    path.write_text(json.dumps({
        "n_individually_audited_DFT_points": 1,
        "evaluations": [{"directory": name}],
    }), encoding="utf-8")


def test_view_links_only_the_historical_audited_points(tmp_path: Path) -> None:
    workdir = tmp_path / "append_only_cache"
    original = workdir / "eval-000001-abcdef123456"
    later = workdir / "eval-000002-abcdef123456"
    for directory in (original, later):
        directory.mkdir(parents=True)
        (directory / "result.json").write_text("{}", encoding="utf-8")
    audit = tmp_path / "audit.json"
    _audit(audit, original.name)
    view = tmp_path / "historical_view"
    try:
        assert create_view(workdir, audit, view) == 1
    except OSError as exc:
        pytest.skip(f"directory symlinks unavailable in this environment: {exc}")
    assert (view / original.name).is_symlink()
    assert (view / original.name / "result.json").read_text(encoding="utf-8") == "{}"
    assert not (view / later.name).exists()
    assert (later / "result.json").is_file()  # original cache remains intact
    with pytest.raises(ValueError, match="view must not exist"):
        create_view(workdir, audit, view)


def test_view_rejects_unsafe_or_missing_names_before_creating_directory(tmp_path: Path) -> None:
    workdir = tmp_path / "cache"
    workdir.mkdir()
    audit = tmp_path / "audit.json"
    _audit(audit, "../outside")
    view = tmp_path / "view"
    with pytest.raises(ValueError, match="unsafe"):
        create_view(workdir, audit, view)
    assert not view.exists()

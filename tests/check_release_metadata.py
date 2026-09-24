"""Guard the public package metadata and deliberate PyPI publish policy."""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXPECTED_DESCRIPTION = (
    "VARiable-cell Nudged Elastic Band code with universal first-frinciples calculators"
)


def _project_field(text: str, field: str) -> str:
    match = re.search(rf"^{re.escape(field)}\s*=\s*\"([^\"]+)\"\s*$", text, re.MULTILINE)
    if match is None:
        raise AssertionError(f"missing [project] field: {field}")
    return match.group(1)


def main(argv: list[str] = ()) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tag", help="Release tag that must match package version")
    args = parser.parse_args(argv)
    pyproject = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    assert _project_field(pyproject, "name") == "varneb"
    version = _project_field(pyproject, "version")
    assert re.fullmatch(r"\d+\.\d+\.\d+", version)
    version_source = (ROOT / "vcneb" / "version.py").read_text(encoding="utf-8")
    assert f'__version__ = "{version}"' in version_source
    if args.tag is not None:
        assert args.tag == f"v{version}", f"tag {args.tag} does not match v{version}"
    assert _project_field(pyproject, "description") == EXPECTED_DESCRIPTION
    assert _project_field(pyproject, "license") == "GPL-3.0-or-later"

    workflow = (ROOT / ".github" / "workflows" / "publish-pypi.yml").read_text(
        encoding="utf-8"
    )
    assert 'tags: ["v*"]' in workflow
    assert "  release:" not in workflow
    assert "workflow_dispatch:" in workflow
    assert "inputs.publish == true" in workflow
    assert "--tag \"${GITHUB_REF_NAME}\"" in workflow
    assert "on:" in workflow
    assert "push:\n    branches:" not in workflow

    print("release_metadata_regression=ok")


if __name__ == "__main__":
    main(sys.argv[1:])

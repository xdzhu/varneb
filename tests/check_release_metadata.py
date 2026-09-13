"""Guard the public package metadata and deliberate PyPI publish policy."""

from __future__ import annotations

import re
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


def main() -> None:
    pyproject = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    assert _project_field(pyproject, "name") == "varneb"
    assert re.fullmatch(r"\d+\.\d+\.\d+", _project_field(pyproject, "version"))
    assert _project_field(pyproject, "description") == EXPECTED_DESCRIPTION

    workflow = (ROOT / ".github" / "workflows" / "publish-pypi.yml").read_text(
        encoding="utf-8"
    )
    assert 'tags: ["v*"]' in workflow
    assert "release:" in workflow and "types: [published]" in workflow
    assert "workflow_dispatch:" in workflow
    assert "inputs.publish == true" in workflow
    assert "on:" in workflow
    assert "push:\n    branches:" not in workflow

    print("release_metadata_regression=ok")


if __name__ == "__main__":
    main()

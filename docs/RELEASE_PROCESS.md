# VARNEB release process

The package is published by GitHub Trusted Publishing. An ordinary push to
`main` does **not** publish to PyPI. Only a pushed `vX.Y.Z` tag or an explicit
workflow dispatch with `publish=true` starts the publish job. Creating a
GitHub Release from an existing tag does not start a second publish job.

Before tagging:

1. Confirm `pyproject.toml`, `vcneb/version.py`, and the intended tag agree.
2. Run `python tests/check_release_metadata.py --tag vX.Y.Z`,
   `python -m pytest -q`, and the manuscript/figure contract checks.
3. Build an sdist and wheel with `python -m build`; inspect their contents and
   install the wheel in an isolated environment. Build files are local and
   ignored by Git.
4. Review `git diff` and push the reviewed commit to `main`.
5. Push the annotated tag. Verify the publish workflow succeeds and the new
   version appears on PyPI before announcing the release.

The PyPI Trusted Publisher is bound to `xdzhu/varneb` and
`.github/workflows/publish-pypi.yml`. The publishing job uses the `pypi`
GitHub environment and OIDC `id-token: write`; no PyPI API token belongs in
this repository. PyPI versions are immutable, so fix a broken release with a
new version instead of overwriting one.

# Local VARNEB regression snapshot — 2026-09-28

This is a local Windows/PowerShell test run, **not** GitHub CI or a material
DFT validation. The code baseline was Git commit
`4e1b053f628e07c93dc0166f43f701c41595b297`; unrelated research files
were present in the working tree and were not cleaned or staged.

| Check | Command | Result |
| --- | --- | --- |
| All tests visible under `tests/`, including untracked research tests | `python -m pytest -q tests` | 568 passed, 1 skipped, 209 warnings; 84.54 s |
| Git-tracked `tests/test_*.py` only | `$testFiles = @(git ls-files 'tests/test_*.py'); python -m pytest -q -- $testFiles` | 549 passed, 1 skipped, 209 warnings; 40.25 s |

Environment: Python 3.10.9, ASE 3.28.0, NumPy 1.23.5, SciPy 1.10.1,
pytest 7.1.2. Both runs exited with status zero. The 19 additional passing
tests in the first run came from untracked research files, so only the second
line restricts collection to tracked test files. Because both runs used a
dirty worktree, neither is proof of a pristine checkout. The tests exercise
parsers, optimizer logic, backend contracts, figure evidence, and analysis
code; they do **not**
replay every remote SCF or certify any TS or conditional PES.

Warnings comprise an ASE optimizer notice about passing precomputed forces to
`step()` and upstream spglib/Phonopy deprecations. The precomputed-force call
currently avoids an extra expensive image evaluation; changing it requires a
separate evaluation-count regression rather than warning suppression. No CI,
release, PyPI publish, or new cluster calculation was started for this check.

## Local distribution smoke

After the explicit `MANIFEST.in` and README-link fix, an isolated
`python -m build --outdir tmp/package_qa_20260928b` produced the 0.0.2 sdist
and wheel. `python -m twine check` passed both. Archive inspection found the
manual and analytic toy script in the sdist, no `cluster/`, `outputs/`,
`validation/`, `benchmarks/`, or `paper/` research trees, and only the
`vcneb` package plus metadata in the wheel. A temporary venv installed the
wheel with `--no-deps` and `--system-site-packages`; import resolved to that
venv's `site-packages`, `varneb --version` reported `0.0.2`, and
`varneb backends --json` returned the declared capability matrix. The toy
script returned a `0.250004 eV` barrier both via the installed wheel and
from the unpacked sdist. `python tests/check_release_metadata.py --tag
v0.0.2` passed. These checks do not constitute a clean dependency install,
release tag, or PyPI upload.

An initial `--no-isolation` build failed because this machine has setuptools
65.6.3, below the project's declared `>=68` minimum; the isolated build
resolved that environment mismatch without changing package metadata or
publishing anything. The local final archive SHA-256 values were
`5631230e9b4ec2935590df0c115e601f2d2eb2606858dd46753dd7331a340646`
(sdist) and
`71b8523f53a5e6d066052d2bf34cb84af0a6a40c0d84a754f10f2c0069750a17`
(wheel). Both artifacts remain under ignored `tmp/` and are not release
artifacts.

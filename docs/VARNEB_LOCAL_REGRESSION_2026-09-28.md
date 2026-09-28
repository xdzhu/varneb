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

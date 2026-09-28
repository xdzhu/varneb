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

## Clean-dependency install and current tracked tests

A later isolated `python -m build --outdir tmp/package_smoke_20260928_01`
at source commit `f7e81cf` produced another 0.0.2 sdist/wheel. The subsequent
`3ccddcc` figure commit changed no packaged code, metadata, README, manual, or
toy-example source. The new sdist SHA-256 is
`ac2ec9bcbfd189f28667f9fc943e720e879f6cebac24572896f26b30015803a7`;
the wheel SHA-256 is
`abbddaddee5e64d2a36f8edc6dd09f72cafd9c7f24770dd6e81ab10900b86b7c6`.
The sdist has 46 files, including the required README, manual and toy example,
with no research-output, benchmark, cluster, test, or paper tree. The wheel
has 37 entries consisting only of `vcneb/` and distribution metadata.
`python -m twine check` passed both.

Unlike the earlier system-site-packages smoke, a fresh Python 3.10.9 venv
without inherited packages installed this wheel and resolved ASE 3.29.0,
NumPy 2.2.6, SciPy 1.15.3, and Matplotlib 3.10.9. From a directory outside
the source checkout, `vcneb.__file__` resolved inside that venv's
`site-packages`; `varneb --version` returned `0.0.2`, and both
`varneb backends --json` and `varneb optimizers --json` completed. The
unpacked sdist's analytic toy example returned `barrier_eV=0.250004` and
`delta_eV=0.000000`, and wrote a final chain and band plot. At `3ccddcc`,
`python tests/check_release_metadata.py --tag v0.0.2` passed, and the
Git-tracked `tests/test_*.py` set gave **562 passed, 1 skipped** in 36.74 s
(209 warnings). The tested working tree still contains unrelated user
research files; this is not a pristine-checkout claim. No release tag, CI,
PyPI upload, DFT calculation, or other cluster job was started.

## Config-input safety after the distribution smoke

The isolated 0.0.2 wheel also completed `varneb init` → `validate-config` →
`prepare` from a separate run directory using the toy run's two endpoints and
the generic ASE backend. The preflight reported `calculator_attached=false`
and wrote seven images. A second identical `prepare` was idempotent;
`varneb run` without `--execute` returned status
2 before any calculator could start.

That exercise exposed a source-level risk not yet corrected in the tested
0.0.2 wheel: JSON `"false"` had been coerced to Python `True`, a fractional
`n_images` to an integer, and unknown keys (for example `pressure_GPa`) were
ignored. The development source now rejects these with named-field errors at
`validate-config`, before endpoint preparation or DFT. It also rejects
unknown calculator keys and non-object JSON. The focused config/CLI tests gave
24 passes; the Git-tracked suite after this change gave **575 passed, 1
skipped** (209 warnings, 28.91 s). This is a source-tree safety improvement
for a subsequent release, not a retroactive claim about the tested wheel.

## Current-source clean-wheel and sdist smoke

The post-config-fix development source at base commit `3661880` plus the
selective README, manual, backend-note, and test changes was rebuilt with
isolated `python -m build` under ignored `tmp/package_smoke_20260928_03`.
`python -m twine check` passed both archives, and
`python tests/check_release_metadata.py --tag v0.0.2` passed. The 51-entry
sdist includes the README, detailed manual, and analytic toy example; it
contains no `cluster/`, `outputs/`, `validation/`, `benchmarks/`, `paper/`, or
`tests/` tree. The 37-entry wheel contains only `vcneb/` and its metadata.
The local test-archive SHA-256 values are
`293cc04668b3dcb7d22b87ca961fd602e11400a707bfd01cb26e341af74fd844`
(sdist) and
`0cc028912fd963597c881966ed5cab1d7ea23843e73003b5c6e6958bbc627fee`
(wheel); neither is a tagged or published release artifact.

A fresh Python 3.10.9 virtual environment was verified to have
`include-system-site-packages = false`. The first dependency resolution
hit transient PyPI timeouts for `python-dateutil`; an independent HTTP HEAD
check later returned 200. Pin-installing `python-dateutil==2.9.0.post0`
into that *same clean venv* allowed the declared `wheel[plot,mode]`
dependencies to resolve from PyPI (including ASE 3.29.0, NumPy 2.2.6,
SciPy 1.15.3, and Matplotlib 3.10.9). This was a network/retry workaround,
not a change to package requirements or a system-site-packages fallback.
After installing the final rebuilt wheel, `pip check` passed and, from a
separate run directory, `vcneb.__file__` resolved to that venv's
`site-packages`; `varneb --version` reported `0.0.2`.

The unpacked *final* sdist's analytic toy run again returned
`barrier_eV=0.250004` and `delta_eV=0.000000`. On the final installed wheel,
`varneb init` → `validate-config` → `prepare` succeeded for the toy endpoints
and generic ASE backend. The seven-image calculator-free preparation was
idempotent (trajectory SHA-256 unchanged on the second call). Three invalid
configuration probes—string `"false"` for `mic`, fractional `n_images=7.5`,
and an unknown `pressure_GPa` key—each failed with exit code 2; so did
`varneb run` without `--execute`. No DFT was launched.

The clean ASE 3.29.0 Windows installation did not provide
`ase.calculators.abacus`, and `varneb doctor --backend abacus` correctly
reported that the optional adapter or an explicit custom factory is needed.
This is a local installation prerequisite, not a contradiction of the
archived material-level ABACUS validation. The static `validated` label and
the README/manual now distinguish those two meanings.

Post-change local tests: all visible `tests/`, including unrelated untracked
research tests, yielded **599 passed, 1 skipped, 209 warnings** in 32.86 s;
the Git-tracked `tests/test_*.py` set yielded **580 passed, 1 skipped, 209
warnings** in 34.50 s. Both runs used a dirty working tree, so they do not
prove a pristine checkout, remote material reproducibility, or a new release.
No CI, release tag, PyPI upload, DFT run, or cluster job was started by this
**packaging smoke**. Subsequent separately audited BTO/GaN static-grid jobs
are recorded in `paper/VARNEB_CPC/STAGE_GATE_2026-09-27.md`.

## CPC figure-evidence and manuscript-gate regression

After the 81-point BTO and GaN figures entered the paper, a read-only test
now checks each committed figure source CSV against its SHA-256 figure-QA pin,
the independent 22/56-new-point raw audit, and a complete unique 9×9
coordinate grid. It also checks that the corresponding PDF is included in
the intended main or supplementary TeX source. The older manual manuscript
claim guard had gone stale because it matched line wrapping and an earlier
figure location; it now normalizes prose whitespace, checks the current
scientific limits, and runs under pytest. The focused checks give **3 passed**.

The full visible local suite after these changes gave **609 passed, 1 skipped,
209 warnings** in 74.94 s. This run did not launch any calculator, CI, or
release workflow. The previous isolated wheel/sdist smoke above is still
useful but predates the latest README/manual text and these repository-only
paper tests; it is not evidence that an exact-current-revision package has
been rebuilt or tested in a pristine checkout.

## Committed-package-source smoke after the optimizer JSON fix

At commit `06baea6`, `git archive` exported only the tracked packaging inputs
(`pyproject.toml`, `MANIFEST.in`, licence, both READMEs, `vcneb/`, the detailed
manual, and the analytic toy example) into a separate ignored directory.
An isolated `python -m build` produced both `varneb-0.0.2` archives from that
clean source snapshot; `twine check` passed both. The sdist contains the
README, manual, toy and package source, not `tests/`, `cluster/`, `outputs/`,
`benchmarks/` or `paper/`. Local artifact SHA-256 values are:

- sdist: `b7fe06534a692b77eb4d9f2239eb567934e86e1a3bb84518052d0bdcac3fedca`;
- wheel: `6a67c181630246bea7fb0b4eb7c6017254b3b6fd04b87cbfeb309ce980f5bc33`.

The clean Python 3.10.9 venv created for the previous committed snapshot has
`include-system-site-packages = false`. Its dependencies were installed from
PyPI, then the `06baea6` wheel replaced the prior local wheel using
`--force-reinstall --no-deps`. From a separate run directory, `vcneb.__file__`
resolved to that venv's `site-packages`; `pip check`, `varneb --version`,
`backends --json`, `optimizers --json`, `varneb init`, and `validate-config`
passed. The optimizer registry now emits eight JSON boolean `true` values
rather than string `"true"`, consistent with the backend matrix. A copied
analytic toy example ran against the installed wheel and returned
`barrier_eV=0.250004`, `delta_eV=0.000000`. The full local source suite after
the change gave **609 passed, 1 skipped, 209 warnings** in 56.69 s.

This is a clean *package-source* snapshot and installed-wheel smoke, not a
full-repository pristine checkout, a new tagged release, a PyPI publication,
or material-level DFT reproducibility. The version remains `0.0.2`; these
local artifacts are not the already-published PyPI files. No CI or DFT was
started by this check.

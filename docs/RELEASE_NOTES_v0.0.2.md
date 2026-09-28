# VARNEB 0.0.2

This release packages the current calculator-independent VCNEB core and
documents the completed five-backend GaN B4→B1 validation at 45.7 GPa.
The CP2K path converged at step 52 with 29 total images and a maximum
generalized force of 0.09409 eV/Å; its 0.29276 eV/GaN barrier and source
data are included in the CPC manuscript package. The manuscript now contains
the converged cross-backend figure, BTO mode analysis, and explicit evidence
and limitation statements.

The case library gains a compact CP2K GaN fixture and a production caution:
multiple independent CP2K MPI jobs must have verified, disjoint CPU affinity.
The Slurm template defaults to one CP2K worker and rejects unverified
multi-worker launches. This does not change the physics or silently alter an
existing run; users must recheck endpoint and static-force/stress contracts
on their own installation.

The release workflow publishes only on an explicit version tag or manual
dispatch. Ordinary source pushes and GitHub Release creation do not publish.

## Local artifact smoke, 2026-09-28

From clean commit `c88f5d3`, an isolated `python -m build` produced the
0.0.2 sdist and wheel; `python -m twine check` passed for both. The sdist
contains the README, detailed manual, and analytic toy example, while the
wheel contains the import package and console entry points rather than raw
research outputs. A fresh Windows virtual environment without
`--system-site-packages` resolved ASE 3.29.0 and NumPy 2.2.6, installed the
wheel, and ran the installed `varneb --version`, `backends --json`, and
`optimizers --json` commands. The unpacked sdist toy returned a 0.250004-eV
barrier. With two generated Ar endpoints, installed `varneb init`,
`validate-config`, and `prepare` wrote a seven-image, two-atom-per-image
calculator-free chain; a repeated `prepare` preserved its SHA-256. Bare
`varneb run` refused execution without `--execute`. This is local artifact
validation, not proof of the currently hosted PyPI file or any DFT backend,
and it did not trigger publication.

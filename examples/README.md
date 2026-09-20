# VARNEB examples and validation cases

The repository separates calculator-free demonstrations from licensed
first-principles workflows.  No VASP binary or PAW dataset is distributed.
For every VCNEB driver, `n_images` counts both fixed endpoints; only the
`n_images - 2` interior images are dispatched to workers.

The stable public entry points are [`../README.md`](../README.md),
[`quickstart/`](quickstart/), and [`cases/`](cases/).  The parent-level
scripts below remain compatibility entry points for the existing validation
suite.

## Calculator-free examples

The top-level Python scripts demonstrate the analytic model, fixed-cell ASE
reduction, convergence matrix, mode guidance, directional constraints, restart,
and release-and-refine interfaces.  These are suitable for local testing and do
not require a DFT executable.

## Material cases

| Case | Purpose | Status and documentation |
| --- | --- | --- |
| `batio3_vasp_pbe_paw_static/` | Minimal VASP static-input fixture for backend validation | Input-only example; licensed POTCAR is not committed |
| `gan_qian_2013/` | GaN B4/B3 to B1 literature-path studies at finite pressure | Includes auditable seeds and accepted B3/B1 endpoints; see `docs/GAN_QIAN_2013_REPLICATION_PLAN.md` |
| `cdse_sheppard_2012/` | CdSe rock-salt to wurtzite cell- and atom-mapped G-SSNEB comparison | Preparation/launch description only while production paths remain under review; see `docs/CDSE_SHEPPARD_2012_REPLICATION_PLAN.md` |
| `prl2023_hfo2_fig2a/` | HfO2 polymorph paths from the 2023 PRL figure-2a study | Reproduction contract and failure-safe VASP workflow; not yet a manuscript result |

LAMMPS and CP2K adapters use the same isolated-image factory as the DFT
backends.  Their adapter status is intentionally not a material-validation
claim; use `varneb doctor` and the HF module list before selecting a potential
or basis/pseudopotential set.

Material directories contain compact structures, preparation metadata, or
documentation.  Scheduler templates live in `cluster/`; run outputs and restart
files remain outside the distributable example package.

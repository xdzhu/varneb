# Cluster launch templates

The Slurm templates in this directory are for the Hefei cluster.  The ordinary
templates run ABACUS through `srun` while the Python driver evaluates images
sequentially; the opt-in `hf_batio3_vcneb_parallel.slurm` template uses
exclusive concurrent image job steps.
The production BaTiO3 and HfO2 templates default to 32 MPI tasks (one CPU per
task); these are the project-wide settings for subsequent runs.  The
templates assume the shared `/public/home/iai806` layout and the
`hfacnormal01` partition.  The BaTiO3 template uses ABACUS Dojo-NC-FR at
`ecutwfc=100 Ry` with the `Orb-DZP-10au` basis directory; Ba, Ti and O all
use their 10 au DZP orbitals.

The HfO2 endpoint template uses ABACUS-native `calculation cell-relax` with
`relax_method bfgs`, `relax_nmax`, `force_thr_ev`, and `stress_thr`; ASE is
used only to convert the input endpoint to `STRU` and to export the final
`OUT.ABACUS/STRU_ION_D` as a VASP-readable `CONTCAR`.  If native step control
becomes unstable, the ASE-driven `examples/relax_abacus_endpoint.py` is an
equivalent fallback on the same energy surface and degrees of freedom; its
`--stress-kbar` option enforces the same endpoint gate.

The ready-to-submit conservative fallback is
`cluster/hf_hfo2_endpoint_ase.slurm` (32 MPI, 100 Ry, Orb-DZP-10au,
`BFGS`, default `maxstep=0.005`). Set `STRUCTURE` and `WORKDIR` when
submitting one endpoint; accept the endpoint only when both `force_converged`
and `stress_converged` in `relax_summary.json` are true.
Use `scripts/promote_hfo2_endpoints.py` to atomically publish the two accepted
`CONTCAR` files into `relaxed_T/` and `relaxed_PO/` before starting VCNEB.

From `ssh hf`, after checking `sinfo` and `squeue`, sync the repository and
submit the two endpoint relaxations independently:

```bash
sbatch --export=ALL,ENDPOINT_TAG=pbe100_dzp10au cluster/hf_batio3_endpoint.slurm cubic
sbatch --export=ALL,ENDPOINT_TAG=pbe100_dzp10au cluster/hf_batio3_endpoint.slurm tetragonal
```

After both `CONTCAR` files exist, the default T→C no-climb preconvergence run is:

```bash
sbatch --export=ALL,DIRECTION=tetragonal_to_cubic,N_IMAGES=7,STEPS=300,FMAX=0.10,MAXSTEP=0.02 \
  cluster/hf_batio3_vcneb.slurm
```

Set `DIRECTION=cubic_to_tetragonal` explicitly for the reverse calculation.

The BaTiO3 template uses symmetric-positive `log_strain` cell interpolation,
automatic same-species endpoint mapping, and minimum-image atom displacements
so that the initial path follows the tested variable-cell geometry.  It also
rejects a path below `1.6 Angstrom` minimum separation or above `0.10` cell
deformation before launching any image calculator.

The scripts use `${SLURM_NTASKS}` in the `srun` command, so the actual
allocation controls the MPI width.  The default is now 32 tasks:

```bash
sbatch --export=ALL,N_IMAGES=7,STEPS=300 \
cluster/hf_batio3_vcneb.slurm
```

For image-level concurrency, use `hf_batio3_vcneb_parallel.slurm`.  This is a
different resource model: the Python process remains the VCNEB controller and
launches concurrent `srun --exclusive` steps, one per image calculator.  The
default allocation is 128 tasks, allowing four 32-MPI image workers on one
128-CPU node:

```bash
sbatch --export=ALL,N_IMAGES=7,STEPS=300,IMAGE_WORKERS=4,IMAGE_MPI=32 \
  cluster/hf_batio3_vcneb_parallel.slurm
```

`IMAGE_WORKERS * IMAGE_MPI` must not exceed the allocation.  The 32-MPI value
is therefore the width of each ABACUS worker, not the total controller job
width.  The parallel template is opt-in; the ordinary template remains the
serial-image reference for reproducibility and debugging.
Each parallel run also writes an image-worker JSONL manifest under its work
directory, so failed workers and retries remain visible beyond Slurm stdout.

`hf_batio3_vcneb_qe_distributed.slurm` is the corresponding QE BTO T→C
validation gate. It is intentionally fixed to seven total images and five
32-MPI interior workers (160 ranks over two nodes), with `--no-climb` because
the established ABACUS reference is barrierless. By default it runs only the
driver's `--validate-only` preflight; set `RUN_DFT=1` only after supplying the
reviewed `QE_ENV_SCRIPT`, `ESPRESSO_PSEUDO`, and Ba/Ti/O UPF names. The QE
cutoffs are explicit QE/Ry parameters and must be independently converged;
they are not a mechanical conversion of VASP settings.

`hf_batio3_vcneb_vasp_distributed.slurm` provides the equivalent VASP BTO
execution gate: it is fixed at seven total images and five 32-MPI interior
workers (160 ranks over two nodes), keeps endpoints cached at the manager, and
uses ordinary no-climb NEB. It refuses to start unless `RUN_DFT=1`, a licensed
`VASP_BIN` visible on `hfacnormal01`, and `VASP_INITIAL_DIR`/`VASP_FINAL_DIR`
are provided. The initial directory must contain the frozen BTO `CONTCAR` plus
reviewed `INCAR`, `KPOINTS` and `POTCAR`; it is never inferred from ABACUS
inputs. Run `examples/run_vcneb_vasp.py --validate-only` on those directories
before requesting the 160-rank allocation.

`hf_batio3_gamma_phonon_abacus.slurm` seeds the separate first-principles
Gamma-mode analysis at the already audited cubic or tetragonal BTO endpoint.
It uses a `1×1×1` finite-displacement cell and static 32-MPI ABACUS SCFs. The
default `RUN_DFT=0` only generates and records displacement folders; `RUN_DFT=1`
is required to evaluate forces and assemble `FORCE_SETS`. This first workflow
is a Gamma-mode decomposition, not a converged phonon dispersion or QHA study.

For exact-state recovery across a new work directory, set
`IMAGE_CACHE_DIR=/path/to/image-cache` and optionally
`IMAGE_CACHE_NAMESPACE=ecut100-dzp10au-k2x2x2-scf1e-8`.  The cache is keyed by
image index plus exact species/positions/cell/PBC bytes, rejects a namespace
mismatch, and records hits/misses in the worker manifest.  Keep one cache
namespace per calculator parameter set; it is opt-in and does not replace the
fixed endpoint cache in the VCNEB controller.
For ASE's ABACUS adapter, keep atoms grouped by species in the input order;
the project patch suppresses the otherwise process-global `ase_sort.dat` for
this identity ordering, preventing concurrent parser races.  Use the serial
template for non-identity ordering unless the calculator supplies a
directory-local sort-file implementation.

For seven total images (five interior images) in one wave, use
`hf_hfo2_vcneb_distributed.slurm`.  It requests two 128-CPU nodes (160 worker
tasks),
keeps one Python VCNEB manager, and launches five isolated 32-MPI ABACUS
workers through `srun --exclusive --nodes=1`; the two fixed endpoints are
evaluated once and cached by the manager.  This removes the 4+3 worker waves
of the one-node template while preserving the same calculator settings:

```bash
sbatch --export=ALL,N_IMAGES=7,IMAGE_WORKERS=5,IMAGE_MPI=32,STEPS=300,\
RESUME=1,RESUME_TRAJECTORY=/path/to/vcneb.traj,RESUME_STEP=-1 \
  cluster/hf_hfo2_vcneb_distributed.slurm
```

`N_IMAGES` always includes both endpoints, so the number of worker images is
`N_IMAGES-2` (the template derives this automatically when `IMAGE_WORKERS` is
omitted).  The distributed template is the preferred production layout for
seven total images;
the one-node four-worker template remains useful for smoke tests and limited
allocations.  `--nodes=1` on each worker step is intentional: it prevents a
single 32-rank image calculation from being split across both nodes.

The production default is `CELL_INTERPOLATION=log_strain`; to perform an
independent initial-path check without changing the endpoints, set
`CELL_INTERPOLATION=linear` and use a distinct `WORKDIR`.  `MAPPING=auto` and
`ALIGN_TRANSLATION=1` remain the default audited endpoint gauge.
When the shared source directory is not a Git checkout, pass
`VCNEB_GIT_REVISION=<commit>` in `--export` so the preflight and summary retain
the exact source revision; otherwise they record `remote-sync-unknown`.

Do not request `--exclusive` for these small tests.  For the 12-atom HfO2
fixture, after endpoint relaxations complete, a staged no-climb preconvergence
run is:

```bash
sbatch --export=ALL,N_IMAGES=7,STEPS=300,FMAX=0.10,NO_CLIMB=1 \
  cluster/hf_hfo2_vcneb.slurm
```

Record the 32-rank allocation in each run manifest and compare wall time against
the earlier 16-rank endpoint jobs; do not alter an already-running allocation.

The HfO2 template defaults to a `0.25` deformation preflight threshold because
the independently relaxed 12-atom endpoints give a measured path deformation
of about `0.124`; override `MAXIMUM_DEFORMATION` downward when the physical
case justifies a stricter path gate.

An independent nine-image run can use a distinct work directory:

```bash
sbatch --export=ALL,N_IMAGES=9,STEPS=300,FMAX=0.10,MAXSTEP=0.02,\
WORKDIR=validation/batio3_cubic_to_tetragonal/vcneb_n9_fire_pbe100 \
  cluster/hf_batio3_vcneb.slurm
```

Use `NO_CLIMB=0` only after the no-climb path is stable.  Use `RESUME=1` with
the same `WORKDIR` after an interrupted run, or provide
`RESUME_TRAJECTORY=/path/to/previous/vcneb.traj` with a new `WORKDIR` for a
recoverability test or a conservative-step restart.  Every job writes its Slurm
stdout/stderr under the shared run root and the VCNEB driver writes the
per-image summary under its work directory.

The BTO image-count and CI acceptance gates are specified in
`docs/batio3_validation_protocol.md`.

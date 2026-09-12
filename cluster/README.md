# Cluster launch templates

The Slurm templates in this directory are for the Hefei cluster.  They use
one 16-task allocation by default and run ABACUS through `srun`, while the
Python driver evaluates the images sequentially.  The 128-CPU nodes are not
reserved exclusively for these small five-atom tests.  The templates assume the
shared `/public/home/iai806` layout and the `hfacnormal01` partition.

From `ssh hf`, after checking `sinfo` and `squeue`, sync the repository and
submit the two endpoint relaxations independently:

```bash
sbatch cluster/hf_batio3_endpoint.slurm cubic
sbatch cluster/hf_batio3_endpoint.slurm tetragonal
```

After both `CONTCAR` files exist, a no-climb preconvergence run is:

```bash
sbatch --export=ALL,N_IMAGES=7,STEPS=300,FMAX=0.03,MAXSTEP=0.02 \
  cluster/hf_batio3_vcneb.slurm
```

The BaTiO3 template uses symmetric-positive `log_strain` cell interpolation,
automatic same-species endpoint mapping, and minimum-image atom displacements
so that the initial path follows the tested variable-cell geometry.  It also
rejects a path below `1.6 Angstrom` minimum separation or above `0.10` cell
deformation before launching any image calculator.

The scripts use `${SLURM_NTASKS}` in the `srun` command, so the actual
allocation controls the MPI width.  For this small cell start with 16 tasks;
after a timing check, a larger run can request 32 tasks without changing the
workflow:

```bash
sbatch --ntasks=32 --export=ALL,N_IMAGES=7,STEPS=300 \
  cluster/hf_batio3_vcneb.slurm
```

Do not request `--exclusive` for these small tests.  For a larger HfO2 cell,
select the MPI width from an explicit scaling test and record it in the run
manifest.

An independent nine-image run can use a distinct work directory:

```bash
sbatch --export=ALL,N_IMAGES=9,STEPS=300,FMAX=0.03,MAXSTEP=0.02,\
WORKDIR=validation/batio3_cubic_to_tetragonal/vcneb_n9_fire_pbe100 \
  cluster/hf_batio3_vcneb.slurm
```

Use `NO_CLIMB=0` only after the no-climb path is stable.  Use `RESUME=1` with
the same `WORKDIR` after an interrupted run.  Every job writes its Slurm
stdout/stderr under the shared run root and the VCNEB driver writes the
per-image summary under its work directory.

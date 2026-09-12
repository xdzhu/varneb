# Cluster launch templates

The Slurm templates in this directory are for the Hefei cluster.  They use
one exclusive 40-task node per job and run ABACUS through `srun`, while the
Python driver evaluates the images sequentially.  The templates assume the
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

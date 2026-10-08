# Common-substrate HfO2 endpoint starters (G2, DFT pending)

Ten geometry starters cover five registered phases/ordered variants at the two
preregistered substrate strains, epsilon=0 and +0.01. The independent +0.005
holdout is intentionally absent. These are **not relaxed endpoints or energies**.

Reproduce in a fresh namespace:

```powershell
python -m scripts.prepare_hfo2_clamped_endpoints --output path/to/fresh/seeds
```

`clamped_seed_manifest.json` hashes all sources and per-seed manifests. The
transformation is a proper cyclic rotation: new Cartesian x/y/z=old z/x/y;
cell rows are also cyclically reordered. Atom order is never changed. The
source T long axis is therefore research x; substrate rows 0/1 are research
x/y. All phases use **the same T substrate**, not their own relaxed plane.
The third vector and all atomic positions subsequently relax, including two
out-of-plane tilts; no symmetry or mode is forced to retain a phase.

`free_phase_in_plane_length_changes` reports how far each starter is distorted
from its own free phase. A zero T-referenced substrate strain can impose a
several-percent strain on PO/M. This is a different ensemble from the P=0
free-cell calculation. No rotated/strained raw E/F/stress or Berry result is
cached from that calculation.

After G1 and a reviewed bounded canary, run a selected endpoint through:

```bash
python -m scripts.relax_clamped_ase_endpoint \
  --seed-manifest path/to/seeds/strain_0000/PO_plus/endpoint_seed.json \
  --workdir path/to/fresh/clamped-endpoint \
  --factory examples.hfo2_fixed_input_factory:make_factory \
  --parameters benchmarks/hfo2_channels/20261008/M_endpoint_factory_parameters.json \
  --command 'mpirun -np 32 /public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/app/bin/abacus' \
  --fmax 0.03 --stress-kbar 2 --maxstep 0.02 --steps 4
```

The factory byte-checks the unchanged INPUT/KPT/pseudopotential/full 10-au DZP
contract before every real SCF. `ecutwfc=100 Ry` remains fixed. This command
belongs in a 32-CPU HF/hfacnormal01 allocation, **not on a login node**; the
matching Slurm template is `cluster/hf_hfo2_clamped_endpoint_20261008.slurm`.
It does not submit itself or expand the G2 production matrix.

The BFGS gate uses atomic force <0.03 eV/A and open-traction norm <2 kbar.
For the tilt-released boundary the traction is `(sigma+P I) @ n`; substrate
reaction stresses may be large and are recorded separately. Geometric guards
run before a proposed structure can mutate the atoms or reach DFT.
`endpoint_relax_summary.json` distinguishes step-limit from numerical
convergence, with **phase/variant audit still pending** even after convergence.
Every raw SCF has a separate immutable directory. No unreviewed automatic
restart, phase-restoring constraint, CI, or physical-parameter adjustment occurs.

Current evidence: preparation/ASE-EMT tests only. No HfO2 clamped SCF has been
submitted and no strain-dependent barrier or phase-stability conclusion follows.

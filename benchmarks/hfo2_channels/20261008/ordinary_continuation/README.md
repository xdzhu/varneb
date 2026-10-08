# Terminal-run exact-cache continuation (2026-10-08)

The two R12 allocations finished normally at their ten-step health-test limit,
not at the ordinary 0.10 eV/A convergence threshold. Their latest complete
step10 chains become the seeds of bounded R14 ordinary continuations:

| Chain | Source allocation | Source residual (eV/A) | Total / interior | Continuation |
|---|---|---:|---:|---|
| T to PO, lifted sampling repair (`gap`) | 28274895, COMPLETED | 0.28499524288774647 | 10 / 8 | 28300425 |
| PO to M | 28275259, COMPLETED | 0.2933253052250083 | 9 / 7 | 28298794 |

`gap` is the historical sampling-repair label, not an electronic-gap claim.
The new allocations are ordinary same-input relaxations, not additional
scientific channels, saddle certification, or an acceleration benchmark.

## Preparation and evidence

`scripts/prepare_hfo2_observation_restart.py` refuses a live source allocation,
an older complete snapshot, an already converged chain, an existing output,
or changed input/geometry/raw SCF evidence. Every image is checked against
its original ordered periodic structure and freshly parsed energy, forces,
and stress. The seed has no stored calculator; its immutable factory
parameters point to the exact audited raw SCFs, including both endpoints.
No atom remapping, coordinate relifting, source overwrite, or DFT is performed.

For both cases the actual R12 production factory was independently exercised
on HF ASE3.23.1b1 while its DFT call method was replaced *in memory* with a
hard failure. All 19 seed evaluations were obtained from exact geometry
caches, and ordinary force replay reproduced the terminal residuals. Separate
public-CLI geometry/lift preflights passed. See the four `*_preflight.json`
reports, each restart's `manifest.json`, and the parent
`submission_handles_r14.json`. Local numeric evidence can be tested without
HF, ABACUS, or access to the licensed inputs:

```bash
python -m pytest -q tests/test_hfo2_observation_restart.py
```

This numerical replay does not replace the recorded live raw-SCF audits.
Raw SCF directories are on HF; only their hashes and E/F/stress records are
published here. Exact `seed.traj` and POSCAR bytes are versioned explicitly
so that Git newline conversion cannot silently invalidate the evidence.

## Unchanged physics and bounded execution

Both continuations use the immutable R12 production archive (`cc536a2`),
ABACUS/PBE/LCAO 100 Ry with complete 10-au DZP, original 2x2x2 Gamma mesh,
Hf4O8, P=0, E=0, ordinary NEB 0.10 eV/A, k=0.2, FIRE maxstep=0.02.
Only the execution segment grows to **at most 80 steps / 24 hours**.
FIRE momenta/timestep history are **not** restored: this is a fresh-optimizer
geometry/cache continuation and must not be described as state restoration
or evidence of faster convergence.

Each allocation requests one node, 32 CPUs and launches one `mpirun -np 32`
SCF at a time with one thread per rank; endpoints remain cached. Both jobs
depend on termination of *both* R13 switching-health allocations
(`afterany:28288045:28288063`) to retain at most two active study chains.
The dependency is a capacity gate, not proof of predecessor convergence.
No automatic resubmission, requeue, CI, physical-parameter adjustment, or
single-force-rebound cancellation is introduced.

Recent source steps took roughly 11--21 minutes. The 80-step cap can
therefore exceed the 24-hour wall cap; whichever is reached first ends the
segment. There is no defensible exact convergence time from these pilots.
Review a complete terminal segment before any subsequent continuation.

# HfO2 polymorph-path reproduction case

This case reproduces selected intrinsic paths associated with Fig. 2a of the
2023 HfO2 polymorphism study.  The detailed literature audit, unresolved input
choices, path-label corrections, and acceptance criteria are documented in:

- `docs/prl2023_fig2a_vcneb_reproduction_notes.md`
- `docs/PRL2023_HFO2_FIG2A_VARNEB_PLAN.md`
- `docs/vasp_input_contract.md`

The preparation entry point is `scripts/prepare_prl2023_hfo2_fig2a.py`.
Hefei Slurm templates are named `cluster/hf_prl2023_hfo2_*.slurm`.

No VASP binary, POTCAR, author-owned archive, or active run directory is
distributed here.  The current workflow uses 20 total images, fixed cached
endpoints, isolated interior workers, exact POSCAR cell serialization, and a
native VASP lattice-consistency probe before electronic calculations.  A
completed scheduler job is not accepted as a reproduced barrier until force,
stress, phase identity, energy normalization, and path mechanism are audited.

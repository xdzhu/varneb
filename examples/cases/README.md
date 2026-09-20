# Compact case index

The reviewed material fixtures remain in their historical, provenance-rich
directories.  This index is the stable entry point for new users:

| Case | Backend(s) | Purpose |
| --- | --- | --- |
| BaTiO3 | VASP, QE candidate | simple cross-backend static/preflight comparison |
| GaN B3/B1 | VASP | literature path and finite-pressure barrier |
| CdSe RS/WZ | VASP | G-SSNEB mapping and small-barrier audit |
| HfO2 T/PO/M | ABACUS, VASP | variable-cell production and acceleration benchmark |
| Optional backend smoke (H2) | QE, LAMMPS, CP2K, ABINIT | calculator energy/force/stress contract only |

Open the corresponding directory under `examples/` and the linked plan under
`docs/` before submitting a job.  Do not create a new material case merely to
test a backend; use the compact BaTiO3 or toy case unless a literature
comparison requires otherwise.

# BaTiO3 VASP static-image input template

This template is for a fixed-geometry BaTiO3 T-to-C VARNEB image.  It uses an
explicit `ENCUT = 600 eV`, `4×4×4` Gamma-centred mesh, `EDIFF = 1E-8`, and the
static VASP contract (`IBRION=-1`, `NSW=0`, `ISIF=2`, `ISYM=0`).  `600 eV` is a
declared starting point for the reviewed `Ba_sv/Ti_sv/O` PBE PAW set whose
largest `ENMAX` is 400 eV; it is not a conversion of the ABACUS cutoff.

`POTCAR` is intentionally absent.  Use
`scripts/setup_batio3_vasp_static_case.py` to copy a user-licensed Ba/Ti/O PBE
PAW file into a new work directory and to reject endpoints that do not exactly
match the accepted ABACUS reference.

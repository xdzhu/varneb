# Original GaN/VASP production endpoints at 45.7 GPa

These are the *original* image-0 B4 and image-28 B1 static inputs and raw
outputs paired with the 29-image VASP chain used for the manuscript. They
come from the completed hf migration/static job 27719609, immediately before
the completed production VCNEB job 27719610. Their source directories were
`.../gan_b4_b1_vasp_pbe_paw_qian/hf_static_initial/00` and
`.../gan_b4_b1_vasp_pbe_paw_qian/hf_static_final/28`. The archived OUTCAR
SHA-256 values are
`3c6795acf28dad6d00dd33a6f99cd755bf0b0cf8c9c4490bcc1fb3de3005dc87`
and `5b69743b754b170807ada26d9ed147de3deb83fc6bc7b8749434b2f7cf98a671`.

The two statics have the same `INCAR` and `KPOINTS` bytes: VASP/PBE,
`ENCUT=600 eV`, `EDIFF=1e-7`, `ISYM=-1`, `SYMPREC=1e-4`, `PREC=Accurate`,
and a Gamma-centered 8x8x6 mesh. Their remote POTCAR SHA-256 values both
equal the production image-15 POTCAR hash
`f94781ce6cf9b9454c093353320c5e80a2a0aceea6fb8353b33b01ac37b95168`.
POTCAR is not redistributed. The static `INCAR` has `IBRION=-1` and `NSW=0`;
45.7 GPa is the *analysis* pressure for `H=E+PV` and the raw stress residual,
not a `PSTRESS` setting in these static calls.

Reproduce the hash-bound audit from repository root:

```bash
python -m scripts.audit_gan_45p7_vasp_endpoints \
  --output /path/to/fresh/endpoint_audit.json
```

The committed result is
`benchmarks/numerical_integrity/gan_45p7_vasp_endpoint_static_audit_20260930.json`.
It independently parses the two completed/SCF-converged OUTCARs and checks
their exact input POSCAR geometries, energies, forces, and stresses against
the archived evaluated chain. Both endpoint atomic-force maxima are below
0.02 eV/A. The **production** B4/B1 raw-stress residuals against 45.7 GPa
are **2.91154/0.23863 kbar**: B4 misses the diagnostic 2-kbar gate; B1
passes. The separate signed basin-return structures have different residuals
(1.743/2.528 kbar) and must not be substituted for these original endpoints.
No endpoint was reoptimized, no production parameter was changed, and this
audit alone cannot predict a changed VCNEB barrier.

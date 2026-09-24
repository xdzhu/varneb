# VARNEB 0.0.2

This release packages the current calculator-independent VCNEB core and
documents the completed five-backend GaN B4→B1 validation at 45.7 GPa.
The CP2K path converged at step 52 with 29 total images and a maximum
generalized force of 0.09409 eV/Å; its 0.29276 eV/GaN barrier and source
data are included in the CPC manuscript package. The manuscript now contains
the converged cross-backend figure, BTO mode analysis, and explicit evidence
and limitation statements.

The case library gains a compact CP2K GaN fixture and a production caution:
multiple independent CP2K MPI jobs must have verified, disjoint CPU affinity.
The Slurm template defaults to one CP2K worker and rejects unverified
multi-worker launches. This does not change the physics or silently alter an
existing run; users must recheck endpoint and static-force/stress contracts
on their own installation.

The release workflow publishes only on an explicit version tag or manual
dispatch. Ordinary source pushes and GitHub Release creation do not publish.

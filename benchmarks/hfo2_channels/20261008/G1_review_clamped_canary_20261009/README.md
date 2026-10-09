# E045: G1 candidate review and one clamped PO+ canary

This finite experiment advances the registered G1-to-G2 gate; it does not
extend the scientific channel/material/boundary matrix. Four terminal
candidate edges share the same ordered PO+ initial well and original six
physical-input hashes. The review replays all forty E/F/stress observations,
endpoint physical screens and the three declared symmetry tolerances,
mode/strain/reconstruction diagnostics, forward-minus-reverse energy
identities, thirteen existing peak statics and fourteen-interior Berry audits.
No extra DFT is required for that review. HF additionally rechecks the original
forty raw logs, all six physical bytes and parsed E/F/stress.

The review is **candidate connectivity and evidence readiness**, not a
continuous activation-barrier error bound, certified TS, complete stability
test, distinct-MEP/winding classification or exhaustive escape network.
The lower Pbcn snapshot remains a coverage/stability caveat. T's existing
two-amplitude Gamma calculation has no negative optical modes in its sampled
atomic subspace; it does not certify all-q or joint atomic/cell stability.
These limitations remain research questions for the registered G2/G3 tests,
not reasons to invent an additional endpoint or Hessian matrix.

## Preregistered canary

Only `clamped_endpoint_seeds/strain_0000/PO_plus/endpoint_seed.json` is selected.
It is a fresh, unrelaxed geometry: no free-cell E/F/stress cache is reused.
The proper cyclic new-x/y/z=old-z/x/y rotation and T-referenced substrate
were registered before this experiment. The first two cell rows are fixed;
the third row's length and two tilts and all atomic positions are released.

- BFGS, exactly at most four steps; at most five fresh SCFs.
- Atomic gate 0.03 eV/A, open-traction-norm gate 2 kbar, maxstep 0.02 A.
- Fixed PBE/100 Ry/full10auDZP, original INPUT/KPT/Hf/O pseudopotential and
  orbital bytes. The electronic settings are never changed to reach a phase.
- One HF/hfacnormal01 allocation, 32 genuine MPI ranks, one thread per rank,
  direct `mpirun`; 30-minute resource cap, no other batch submitted automatically.
- A step-limit canary with complete raw, byte and boundary audits is a healthy
  bounded test, **not** an endpoint-convergence result. Unexpected exceptions
  or incomplete summaries still fail. Clamped reaction stress is recorded but
  excluded from the open-traction convergence criterion.

The generic relaxation entry's nonconverged exit1 is retained. The dedicated
canary Slurm wrapper accepts it only after independent replay of the complete
step-limit summary and all actual SCFs. It does not hide an electronic failure.
No phase-restoring symmetry or mode constraints, CI, restart, +0.5% holdout,
new peak statics or new scientific channel follow automatically.

## Offline commands

Use fresh output paths. HF raw review requires the recorded source directories;
local review uses frozen numeric observations and exported raw sampling/Berry data.

```sh
python -m scripts.audit_hfo2_G1_gate --output /new/G1_review.json
python -m scripts.audit_hfo2_G1_gate --raw-on-HF --output /new/G1_review_HF.json
python -m pytest -q tests/test_hfo2_G1_gate.py tests/test_hfo2_clamped_canary.py tests/test_hfo2_clamped_endpoints.py tests/test_epitaxial_boundary.py tests/test_hfo2_fixed_input_factory.py
```

The synthetic refusal fixtures test the protocol, not HfO2 physical results.
Submission/completion receipts are added only after observed execution. All
production and analysis snapshots remain separate and immutable. A canary
pass must be reviewed before a bounded endpoint continuation or the existing
G2 matrix; the reserved +0.5% path labels remain unseen.

## Actual G1 review and submission

`G1_review_local.json` is the completed local numeric replay. Fresh
`G1_review_HF.json` additionally verifies forty original terminal raw logs and
all six physical-file bytes; it passes the **candidate** gate with zero DFT.
`canary_preflight.json` rechecks ten registered seeds and the actual physical
source without creating/evaluating calculators. The data are not predictions.

The clean production snapshot `fd738b1734` passed 1316 tests with two skips
before submission; focused tests passed 63. Runtime archive/source and actual
receipt hashes are recorded in `validation_delivery.json` and `submission.json`.
Job **28446324** was submitted at 21:54:01 CST and observed RUNNING on node6
with 32 allocated CPUs at 21:54:39. This observation is not completion evidence.
No other G2/G3/holdout calculation is automatically launched.

`G1_review_archive.json` is the fresh clean-archive numeric replay; it equals
the HF report except its explicitly narrower raw-check scope. The older
`G1_review_local.json` is preserved. Its Gamma audit digest binds the working
LF file, whereas the Git archive retains that historical JSON's CRLF bytes.
The exact pair differs in 2255 carriage returns and no numeric/text content;
both hashes and the failed initial strict comparison are recorded in the
receipt. No generic hash exception or physical-input mutation is introduced.

## Actual bounded canary result

28446324 completed0:0 at22:07:34CST after812allocation seconds. Five genuine
32MPI SCFs cost756.647228 measured subprocess seconds/6.725753coreh; the
allocation cost is7.217778coreh. All six physical bytes and complete E/F/stress
were audited on HF. The final atomic max is0.037333755eV/A, open traction
6.948494kbar, so the correct physical status is **step_limit, not converged**.
All five geometries retain Pca2_1 at .001/.01/.05A/1degree tolerances and the
fixed plane drift is at floating-point round-off. This is healthy boundary/
transport evidence, not Hessian/phase stability, electronic polarization or
G2 channel evidence.

`completed_HF/` contains exactly49 observable files, no pseudo/orbital/charge
redistribution. Original HF audit/summary/raw bytes are retained. Stock ASE's
missing ABACUS-I/O registration exposed a local replay limitation; a separate
post-DFT analysis change adds a strictly case-writer-specific Direct/Hf4O8
reader. It reproduces the actual data without changing/repeating DFT or
overwriting the original runtime. The first failed local replay is recorded,
not counted as a passing check. The original global ignore exception omitted
the endpoint subdirectory; only this exact case's observable paths are now
allowed, not arbitrary calculation outputs.

Next registered finite experiment is [one geometry continuation](../clamped_PO_continuation_E046_20261009/README.md):
max20new steps/SCFs, exact last-clamped-SCF reuse, fresh BFGS Hessian and
independent raw/phase review before any wider G2 work.

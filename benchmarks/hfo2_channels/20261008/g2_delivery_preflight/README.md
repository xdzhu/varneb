# G2 delivery preflight — E036, 2026-10-09

This is an **inert deployment check**, not a clamped HfO2 calculation, phase
certificate, scientific prediction or permission to skip G1. Ten previously
registered starters remain unrelaxed; no endpoint energy/force/stress is reused
from free-cell data. The +0.5% holdout is still absent.

## Frozen delivery

The source was archived from commit `1e9259cddcb540b5c13accb01bc3fb19ec080528`,
not the dirty checkout. The 22,725,666-byte `source.tar.gz` has SHA256
`14420324da90df38d9d5bf65321e2ff117abfeb57bc0dc29d2a2dbdc23437dac`.
It was extracted locally outside Git and into the new HF namespace
`/public/home/iai806/abacus/agent-runs/20261009-varneb-G2-delivery-preflight-1159/source`.
No running source or existing case was overwritten. The standalone
[validator](check_preflight.py) was transferred separately with SHA256
`c61b94316b0173104b2d7c196227a061501da32ae4a87f129f1a70fd0637099a`;
it is not silently claimed to be in that older source archive.

The archive's runtime/test directories (`vcneb`, `scripts`, `examples`,
`tests`, `cluster`) are unchanged from the E032 commit with the recorded
1,162-pass/2-skip full clean regression. That older full suite is not rerun or
claimed as a new result. The present archive independently passes **48**
existing endpoint/epitaxial/input tests in **1.64 s** (JUnit 1.642 s), with no
failures, errors or skips. These use analytic/synthetic/EMT fixtures, not DFT.

## Actual HF checks

The [captured HF report](preflight_HF.json) and
[delivery receipt](validation_delivery.json) preserve the source and input
hashes. Two read-only executions use the existing ICU Python/ASE 3.23.1b1 /
NumPy 1.26.4 environment; no software is installed. They verify:

- All ten seed and manifest checksums, ordered Hf4O8, valid periodic geometry,
  common substrate at each condition, tilt-released third vector and P=0.
- The +1% substrate is the declared scaled zero-strain substrate. Maximum
  plane drift is about `5e-17 A`, a floating-point comparison, not a DFT
  precision or strain-error estimate.
- All six original INPUT/KPT/UPF/full 10-au DZP hashes. Calculator construction
  leaves results empty and never creates the reserved endpoint directory.
- The existing Slurm template passes `bash -n`; it is **not executed**.

The helper never calls energy/force/stress evaluation or `sbatch`. Small
read-only validation does not launch ABACUS on a login node. A future real
evaluation still requires a reviewed HF/hfacnormal01 allocation with 32 MPI,
one thread per rank, unchanged PBE/100 Ry/full 10-au DZP and original KPT.

## Reproduce without starting DFT

Archive/extract the named commit in a fresh namespace, then use the separately
versioned helper. Its physical-source parameter file points to the existing
audited HF inputs; those inputs are not redistributed in this case.

```bash
env PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  /public/home/iai806/.conda/envs/icu/bin/python check_preflight.py \
  --source /path/to/fresh/extracted/1e9259c/source \
  --probe-root /path/to/a/nonexistent/reserved_endpoint
```

Local focused regression in the extracted archive:

```bash
python -m pytest -q tests/test_hfo2_clamped_endpoints.py \
  tests/test_epitaxial_boundary.py tests/test_hfo2_fixed_input_factory.py \
  tests/test_relax_ase_endpoint.py
```

The initial local launcher used an as-yet nonexistent working directory and
was rejected **before process creation**. Creating/extracting from the repository
first, then invoking pytest in the clean source, fixes that delivery ordering.
It was not a DFT, numerical convergence or source-code failure.

## Next authorized step, still gated

After every G1 candidate has a valid terminal path and its common-initial,
input/SCF, endpoint/variant and sampling audit, review a four-step BFGS canary
for the existing zero-strain PO+ starter. Its template is already in this
archive. No canary, G2 matrix, +0.5% holdout, CI or new channel is submitted by
this milestone. A phase that disappears after real relaxation must be reported,
not restored with a mode/symmetry constraint. This check reduces deployment
risk; it supplies no JCTC material or prediction verdict.

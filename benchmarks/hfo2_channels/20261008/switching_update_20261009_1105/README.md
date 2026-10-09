# E034: a lower central snapshot, not a certified intermediate

2026-10-09. Two fully evaluated observations of the existing G1 switching
allocations are frozen here: preserving step59 and reversing step22. No DFT,
restart, cancellation, endpoint optimization or input change was performed.
The original production archive remains byte-identical. Both jobs are live.

The common ordered PO+ energy is -9783.249675811956 eV/Hf4O8 cell, four formula
units, P=0/E=0. ABACUS/PBE/100Ry/full10-auDZP, original mesh/pseudos/orbitals,
ordinary 0.10eV/A, no CI and the production periodic lift are unchanged.
Every exported image matched a completed genuine DSIZE=32 raw SCF with finite
energy/forces/stress and the exact physical-input and ordered-geometry hashes.
The numeric-only replay reproduced the production force and energy log.

| Common-PO view | Frozen step | Ordinary residual (eV/A) | Sampled maximum (meV/f.u.) | Ordinary pass |
|---|---:|---:|---:|---|
| PO to T, reused reverse view | 6 | 0.059881616 | 115.210164 | yes |
| PO to M, reused terminal observation | 39 | 0.097904490 | 71.582207 | yes |
| T-pattern-preserving flip | 59 | 0.156442033 | 35.045730 | no |
| T-pattern-reversing flip | 22 | 0.310501161 | 394.632486 | no |

These are not final competing barriers, a complete G1 pass, stationary saddle
certificates or matched-boundary predictions. Nineteen earlier T--PO/M image
records are reused; the18 newly frozen records include four cached endpoints.
The37 network records are neither37 new SCFs nor independent replicates.

## What changes the next interpretation

The preserving chain retains two sampled maxima at images3/5, approximately
35.045730/35.045728meV/f.u. Its central image4 is now **-1.605722meV/f.u.**
relative to the identical PO+ endpoint, compared with +17.440944 at step48.
It is Pbcn at each registered symmetry tolerance (0.001/0.01/0.05A, angle1deg).
This is an actual lower-energy snapshot, not a stable nonpolar intermediate:
central atomic/cell residuals are0.111271/0.156442eV/A and transverse stability
is unmeasured. Its raw stress_xx is-0.005887307eV/A3. Near-zero tangential
force due to the symmetric chain does not remove these perpendicular forces.
The old central peak must not become the preselected TS or initial state of
an extra path. Final basin/branch and sampling audits remain required.

Its deformation relative to the registered original T chart has F_xx=1.110022
and Green strain_xx=0.116075. The already registered clamped plane contains
that original x direction (new x/y/z = old z/x/y). A finite, prospective
interpretation to check in the existing G2 matrix is therefore whether the
same substrate suppresses/replaces this expanded central branch, rather than
only shifting a single saddle. This observation alone neither proves that
expansion causes the basin nor establishes an independent forecast. It does
not add a Pbcn endpoint, another path, strain or Hessian to the finite budget.

The reversing chain is still atomic-dominated at image2. Its sampled central
maximum is Pbca at0.001/0.01A but Pa-3 at0.05A. Keep the complete sensitivity
grid: neither selecting the desired tolerance nor a vanishing T-triplet
projection identifies a phase or a local soft mode. The central residual is
only0.050829eV/A, yet the whole chain is0.310501 and no joint negative-direction
count exists. It is not a full-variable first-order TS certificate.

## Reproduction and validation

```sh
python -m scripts.analyze_hfo2_chain_observations \
  --root benchmarks/hfo2_channels/20261008/switching_update_20261009_1105 \
  --variants benchmarks/hfo2_channels/20261008/reference_variants \
  --gamma benchmarks/hfo2_channels/20261008/gamma_analysis/T_d0.01.npz \
  --output /new/path/analysis.json
python -m scripts.analyze_hfo2_network_progress \
  --specification benchmarks/hfo2_channels/20261008/switching_update_20261009_1105/network_specification.json \
  --output /new/path/network.json
python -m pytest -q tests/test_hfo2_chain_observation.py tests/test_hfo2_observation_analysis.py tests/test_hfo2_network_update.py tests/test_hfo2_network_progress.py
```

No implementation or tests changed. The clean E032 source passes the35
existing focused regression tests (58.74s) and reproduces both new local
reports byte-for-byte. HF/local two-frame replay matches1298 floating fields
within2.842171e-14 and249 nonfloating fields exactly; four-channel replay
matches697 floating and882 nonfloating physical/source fields with the same
maximum difference. The original reports retain three package-version
differences and38 structure-warning fields; these are not physical agreement.
No remote pytest or new full-suite run is claimed. Raw bytes are protected
against Git newline conversion;27numeric evidence files match their raw Git
blobs and an extracted staged archive byte-for-byte. Later delivery metadata
does not change those files. See [the receipt](validation_delivery.json).

At11:17:53CST both same job handles remain RUNNING on hfacnormal01: preserving
28319570/node11 has a later log step60/0.151932; reversing28319571/node26 has
step23/0.300858. Those rounded live logs do not replace the frozen observations.
Seven recent step intervals are roughly11--15min. Continuing recent decay
would put preserving near the ordinary threshold around13--15h and reversing
around15--18h CST today; this is a conditional trend estimate, not a convergence
promise. The preserving segment has at most20 steps left after60 and must be
audited if its80-step cap is reached without convergence. Do not stop it on a
single rebound or restart a still-live job. G2/holdout remain unsubmitted.

## Bounded literature follow-up

The [Behara publisher entry](https://journals.aps.org/prmaterials/abstract/10.1103/PhysRevMaterials.6.054403)
explicitly requires a subscription for supplemental material. Its linked SI
endpoint fails and the author lab entry provides no SI link in the retrieved
page. No SI/full-paper completion claim, credential, access bypass or author
contact is made. Nature's two Lee article/PDF URLs also failed direct access;
their DOI/version relationship remains unresolved. These limitations do not
block the live calculations and are not evidence of novelty.

# HfO2/BaTiO3 validation error budget and resource efficiency

This is a reporting-level budget, not a statistical uncertainty estimate. It
separates numerical path convergence, electronic-structure sensitivity, and
endpoint residuals so that a quoted barrier is not given more precision than
the evidence supports.

## HfO2 T-to-PO

| contribution | evidence | scale |
| --- | --- | ---: |
| ordinary image-count spread | 7 vs 9 total images, both audited | 0.0005417 eV (0.35%) |
| cell interpolation control | 7-image linear vs log-strain, same endpoints/settings | 0.0029263 eV (1.87%) |
| CI refinement shift | ordinary 7-image to CI 7-image | 0.0275787 eV |
| endpoint residual | T/PO independent BFGS endpoints | max force 0.000499 eV/A; max stress 0.07076 kbar |
| geometry validity | independent path audit | min distance 2.026--2.036 A; no audit issues |
| five-image control | 5 total images, extended rebound window | not converged; excluded from barrier comparison |

The ordinary barrier is therefore reported as 0.1562--0.1568 eV at the
tested image resolution; the CI value 0.1291722 eV is a separate saddle
refinement, not an uncertainty bar around the ordinary value.

The completed production runs used 32 MPI ranks per interior image and fixed
endpoint caching:

| run | Slurm job | total/interior images | workers | wall time | result |
| --- | ---: | ---: | ---: | ---: | --- |
| ordinary 7 | 27678218 | 7 / 5 | 5 | 01:08:38 | 0.1567509 eV |
| ordinary 9 | 27678407 | 9 / 7 | 7 | 03:47:07 | 0.1562092 eV |
| CI 7 | 27678507 | 7 / 5 | 5 | 03:12:03 | 0.1291722 eV |

These wall times are allocation-level elapsed times, not a claim of perfect
parallel efficiency. Worker manifests are the authoritative record of image
ownership and successful evaluations.

## BaTiO3 T-to-C

The ordinary forward 5/7/9-image paths are all barrierless and give the same
reaction enthalpy, +0.0871289 eV. A tighter 6x6x6, SCF 1e-9 endpoint/path
check shifts the absolute reaction energy by 0.01476 eV but preserves the
barrierless topology. The CI gate is therefore withheld.

## Reporting rule

Do not combine the CI shift, image-count spread, and electronic-precision
shift into one formal error bar: they probe different calculations and
different physical/path choices. Report them separately and retain the raw
summaries, audits, and worker manifests.

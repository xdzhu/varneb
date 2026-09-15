# HfO2 T-to-PO VCNEB Closure

This note closes the current HfO2 T-to-PO production VCNEB pass at 100 Ry with
full 10 au DZP Hf/O orbitals.  The endpoints are fixed 12-atom Hf4O8
conventional-compatible cells and are evaluated once/cached during each VCNEB
run; only interior images are dispatched to ABACUS image workers.

## Accepted production evidence

| Branch | Job | Images | Resource model | Final force (eV/A) | Barrier (eV) | Audit |
|---|---:|---:|---|---:|---:|---|
| Ordinary log-strain | 27678218 | 7 total / 5 interior | 5 x 32 MPI workers, 2 nodes | 0.0455155873 | 0.1567509 | ok |
| Ordinary log-strain | 27678407 | 9 total / 7 interior | 7 x 32 MPI workers, 2 nodes | 0.0497190884 | 0.1562092 | ok |
| Ordinary linear-cell control | 27678924 | 7 total / 5 interior | 5 x 32 MPI workers, 2 nodes | 0.0495086804 | 0.1596772 | ok |
| CI refinement | 27678507 | 7 total / 5 interior | 5 x 32 MPI workers, 2 nodes | 0.0295475450 | 0.1291722 | ok |

The ordinary 7/9-image barrier spread is 0.0005417 eV, so the ordinary path is
image-count stable at the present resolution.  The CI refinement is the
accepted climbing-image result for the same endpoint convention and production
calculator settings.

## Rejected diagnostic branch

The endpoint-displacement mode-guided branch is retained only as an
unconverged mechanism diagnostic:

| Job | Resume source | Final force (eV/A) | Provisional barrier (eV) | Audit |
|---:|---:|---:|---:|---|
| 27679830 | initial mode-guided path | 0.1545267023 | 0.0627620 | geometry valid, force failed |
| 27683214 | step-43 continuation | 0.0882049964 | 0.0043842 | geometry valid, force failed |
| 27687189 | continuation of 27683214 | 0.0866792767 | 0.0038557 | geometry valid, force failed |

For job 27687189, the independent audit found no geometry problem
(`minimum_path_distance_A = 2.0349939195`, `maximum_deformation = 0.0492007537`,
`maximum_stress_kbar = 1.7036551219`) and failed only because the generalized
force remained above the 0.05 eV/A target.  Its final 25 FIRE steps were
effectively plateaued near 0.087--0.089 eV/A, so the small apparent barrier must
not be promoted to a physical HfO2 barrier.

## Conclusion

For current HfO2 T-to-PO reporting, use:

- ordinary VCNEB barrier: 0.1562--0.1568 eV from the converged 7/9-image
  log-strain pair;
- CI-refined barrier: 0.1291722 eV from job 27678507;
- reaction enthalpy: -0.3252848 eV for the shared endpoint convention.

Do not cite the mode-guided 0.0039 eV branch as a converged mechanism.  Further
work on true material mode constraints should start from a physically motivated
mode basis or release-and-refine protocol, not by extending the plateaued
endpoint-displacement diagnostic indefinitely.

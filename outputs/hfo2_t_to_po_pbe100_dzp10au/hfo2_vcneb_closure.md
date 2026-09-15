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

## Loose-converged diagnostic branch

The endpoint-displacement mode-guided branch does not satisfy the stricter
0.05 eV/A target used for the main ordinary VCNEB comparison, but its final
continuation does satisfy a conventional loose NEB threshold of 0.10 eV/A:

| Job | Resume source | Final force (eV/A) | Barrier (eV) | Audit |
|---:|---:|---:|---:|---|
| 27679830 | initial mode-guided path | 0.1545267023 | 0.0627620 | geometry valid, force failed |
| 27683214 | step-43 continuation | 0.0882049964 | 0.0043842 | geometry valid, below 0.10 eV/A |
| 27687189 | continuation of 27683214 | 0.0866792767 | 0.0038557 | geometry valid, below 0.10 eV/A |

For job 27687189, the independent audit found no geometry problem
(`minimum_path_distance_A = 2.0349939195`, `maximum_deformation = 0.0492007537`,
`maximum_stress_kbar = 1.7036551219`).  Its final 25 FIRE steps were effectively
plateaued near 0.087--0.089 eV/A, which is above the stricter 0.05 eV/A target
but below a loose 0.10 eV/A NEB criterion.  The resulting 0.0038557 eV barrier
may be cited only as a loose-threshold endpoint-displacement mode-guided result.

## Real-material strict-subspace to full-space release diagnostic

Job 27693064 first found a 7-total-image strict-subspace path; its small force
is a projected-force quantity and is not comparable to an unconstrained force.
Job 27693085 then resumed that chain with `constraint_mode=none`, fixed/cached
endpoints, five 32-MPI interior workers and all atomic-plus-cell degrees of
freedom released.  The 300-step full-space continuation is geometrically valid
(`minimum_path_distance_A = 2.0222188115`, `maximum_deformation = 0.0491735878`,
`maximum_stress_kbar = 4.3742647607`).

The full-space generalized force reached its global minimum of `0.097528 eV/A`
at step 275.  Following the rebound policy, the run continued for 25 more
steps, ending at `0.108376 eV/A`; it therefore fails the strict `0.05 eV/A`
target.  Under a conventional ordinary-NEB `0.10 eV/A` criterion, the archived
step-275 chain is a loose acceptance snapshot with a `0.1380114 eV` barrier
(`34.50 meV/f.u.`) and highest image 2.  It is not mixed into the strict
ordinary/CI table.  The full trajectory, log, selected snapshot and audit are
stored in `vcneb_n7_mode_subspace_release_job27693085/`.

## Conclusion

For current HfO2 T-to-PO reporting, use:

- ordinary VCNEB barrier: 0.1562--0.1568 eV from the stricter 7/9-image
  log-strain pair;
- CI-refined barrier: 0.1291722 eV from job 27678507;
- endpoint-displacement mode-guided diagnostic barrier: 0.0038557 eV from job
  27687189, only under a loose 0.10 eV/A NEB force criterion;
- strict-subspace-to-release loose snapshot: 0.1380114 eV from step 275 of job
  27693085, only under a loose 0.10 eV/A ordinary-NEB convention;
- reaction enthalpy: -0.3252848 eV for the shared endpoint convention.

Do not mix the 0.0039 eV branch into the stricter 0.05 eV/A ordinary-image
convergence table.  Further work on true material mode constraints should start
from a physically motivated mode basis or release-and-refine protocol, not by
extending the plateaued endpoint-displacement diagnostic indefinitely.

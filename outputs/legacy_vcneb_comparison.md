# Legacy VCNEB implementation comparison

The historical implementation was inspected under
/home/zhuxd/abacus/8.dielec/vcneb during the 2026-06-18 validation run. The
directory is no longer present on the current hf filesystem, so this record
uses the contemporaneous iteration log as provenance and does not invent a
quantitative rerun.

## Recovered comparison

| topic | legacy implementation (recorded) | current VARNEB |
| --- | --- | --- |
| Extended path variables | Fractional coordinates plus cell degrees of freedom | Fractional coordinates plus deformation-gradient cell variables, with explicit cell scale |
| Endpoint mapping | Species-wise assignment and MIC mapping were available | Audited species-preserving mapping, translation alignment, distance/deformation/fold checks |
| Optimizer | ASE-compatible path updates | ASE FIRE/BFGS/LBFGS through one public runner; ordinary path precedes CI |
| Endpoint treatment | Not documented as a separately cached immutable evaluation | Endpoints excluded from optimizer DOFs and evaluated/cached once; manifests prove interior-only dispatch |
| Stress/cell force | Useful design direction was retained | Row-vector UnitCellFilter convention with finite-difference tests, pressure tests, and non-orthogonal-cell coverage |
| Failure/recovery | Historical run notes exposed periodic-MIC path discontinuity and spurious barrier risk | Calculator preflight, exact-state cache, atomic snapshots, sibling preservation after batch failure, and failure reports |
| VASP layout | Per-image VASP directories were used | Static per-image VASP calculators with isolated directories and explicit IBRION=-1, NSW=0, ISIF=2, ISYM=0 |

## Failure-mode evidence

The historical log records a mic=True endpoint reset that reintroduced a
periodic discontinuity and produced a spurious large HfO2 barrier. The current
implementation keeps the shortest fractional branch through the endpoint and
audits the resulting geometry before calculators are launched. This is a
failure-mode comparison, not a claim that the two implementations produce
identical production barriers.

## Boundary

No old production output or executable is currently available for a
like-for-like numerical rerun. The project therefore keeps this comparison
qualitative and leaves the quantitative legacy-result benchmark open rather
than marking it complete.

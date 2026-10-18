# One finite strain transfer of the lower preserving candidate

Frozen before any new target-interior SCF, 2026-10-19.

## Question and scope

Does the transverse-accommodated, connected zero-strain E273 candidate
remain accessible at +0.5% and +1%, and does its sampled competition with
the existing T escape change? This is a mechanism/initialization control,
not a new independent prediction or a search over strain, sign or material.
The original labels and forecast scores remain unchanged.

Use exactly the final nine-image E273 step9 band, original job28949435,
ordinary residual0.0940120033eV/A. Transfer this one selected source to the
two already studied boundaries. Reuse the E268 +0.5% and E145 +1% endpoint
geometries and exact native caches. No endpoint relaxation or recalculation.
Use the existing E268 geometric recipe: change the two substrate rows,
then apply linearly distributed endpoint fractional/open-row corrections.
Preserve ordered atomic identity and the source's continuous periodic lift;
reject ambiguous registration instead of remapping atoms or fitting a gauge.

## Fixed execution

- ABACUS3.10.0LTS/PBE, Hf4O8,100Ry, full10auDZP, the same six fixed assets,
  Gamma2x2x2, P=0/E=0; no electronic-input changes.
- Two independent ordinary bands, nine total/seven free interior images.
  Local F/A/J constraints are not imposed on these free interiors.
- FIRE,k=.2,maxstep=.02,maximum-cell-step=.02,minimum-distance1.6,
  maximum-deformation.25,allow-tilt, original fixed metric5.12968067458423A.
- Ordinary NEB residual target0.10eV/A, no CI. Each chain: at most60updates,
  427fresh interior calls (seven initial plus seven per update),18hours.
- Each allocation: hfacnormal01,one node,one manager task32CPU/task,
  images serial, every new image true mpirun32,one thread/rank.
  At most two study allocations; other projects are untouched.
- New immutable directory/source. Zero-SCF native endpoint/input preflight;
  submit held, check exact identity/resources, then release. No duplicate
  submission, implicit retry, runtime overwrite or automatic budget extension.

## Prespecified interpretation and stop

Report both full sampled profiles, residuals, geometry, transverse structural
proxy, open row and all call costs. Compare each strain against its original
preserving candidate and T escape at that same strain/initial well. Keep the
original two +0.5% initialization controls visible, not just the lowest band.
If the transverse feature vanishes, report relaxation back toward another
candidate rather than insisting it defines a new branch. A lower ordinary
sampled peak is not a certified TS or global MEP. Close differences remain
unresolved without an error bound; residual0.10 is not an energy uncertainty.

No new strain, opposite-sign run, local release, Hessian, fitted model or
forecast follows automatically. A cap/SCF failure supplies a finite outcome,
not permission to modify physics or resubmit. These two chains answer the
remaining across-strain connection question; then write the supported physical
conclusion and explicit limits. E274 remains a completed four-static protocol.

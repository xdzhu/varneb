# Local strain work and endpoint-only response controls

`vcneb.strain_work.configuration_work` contracts any calculator's ASE energy
gradient convention with an explicit configuration direction. It does not
run DFT, relax structures, select a mode, certify a branch or freeze forecasts.

For ASE row-cell matrix H, fractional positions s and r=sH, let t be an
explicit control parameter. With tensile-positive stress sigma and fixed
external pressure P, the derivative is

```text
d(E+PV)/dt = V (sigma+P I) : [(dH/dt)^T H^-T]
            - sum_i F_i dot [(ds_i/dt) H].
```

The cell term is evaluated at fixed fractional coordinates. Supplying
`fractional_direction` adds only internal-coordinate motion; do not also
include the affine cell displacement there or it is counted twice. Cells
are Angstrom, cell directions Angstrom per unit t, forces eV/Angstrom, and
stress and P eV/Angstrom^3. Outputs are eV per unit t. There is no factor-of-two shear,
mass weighting or VARNEB cell-scale multiplier in this physical work.
Symmetric full stress tensors are required; convert an ASE Voigt vector
with `ase.stress.voigt_6_to_full_3x3_stress` first.
Use a continuously lifted fractional-coordinate direction with the registered
atom ordering; differences of independently wrapped positions can introduce
spurious jumps. This function does not choose a periodic lift or remap atoms.

For biaxial strain epsilon, fix the open third vector for the local partial
and set dH[0:2]/d epsilon = H[0:2]/(1+epsilon), dH[2]=0. The general matrix
formula also works for oblique or globally rotated cells; simply summing
two chosen Cartesian diagonal stresses does not in general define this
intervention. If the released variables are exactly stationary, this partial
is the envelope derivative of that stationary branch. Finite force/open
traction residuals must instead be reported, not silently set to zero.

```python
import numpy as np
from ase.stress import voigt_6_to_full_3x3_stress
from vcneb.strain_work import configuration_work

dh = np.zeros((3, 3))
dh[:2] = atoms.cell.array[:2] / (1 + epsilon)
work = configuration_work(atoms.cell.array,
    voigt_6_to_full_3x3_stress(atoms.get_stress()), dh)
```

This example assumes an already evaluated atoms object; this function itself
never calls its calculator. At nonzero fixed P use `pressure_eV_A3=P` and
compare enthalpy, not raw energy alone. Varying external P would require an
additional V dP/dt term, which is outside this API.

## Actual HfO2 training analysis, not a holdout prediction

[E059](../benchmarks/hfo2_channels/20261008/endpoint_strain_work_E059_20261010/README.md)
reuses all ten original raw-audited endpoint representations at0/+1%, P=0,
the common T substrate and unchanged ABACUS100Ry/full10auDZP. It computes:

- actual well-energy changes and local imposed-plane partials;
- two-endpoint stress-work trapezoids and their observed defects;
- a configuration-chord diagnostic including released-cell and atomic
  residual work, using the existing ordered endpoint lifts;
- the two already registered B1 well-only *training response increments*.

A trapezoid defect mixes unmeasured curvature, finite-relaxation residuals
and numerical effects. It is not a rigorous uncertainty or barrier error bar.
The chord is not an optimized intermediate phase branch. Neither diagnostic
licenses tighter electronic parameters or new intermediate DFT calls.

B1_fixed gives delta B=-delta E_PO+. B1_follow gives
delta B=delta E_FS-delta E_PO+. Without a matched audited anchor barrier,
these are null response implications, not absolute barriers or prospective
forecasts. B0 and B2--B5 remain unavailable here; the complete prediction
freeze and unopened+0.5% holdout rules in
[v1](HFO2_PREDICTION_PROTOCOL.md)/[v2](HFO2_PREDICTION_PROTOCOL_V2_2026-10-09.md)
are not bypassed. The existing [prediction-control API](PREDICTION_CONTROLS.md)
still performs its stricter full-network/anchor/error gates.

Replay from the repository root into a new output file:

```sh
python -m scripts.analyze_hfo2_endpoint_strain_work --output /new/seen_training_endpoint_work.json
python -m pytest -q tests/test_strain_work.py tests/test_hfo2_endpoint_strain_work.py
```

Finite differences of an independent scalar energy with mixed atomic/cell
directions, pressure, oblique-cell rotations and shear test the work formula.
Actual material fixtures check only the seen endpoint data and the explicit
non-forecast boundary, not a missing TS response or method advantage.

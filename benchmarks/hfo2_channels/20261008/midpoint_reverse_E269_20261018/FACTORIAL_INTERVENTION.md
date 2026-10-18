# Minimal frozen atomic–cell intervention

## The question

The initial midpoint preserving configurations differ in both oxygen shuffle
and open-axis stretch. Their energy contrast therefore does not identify
which internal arrangement, cell response or combination accounts for it.
Use a four-state frozen intervention before attributing a difference to a mode.
This is a standard finite interaction contrast, not a new reduction theorem.

For each channel, A is the completed zero-to-midpoint ordinary peak candidate.
B will normally be the upper-to-midpoint terminal peak candidate. An explicit
`--initial` option can prepare diagnostic initial-snapshot decks for review,
but does not admit their execution or promote them to terminal results.

## Four states, only two new SCFs per channel

| Internal fractional pattern | Cell H_A | Cell H_B |
| --- | --- | --- |
| s_A | AA: existing native energy | AB: one new static SCF |
| s_B | BA: one new static SCF | BB: existing native energy |

Both cells have the exact same prescribed substrate within serialization
tolerance. The third vector, including its allowed tilts, is exchanged.
Ordered unwrapped fractions are the internal coordinates in the common parent
chart. Holding them during a cell exchange means affine deformation; it does
not hold laboratory Cartesian positions fixed. Species, atom order and path
registration remain unchanged. A large fractional registration change is
refused, never repaired by atom matching.

All four state energies use the original fixed ABACUS/PBE input, full10auDZP,
100Ry, original six asset hashes, Gamma2x2x2, Hf4O8 and P=E=0. The same common
initial-well energy cancels from every finite contrast. The new mixed states
are frozen SCFs, not relaxation runs or new NEB paths.

At cell A, the internal-pattern contrast is E_BA−E_AA. At pattern A, the
cell contrast is E_AB−E_AA. The finite interaction is

`E_BB − E_BA − E_AB + E_AA`.

The diagonal contrast equals those three signed contributions. Reversing
the order of the two interventions gives different conditional contributions
when the interaction is nonzero. No unique elastic energy, phonon energy or
mode-by-mode Hamiltonian decomposition is inferred. Cell energy in eV is
multiplied by250 for meV/HfO2, because this cell contains four formula units.

## Admission and limits

- The current two E269 allocations remain untouched; no third job is submitted.
- Only if completed E269 results justify the contrast will terminal-source decks
  be generated in a fresh namespace and undergo native HF preflight.
- Existing AA/BB statics are reused only at identical geometries and input bytes.
- At most four new static evaluations across preserving and T; no grid, model
  fitting, derivative stencil, cutoff change or endpoint recalculation.
- A failed or missing corner is unmeasured, not zero. All four measured energies
  are needed for the interaction; no value is inferred from the two diagonals.
- This frozen contrast cannot replace fixed-retained-coordinate F/A/J release,
  demonstrate a stable conditional surface or certify a saddle. If it provides
  no useful distinction, stop this intervention rather than expanding its matrix.

The preparer writes only four small INPUT/KPT/STRU decks plus a manifest.
Pseudopotential/orbital bytes remain on HF. Local native-case-reader geometry
checks are recorded separately from the still-required actual HF/ASE preflight.
No launcher is generated, and `execution_admitted` remains false.

## Actual preparation and replay

Four initial-snapshot diagnostic decks were generated under ignored `build/`
with zero SCFs/submissions. They are not a selected production experiment.
The source code supports terminal exports by default and rejects accidental
use of initial snapshots. Reproduce the diagnostic preparation with:

```sh
python benchmarks/hfo2_channels/20261008/midpoint_reverse_E269_20261018/prepare_factorial.py --exported benchmarks/hfo2_channels/20261008/midpoint_reverse_E269_20261018/initial_complete_chains --output NEW_REVIEW_DIRECTORY --initial
```

Tests exercise actual source geometries and identical INPUT/KPT hashes,
substrate/fraction preservation, no overwrite, missing-energy refusal and
four-formula-unit signed-contrast algebra. Synthetic arithmetic fixtures
are not material results. No crossover energy has yet been calculated.
The latest scoped run passed44tests, including the earlier native comparison
and figure tests. This does not constitute a DFT validation of the crossover states.

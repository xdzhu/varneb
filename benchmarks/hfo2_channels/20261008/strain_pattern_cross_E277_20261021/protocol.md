# Final two strain–pattern cross controls

Author authorization: 2026-10-21. This closes the current HfO2 study; it is
not a new material, model-selection, strain-grid or Hessian campaign.

Question: does the local response to the imposed plane depend on the
retained shuffle state, independently of changing both together?

Let q0 and q1 be the **four** original geometric observables (Qx,Qz,Qy,
mean O-minus-Hf z displacement in the original T metric) at the registered
preserving peaks at strain 0 and +1%. They are not phonon eigenmode energies.
Reuse the previously qualified **minus-sign J** candidates at (0,q0) and
(+1%,q1); the sign is selected before new energies are observed. Compute
only the missing (0,q1) and (+1%,q0) conditions. Keep all four observables
fixed, and release the other 29 atomic and three permitted cell directions.
Three rigid translations remain gauge-fixed. No new A-stage, sign pair,
strain midpoint, complete path, force-constant matrix or stability claim.

Each cross start is the donor diagonal J, with the target prescribed plane
and donor open cell row; fractional atomic coordinates are transported with
one recorded per-atom integer registration from the native F structures.
No atom rematching or image-wise wrapping is allowed. A rigid translation
is removed as a gauge. Independent original-T projections must reproduce
the donor four-observable vector. The open row is an initial guess, not an
extra fixed physical constraint. The new Hamiltonian is unchanged.

ABACUS/PBE: original six asset hashes; 100 Ry; complete 10-au DZP; Gamma
2x2x2; original SCF 1e-8/nmax150/mixing0.3; zero pressure and field.
Use the immutable E268 source and existing 32-MPI mpirun executable/env.
hf/hfacnormal01 only; at most two one-node/32-CPU allocations, 2 h each.
No 235/cu17, CI, parameter adjustment or running-source replacement.

Numerics: safeguarded BFGS; original target-chart branch/on-site numerical
matrix; initial trust radius 0.05 A; local released-gradient target
0.03 eV/A; 1.0 A coordinate-radius guard, minimum distance 1.6 A and maximum
deformation 0.25. At most **40 fresh evaluations per cross**, including
the initial point and failed attempts: 80 total. Cache only exact coordinate
matches. No automatic budget/radius extension or resubmission. An unmet
target is missing evidence, not a fabricated candidate. Numerical guards
and the smaller historical zero-strain guard are disclosed.

Compare U(e,q)=250*[E_candidate(e,q)-E_IS(e)] in meV/HfO2, using each
strain's existing common initial well (four formula units). Define
I=[U(1,q1)-U(1,q0)]-[U(0,q1)-U(0,q0)]. The initial-well shift cancels
in I. Report the two fixed-q strain responses and two fixed-strain pattern
contrasts, plus transverse displacement and open-row geometry. These are
finite local candidate comparisons, not a unique atomic/elastic energy
partition, global conditional minima, certified TSs or new path barriers.

The experiment ends at the two original jobs' terminal results. Positive,
negative or inconclusive results all close this protocol. Then update the
mechanism argument, main figure, reusable case, manuscript and reading PDF,
run relevant tests and push. Do not replace old independent-prediction labels.

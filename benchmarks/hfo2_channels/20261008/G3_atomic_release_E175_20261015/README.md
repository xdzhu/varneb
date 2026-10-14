# E175: prospective atomic-only treatment ablation at the same five Q points

This is the B3 release-treatment part of the registered B2--B5 comparison,
not a replacement for the complete T/Cmma reference or independent-prediction
matrix. Prior E173/E174data and failures are already seen training evidence.

## Question and physical comparison

Do the observed frozen-to-joint response and its original-quadratic mismatch
arise mainly during atomic release, or during additional open-cell release?
For each original T/preserving center use exactly Q=(0,0),(-.1,0),(0,-.1),
(0,.1),(.1,0)Angstrom. Retain the same first-step two-dimensional coupled
plane, input contract and original center. No new Hessian center, whole band,
material, strain condition or physical+.005geometry is introduced.

At each Q, the frozen seed is x0=AQ. Release only the31-dimensional atomic
space satisfying A.T*delta_x=0; the three extra open-cell displacement rows
are identically zero, not merely rounded below a tolerance. Thus the cell is
fixed **to that Q's frozen seed**, not to one common cell across all Q. This
is not a pure atomic-mode coordinate plane or a new fixed-parent reference.
The full physical gradient is pulled back from the unchanged40D chart.
The native map is anchored at the **literal existing frozen full coordinate**:
X(u)=X_frozen+B*(u-u_frozen). Thus the first callback has exactly the cached
geometry even when cross-machine matrix multiplication rebuilds B*u_frozen
with a different last bit. Q/embedding consistency retains the existing
2e-12arithmetic check; cache matching itself remains exact. New displacements
have strictly zero cell rows. This is not a tolerance-based physical cache hit.

Compare the resulting candidate EA(Q) with existing EF(Q) and EJ(Q), keeping
all signed differences and separate origin offsets. Finite local optimizer
candidates need not realize the globally nested-minimum inequality. No sign
clipping, last-trial promotion, post-hoc domain shrink or unique microscopic
cause is assumed. Raw-positive seed blocks at both original probe steps are
only algebraic seed checks, not stability certification after release.

## Fixed finite reservation and execution

- Five additional training-node treatments per existing face, including Q0;
  T29→34/81, preserving30→35/81; remaining47/46. All8independent points per
  face and the physical+.005holdout remain untouched.
- Each node at most12new attempts including failures,60per face,120total.
  The two original centers and eight completed frozen seed SCFs are cached;
  the first callback must hit that exact seed or refuse before DFT.
- Original safeguarded BFGS, original first-step Hessian restricted to31D,
  .03eV/A Euclidean **atomic released** gradient target, .05A initial trust,
  .5A total internal radius and1.6A minimum separation. This is not a full
  joint-gradient pass or changed ordinary NEB. No eigenvalue clipping or
  pseudoinverse. One rebound is not a stopping condition.
- One manually admitted allocation processes both centers serially: HF,
  hfacnormal01,1node/32genuineMPI,mpirun-np32,OMP/MKL/OpenBLAS1,5h maximum.
  At most two committed study allocations including pending/held. No
  auto-resubmission or background monitor. Prior calls near113s suggest
  at most roughly3.8h of SCFs if all120are used, plus preflight/overhead;
  queue delay is separate. This is a cap estimate, not predicted convergence.
- Original ABACUS/PBE/100Ry/full10auDZP/Gamma2³/P0/E0/six hashes/binary and
  old380file source snapshot. Ordinary NEB remains0.10; no CI or INPUT edit.

`worker.py` binds the exact E173controls and native admission reports, the
prepared pure-atomic basis/projector, points, code, physical hashes and caps
in `config.json`. The five-file external bundle is immutable and the old
E161runtime is reused without overwrite. Native preflight reparses all8
frozen SCFs, cached centers/origins and input/runtime identities. The
separate one-shot driver saves intent before sbatch, records the accepted
handle, validates held resources, then releases that same handle. An
ambiguous intent never triggers another submission.

This treatment does not consume independent labels or resolve missing
post-release stability, final predictive representatives/domains, full
strong-reference comparison or numerical barrier-error evidence. Those
remain requirements of the active JCTC goal. No favorable result is presumed.

The first external r1preflight refused before any native DFT or submission
because it demanded bitwise equality from a regenerated frozen embedding.
That source/config and the failed receipt remain archived externally. The
r2bundle uses the explicit affine anchor above, with a1ULP basis-variation
regression. It is the same single E175reservation, not another five-node
treatment or a changed physical seed/input. No running source is overwritten.

The r2native preflight passes: all8original frozen SCFs reparsed, original380
runtime identities checked, and10inert native input decks retain exact seed
cells, prescribed-control zero,17digit STRU round trips and original6asset
hashes. Local/independent4-module checks each pass73tests,16new, not a full
regression. The independent tested tree is79338ebb10a4fc89c854ef5ca3df714c0da36faa.

Actual one-shot job28867177is submitted held, validated1node/32CPUs/1task/
5h/hfacnormal01, then released. Its first queue snapshot is PENDING, not a
failed calculation or evidence of convergence. Intent/accepted/held/released
receipts are retained. No production result from this new treatment is yet
claimed. The original preserving face continues separately, unchanged.

At06:00:19CST the new handle is confirmed PENDING/AssocGrpCpuLimit, with0
new SCFs/points and no failure output. Do not resubmit or stop other jobs
to bypass that queue condition. The original preserving28862441is RUNNING,
23/25processed with17gradientcandidates. No complete-face claim follows.

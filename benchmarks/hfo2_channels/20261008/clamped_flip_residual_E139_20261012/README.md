# E139: native-cache diagnosis of two live +1% switching chains

These are complete evaluated snapshots, **not terminal or converged labels**.
On 2026-10-12 the accepted preserving child28737465 and reversing child28764704
remain RUNNING on hfacnormal01, one node/32CPU each. Only preserving step37 and
reversing step15 are exported; later production updates remain separate.

The unchanged E103 driver (SHA256
6c64c0155800082e28cbbb280511ac761929bdb6ccc5757fa8d68d156b32eac7)
loads its pinned E093 audit support, verifies the original338-file E055
runtime, and checks the accepted journal, actual native identity and submitted
environment. The E103 wrapper hash remains
ab914119aaa7918bd0cb466ed8dd86daafecb20609a077a5d249d89a3228eed0.
Only an audit-process factory registration is added; production source and
parameters are not edited. This is observation export, not a terminal audit.

The exporter checks all nine same-ordered geometries, exact six-file physical
contract, native DSIZE32, electronic convergence and full E/F/stress. Residual
replay agrees with the production log. The portable bundle retains eighteen
raw SCF records, INPUT/KPT/STRU and source audits, both complete trajectories,
and original observation/preflight records. No pseudopotential or orbital is
redistributed. Local replay checks the retained raw records; it does not
claim to re-audit every historical SCF or the undistributed HF assets/runtime.

## Result and interpretation

| Snapshot | Ordinary fmax (eV/Angstrom) | Limiting image/block | Same-image open-cell max | Spring-only max over chain |
|---|---:|---|---:|---:|
| preserving37 | 0.117604576798 | 4 / atoms | 0.059113965811 | 0.005878231381 |
| reversing15 | 0.143674141082 | 4 / atoms | 0.010533619599 | 0.005345109346 |

At the actual limiting vector, signed spring projections are respectively
2.8811661e-7 and3.6067163e-12eV/Angstrom. Perpendicular projections are
0.117604288681 and0.143674141079. Component maxima cannot be added; the
same-vector signed projections reproduce the ordinary residual. Atom vectors
are generalized-coordinate forces conjugate to the registered reference-cell
coordinates, not untransformed Cartesian forces. Cell norms use the original
5.12968067458423Angstrom scale.

This excludes spring dominance **at these snapshots**. It does not demonstrate
that another optimizer, spring constant, redistribution or mode constraint
would accelerate convergence. Nor does it certify a saddle, its Hessian index,
a channel identity or barrier precision. Both fail the ordinary0.10 gate.
No bottleneck/reference/mode selection, G3 or physical+0.005holdout is admitted.

Preserving steps27--37 decrease0.149167→0.117605; reversing steps5--15 increase
0.129821→0.143674 despite decreasing logged maximum path energy. A scalar
energy decrease does not imply decreasing perpendicular NEB residual. These
logs alone do not establish the cause of the force rise or optimizer failure.
They support continuing the finite segments, not stopping on a rebound.

The two native observation hashes are respectively
a4593fcf2c0f1613cced5e884afcf58518b1285c38ffd3aa9336f043e5e4ce31 and
c26a3c8f783c8c758fdf5716bb3aea6dff2885f5f19c5fcde3629ca43f55282f.
No DFT is initiated by this analysis; authorized production continues normally.
No scheduler write, endpoint recomputation, cutoff/SCF/optimizer change or
running-source replacement. No paper numeric table or figure is updated with
these unconverged labels. G2 remains6/8; the broad research goal remains active.

## Replay

```sh
python -m pytest -q tests/test_hfo2_E139_live_residual.py
```

Actual local and independent archive regressions each pass123tests, zero
failures/errors/skips. The verified archive tree is
4fa93092c4ad0c5458104efaa712075eca82d462; all6013tracked files are present
before independent tests start. All115new raw-evidence/JUnit files are byte
equal; the four source/README/attribute files are newline-normalized equal.
The first comparison stops on metadata CRLF/LF before any independent test;
raw evidence is never normalized or changed. The same completely extracted
archive passes the corrected strict raw-byte/source-text comparison.

Archive SHA256:f9170d7ac918c62cd22b12a2be889d3fa05554dd0dc23f2dabe1dc1b040f3ee7.
Local JUnit:eb9b5589b28415733fd9f2f75750be6b8e53e9ca95a7639b7cea18ac74321715,
7.990s. Independent JUnit:15fc307c288802ffa16435b1aa66100ff103e1cb7b4138abac84bc9e5502bf98,
7.635s. This is not a new full-suite or HFpytest run. Only this final validation
metadata, the copied independent JUnit and the session ledger are added after
the tested snapshot; scientific data, analyzer and test source are unchanged.
The114native portable files total1,423,411bytes. Original five dirty tracked
files and old untracked user artifacts remain untouched and unstaged.

Ordinary terminal and all-call auditing must still precede the next material
gate; finite40-update termination is not convergence.

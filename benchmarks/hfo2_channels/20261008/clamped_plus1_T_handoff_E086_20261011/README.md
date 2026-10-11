# E086: prepare the registered +1% T-to-M rolling handoff

This is an operational preparation, not a material convergence result. It
does not submit a new job, change dependencies or physical settings, claim
that job28709788 has converged, or select G3 modes/bottlenecks. The target
is the already registered, user-held +1% PO-to-M starter28709789.

`release_existing_M_once.py` separates two actions:

- `--preflight-only NEW_REPORT` checks the original registered M seed,
  endpoint caches, physical inputs,338 runtime Python files, ABACUS binary,
  pilot and three protected held jobs. No terminal receipt is required or
  claimed; no DFT or scheduler write is made.
- `--reviewed-audit-sha256 HASH` additionally requires an independently
  exported/audited ordinary-converged **actual +1% T terminal** receipt,
  native Slurm COMPLETED0:0 and an actual spare study slot. It exclusively
  creates an intent journal before releasing only28709789 once. Any
  existing accepted or ambiguous journal requires reconciliation, not retry.

The flipped starters28709790/28709791 remain held and keep their original
M dependency. Slurm removing a fulfilled T dependency is permitted; removal
of the still-pending M dependency is refused. No dependency is rewritten.
The companion28723655 may still be running or be normally completed, but
at most two study allocations are allowed. Unrelated jobs are not modified.

The original six physical bytes,100Ry/full10auDZP/PBE/P0/Gamma2x2x2,
32MPI/thread limits, ordinary0.10eV/Angstrom and no-CI policy are unchanged.
The+.005 prospective holdout is neither generated nor read. This helper
cannot authorize further chains or parameter changes.

## Review and test scope

The local related suite includes mock native dispatch, an ambiguous accepted
reply, resource/dependency/capacity races, changed original source/seed,
incorrect boundary, nonfinite or above-threshold force, and an already
existing production directory. These are mocked transaction tests, not new
material calculations or native failure-injection experiments. Native
read-only preflight and independently archived source execution are recorded
separately in `validation.json`; no fresh whole-project regression is claimed.

Deploy this private helper only in its new analysis namespace. It imports
the hash-pinned E083 read-only preflight functions; it must not execute the
old E083 admission routine or overwrite the running E055 runtime. Follow-up
terminal receipts and admission must be reviewed against live state first.

The actual native read-only preflight passed, direct process exit0. It checked
9 geometries, two native cached endpoint SCFs with original six inputs,
DSIZE32/full EFS,338 runtime Python files, binary/pilot and all three held
jobs. `pending-preflight.json` byte-matches the HF receipt. Both local and
actual independent Git-archive related suites passed172tests, including46
new mock handoff tests. Local/deployed versus archive Python bytes differ
only in CRLF/LF serialization; normalized sources agree.

At16:31:39, +1% T remained RUNNING, step4/.111567eV/Angstrom: this is **not**
an ordinary convergence pass. Preserving remained RUNNING, step28/.170654.
Both actual Slurm stderr files were empty. G2 remains3/8; this preparation
made zero releases, submissions, dependency edits or new DFT calls. No
admission journal existed. The material terminal gate still has to pass.

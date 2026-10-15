# E190: real native records expose a provenance-path replay bug

2026-10-15. Previous E189 was the tested/pushed manuscript milestone.
The existing jobs now run on HF: atomic-only28867177 started09:42:10 on
node6, released-curvature28868827 started09:44:39 on node9, each32CPUs.
Native logs confirm genuine DSIZE32 and one thread per MPI rank. No new
allocation, input change, resubmission, cancellation or compute elsewhere.

## Root cause, not an electronic-structure failure

The original E183 prefix reader rejected the first real derivative record.
Tracing that exact function showed identical coordinates, energy, gradient,
stress, volume, minimum distance and maximum force, but a different validation
dictionary. The stored raw path contained `bundle/../released_curvature`,
whereas the reader replayed its canonical spelling. Strict remote resolution
proved both paths identify the same directory. The mismatch was provenance
spelling, not a failed SCF, bad physics or a numerical-tolerance problem.

The older smoke read had constructed its result path from the stored raw
path, so both spellings matched. Scanning new records from the canonical
producer namespace exposed the missing case. The added regression exercises
that exact caller/provenance combination, not only the old smoke route.

E183 and both producers remain unchanged. E190 is a separate external reader
revision: after the SAME strict resolved raw-directory equality and namespace
checks, replay the validator using the stored `raw.parent`. This preserves
the original validation dictionary. Input/binary hashes, exact coordinates,
reservations, energies, stress checks and all2e-12 tolerances are untouched.
No cache/record rewriting, symmetry repair, cutoff increase or rerun.
The new regression test verifies this is the only executable-line change.

Use the corrected deployed reader for later prefix/complete audits:
`/public/home/iai806/abacus/agent-runs/20261015-varneb-G3-released-audit-E190-r1/audit.py`,
SHA51d0452968394ff03d26efface7705ee31c7679215cb2204f6c69e13a5391d4f.
Keep original E161 PYTHONPATH and numerical threads1. Reports require fresh
paths outside the producer namespace. Complete mode still requires all272
native probes and normal one-node32CPU Slurm completion. Old E183 is kept
as provenance, not the recommended live derivative reader.

## What actually passed

The revised HF prefix reader independently reparsed13 real derivative
records, retaining the pending next probe and unstarted preserving block.
Its 10:11 snapshot is not the latest running count or a curvature matrix;
no positive/negative index or stable conditional branch was inferred.
Native prefix SHAb80d7571848f151f326f2c19576b8f9d2446800fcea776808bcc06aa97255c52.
The original B3 reader separately passed the first TQ0 candidate with
one new SCF and exact cache-only optimizer replay. Atomic released norm
0.02147688266eV/A is below .03, but joint norm0.11042982399 is not: the
atomic-only arm is not claimed to be joint-stationary. This is one of ten
points, not a completed batch, stable minimum or predictive success.

`first_probe_native.zip` preserves7 raw result/input/STRU/log/audit files
plus their source manifest, no pseudo/orbital/binary. SHA
d7bd98e90ca97969dd208729739703eba2a00ceed4dc7bb5f902e155e0640e32.
The archive's real E/F/stress/gradient can be rechecked locally; the strict
in-place provenance/namespace audit was performed on HF, not simulated by
relocating the archived record. The raw first-record SHAc35035e2... is unchanged.

71 local/71 independent related checks across4modules,13new, zero failures/
errors/skips; not fullsuite. Synthetic controls reproduce the old false
rejection, accept only the identical resolved directory, and still reject
changed E/gradient/stress/geometry/MPI/coordinates and another real directory.
The real captured SCF is also parsed, rather than relying only on mocks.
Tested tree7435743446e8419b6a2905100311aa0312399fb8/archive
10808dee384327c8f5e065670fb8a8460c8e249de8286060b1918f9659ff8821.
All scratch/checkouts are on C:; no private author bank access.

## Running progress and conditional timing

At10:16:44CST both jobs remain RUNNING, no SCF failure records. The atomic
producer reports6/10 gradient candidates and15 completed SCFs. Its47 remaining
attempt ceiling excludes unused budgets from already-finished nodes: they
cannot roll over. At the observed131s/SCF mean, using that entire remaining
ceiling would finish near12:00; this is a cost scenario, not a convergence
or time guarantee. Its walltime limit is14:42:10.

Curvature has17/272 completed probes, about108s/SCF mean. Uniform remaining
cost gives about17:54; different later probes can change this estimate.
Its walltime limit is21:44:39. Later audits must take a new snapshot.

Original100Ry/full10auDZP/physical boundary/.10NEB/.03conditional target,
272probe budget, trainingT34/P35 commitments, eight independent reserved
points/face computed0 and physical+.005 seal remain unchanged. Manuscript
unchanged. Full native batches, stable released blocks and independent
material forecasts still remain; complete goal active/incomplete.
The ABACUS workflow skill guided raw MPI/EFS/input checks; its generic235
examples did not override the authorized HF-only scope.

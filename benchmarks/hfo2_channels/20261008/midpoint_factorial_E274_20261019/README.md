# E274: four frozen midpoint atomic/cell crossovers

Prepared from the **completed** E268/E269 sampled peaks, not the earlier
unconverged snapshots. Exactly four new static SCFs are planned; none has
yet been run or submitted. See the [finite protocol](protocol.md).

The preserving starting-band contrast is10.820844meV/HfO2; T's is
-0.131935. The open-cell row and internal pattern both change. Exchange
them in the same ordered fractional chart, keeping the two existing AA/BB
energies. AB/BA are frozen geometries, not relaxed candidates. This
distinguishes conditional internal and affine-cell effects without a
unique phonon/elastic energy partition or a new independent forecast.

All four native INPUT/KPT bytes remain identical to the original contract.
No pseudo/orbital files are distributed. The geometry preparer checks the
same substrate, ordered fractional registration, minimum separation and
native STRU roundtrip. Actual HF/ASE preflight, held resource validation
and accepted-handle release remain separate required gates.
The native writer's periodic wrapping is explicitly distinguished from
the continuous analysis lift. Formal calls must match the STRU bytes
registered by the native preflight, not an earlier serialization's hash.
This was caught before any submission or DFT call, without changing physics.

One32CPU hfacnormal01 allocation, one hour, serial mpirun-np32 and one
thread/rank; no CI, new NEB, relaxations, retry, model fit or extended grid.
Run code is frozen in a new namespace; prior runtimes are read-only.
Original labels, unsuccessful releases and preparation costs remain.

`control.py` implements native preflight, duplicate-resistant held dispatch,
the four-call compute-only entry and C0-only portable collection. The
existing package helper copies only this control, run script, protocol and
prepared decks into a byte-pinned bundle. No running source is replaced.

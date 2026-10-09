# Terminal ordinary-force passes, not full G1/TS certification

Checked on HF2026-10-09 16:59:39 CST. Both actual allocations COMPLETED0:0.
Exported complete terminal snapshots are matched to every original ordered
periodic raw SCF, original six physical-file hashes and freshly parsed full
energy/forces/stress. The optimizer log residual and highest energy replay
within the declared formatting tolerance. Original sources remain unchanged.
Local numeric/hash/geometry replay independently passes. Zero new export SCFs.

| Candidate | Job / final step | Total / moving | Ordinary fmax eV/A | Sampled maximum relative PO meV/HfO2 |
|---|---|---|---:|---:|
| T-pattern-reversing flip | 28319571 /45 | 9 /7 | 0.094022845666 | 392.822905197 |
| PO→M sampling refinement | 28392675 /1 | 12 /10 | 0.096964012700 | 82.679957306 |

Common PO energy−9783.249675811956 eV/Hf4O8,4f.u.,P0/E0, ordinary0.10,
noCI. M refinement changes the static reconstructed peak82.747071 to the
ordinary-band sampled82.679957; it does not erase the original nine-image
71.582207 observation. Its single step is not an acceleration benchmark.

R allocation06:29:53--16:32:13,36140s at32CPUs,321.244444coreh. M allocation
16:19:58--16:38:30,1112s at32CPUs,9.884444coreh, below20step/8h cap.
These allocation costs include setup/SCF/optimizer overhead, not only SCF.
The21frozen image records reuse endpoint/seed calculations and must not be
called21new or independent SCFs. HF source/raw audits and actual byte hashes
are recorded in the two observations. Large original orbital/UPF data stay
onHF. Original electronic settings, pressure, image lift and endpoint order
are not changed to match these outcomes.

All four original candidate bands now pass their ordinary-force threshold;
P preserving and T/PO also have earlier static sampling checks. This is **not
full G1**: reversing peak sampling, path/variant/polarization and mechanism
audits remain required. An ordinary residual does not establish stationarity,
transverse stability, a sampling error bound, global channel coverage or an
independent strain prediction. No G2/canary/holdout is submitted here.

In particular, the preserving band's lower-energy Pbcn midpoint is not yet
a certified stable intermediate. If it is a stable nonpolar escape basin,
a comparison limited to M/T does not establish global escape resistance.
The existing branch/stability gates must resolve this limitation; no extra
Pbcn endpoint, Hessian or new path is silently added to the finite matrix.

ObservationSHA256:

- R:`23c1f7e3799ab0f1b388bd2791418c400f26eccf5679a3353c05dff326ef3f08`.
- M:`d2ad5af18b18820cc92a3176d727feda13881d80dd0c91d682ee10805566766f`.

Replay without HF/DFT:

```sh
python -m pytest -q tests/test_hfo2_terminal_G1_update.py
```

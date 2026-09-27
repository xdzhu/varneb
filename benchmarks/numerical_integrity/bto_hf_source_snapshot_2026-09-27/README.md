# Exact hf BTO source snapshot for the September 2026 conditional probes

The hf work root
`/public/home/iai806/abacus/agent-runs/20260927-varneb-bto-soft-gridpilot-qx030/code`
is a frozen **non-Git** source snapshot. It was used for the BTO conditional
canaries/relaxations and the `27787309` Hessian probes. Later local development
changed several files, so the current package revision must not be presented
as the byte-identical production executable. This directory preserves the
five differing Python files at their original relative paths, without
overwriting the current implementation:

| Snapshot file | Source SHA-256 |
| --- | --- |
| `examples/run_bto_q1q2_conditional_abacus.py` | `de0c1c3aeac376b9e4f4ea1f9ff742e85e70f9ba08823824c46c358d547a17bf` |
| `scripts/audit_bto_q1q2_conditional_pilot.py` | `afe51053a1c23f119b256246cbe21e97529c6da3c9e037dec45895b360708f29` |
| `vcneb/mode_surface.py` | `91b70d48a488b5a26f1e563bbce1ba0e75581a01d9145158d57108b8dac4bbfc` |
| `vcneb/__init__.py` | `25576062e4c39ea92e552fa78c9e7dd554ce00d8d7667d8a7404eb58245c1868` |
| `vcneb/core.py` | `c69e4066e314fa5f81b30c8eb3f66e3dfe327d9d62abdd718369a7932f920169` |

The relevant reference, frozen-grid runner, evaluator, ABACUS adapter,
phonon, reference-cell, and modes modules matched the local files at the
time of the targeted SHA-256 audit. The static pseudopotential/orbital
assets and ABACUS binary remain identified by hashes in the preflights and
raw result manifests. Source files and remote JSON artifacts were copied
without edits and checked against the source hashes. Later read-only Hessian
and soft-direction audits used separately hashed scripts under the hf
`audit-code-20260927/` directory, not this frozen production tree.

This snapshot is evidence, not a replacement import package and not a claim
that all later `vcneb` features existed when these DFT points were run.

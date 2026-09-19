# GaN VASP Bravais regression structures

Exact POSCAR geometry from the two cu17 failures (VASP 6.3.2). The original
files contained Cartesian coordinates at 16-digit precision. Only the comment
and whitespace have been normalized; no metric, positions, shear or ordering
have been changed. No licensed POTCAR or VASP source is distributed.

Remote root: `/home/zhuxd/abacus/agent-runs/20260916-varneb-v/gan_b4_b1_vasp_pbe_paw_qian/`.

- image 15: `recovery_isym_minus1/vcneb_b4_to_b1_tetragonal_n29_cu17_serial/15/POSCAR`.
- image 24: `recovery_isym_minus1_symprec1e8/vcneb_b4_to_b1_tetragonal_n29_cu17_serial/24/POSCAR`.

Both must pass a non-mutating geometry gate; neither is an invalid cell.
Calculator-free checks cannot certify VASP's Bravais classifier. See
`docs/vasp_input_contract.md` for the actual initialization probes and limits.

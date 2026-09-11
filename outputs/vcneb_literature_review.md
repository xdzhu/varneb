# VC-NEB Literature and Software Survey

Date: 2026-09-12

## Executive conclusion

Variable-cell NEB has already been implemented in several forms. The most direct
VC-NEB implementation is Qian et al. (2013), integrated into USPEX. Closely
related implementations are the generalized solid-state NEB (G-SSNEB) in TSASE,
FLOS/Siesta, QuantumATK, and ABINIT, plus the finite-deformation FD-NEB code.

The open problem for this project is therefore not to claim the first
variable-cell NEB algorithm. A publishable contribution can instead be a modern,
calculator-agnostic, reproducible toolkit that:

- exposes the cell/atomic coordinate convention explicitly;
- accepts energy, forces, and stress from several external engines;
- provides robust cell-orientation and atom-mapping preparation;
- distinguishes hydrostatic enthalpy VC-NEB from general-stress G-SSNEB and
  finite-deformation FD-NEB;
- provides restartable cluster workflows, provenance, validation fixtures, and
  cross-engine benchmark data.

This positioning is technically honest and has a clearer software contribution
than re-implementing the published equations without a reproducibility layer.

## Terminology and method families

### Qian VC-NEB

Qian et al. define an expanded configuration space containing nine strain/cell
degrees of freedom and 3N fractional atomic coordinates. Under hydrostatic
pressure, the objective is the enthalpy surface H = E + P Omega. Atomic forces
and the stress-derived cell force are transformed into a common metric space,
then the usual NEB projection and spring force are applied. The paper reports
variable spring constants, improved tangents, climbing image, cell rotation
avoidance, and atom-sequence matching.

The implementation was added to USPEX in 2013. The paper says the code was
available on request; current USPEX documentation still exposes VCNEB and lists
VASP, GULP, and Quantum ESPRESSO as enabled engines. The current USPEX download
terms restrict redistribution and reverse engineering, so its implementation is
not a suitable code base for a freely redistributable toolkit.

### G-SSNEB / SSNEB

Sheppard et al. combine atomic and cell variables in a generalized solid-state
NEB coordinate system with a stress/Jacobian scaling. This is closely related in
purpose to VC-NEB but is not mathematically identical to Qian's enthalpy-gradient
formulation. The distinction matters when comparing paths, barriers, metrics, and
large-strain behavior.

TSASE exposes `ssneb`, `fire_ssneb`, and `qm_ssneb`; it supports external stress,
cell degrees of freedom, CI-NEB, and image-level parallelism. Its code is an old
ASE extension distributed through an SVN tree and is GPL-oriented. It is valuable
as a reference implementation and benchmark, but not an ideal modern foundation.

### FD-NEB

Ghasemi, Xiao, and Gao extend G-SSNEB to finite deformation. Their central point
is that external work must use a stress measure that is work-conjugate to the
chosen deformation variable. Their code supports Cauchy, PK1, PK2, and
hydrostatic stress modes, and provides VASP and LAMMPS examples through ASE and
TSASE.

FD-NEB is particularly important for our toolkit because phase transitions can
have large cell changes. A production package should expose the stress measure
and deformation convention rather than silently treating every stress tensor as
the same quantity.

## Existing implementations

| Implementation | Algorithm family | Engines / interface | Code and practical status | Relevance |
|---|---|---|---|---|
| USPEX VCNEB | Qian VC-NEB | VASP, GULP, QE according to current manual | Mature integrated workflow; registration and redistribution restrictions | Closest scientific reference and strongest baseline |
| TSASE `ssneb` | G-SSNEB | ASE calculators; stress supplied by calculator | Public, old SVN-based library; includes FIRE/QM and parallel modes | Reference for cell-aware ASE design and parallel ownership |
| VTST `LNEBCELL` | VASP SS-NEB / G-SSNEB style | Patched VASP | Requires modifying and recompiling VASP; not stock VASP | Practical VASP baseline, but poor multi-engine portability |
| FLOS `VCNEB` | Qian-style VC-NEB | SIESTA through FLOOK/Lua hook | Public MIT library; requires a SIESTA build with FLOOK support | Demonstrates a non-invasive engine hook, but is SIESTA-specific |
| FD-NEB | Finite-deformation G-SSNEB | VASP, LAMMPS examples | Public GitHub repository; examples are limited and based on old ASE/TSASE patterns | Stress-measure and large-strain benchmark reference |
| ABINIT `neb_cell_algo=1` | GSS-NEB | Native ABINIT | Implemented but rare in tests; variable-cell NEB is constrained to steepest descent | Useful independent implementation for cross-checks |
| ABINIT `neb_cell_algo=2` | Qian VC-NEB | Native ABINIT | Documentation explicitly marks it “Not yet usable” | Evidence that the Qian formulation is subtle in production code |
| QuantumATK | Native GSS-NEB / VC-NEB workflow | QuantumATK calculators | Commercial, current documentation supports `optimize_cell`, target stress, restart | Mature commercial competitor, not an external-engine toolkit |
| ASE `ase.mep.NEB` | Fixed-cell NEB | Any ASE calculator | Current ASE can interpolate cells as a helper, but says cell interpolation is not implemented for NEB calculations | ASE should remain our data model/optimizer base, not the VC-NEB engine |
| ABACUS + ATST-Tools | Fixed-cell NEB / AutoNEB | ABACUS and DeePMD-kit | Current public workflow toolkit; no VCNEB feature found in its documented status | Reuse its adapter ideas where possible; add a distinct cell-aware layer |

## Key primary sources

1. Qian et al., “Variable cell nudged elastic band method for studying
   solid-solid structural phase transitions,” Computer Physics Communications
   184, 2111-2118 (2013), DOI: 10.1016/j.cpc.2013.04.004.
2. Sheppard et al., “A generalized solid-state nudged elastic band method,”
   Journal of Chemical Physics 136, 074103 (2012), DOI: 10.1063/1.3684549.
3. Ghasemi, Xiao, and Gao, “Nudged elastic band method for solid-solid
   transition under finite deformation,” Journal of Chemical Physics 151,
   054110 (2019), DOI: 10.1063/1.5113716.
4. Henkelman group, TSASE solid-state NEB documentation and source tree.
5. Henkelman group, VTST `LNEBCELL` documentation.
6. SIESTA project, FLOS documentation and repository.

## Evidence from application papers

VC-NEB is used in real materials studies, not only in method papers. Examples
include:

- HfO2 thin-film polymorph transitions, including T -> polar orthorhombic
  Pca2_1, using 40-image VC-NEB paths with Quantum ESPRESSO forces/stresses and
  USPEX orchestration. The work explicitly reports strong dependence on
  mechanical boundary conditions.
- Charged oxygen-vacancy-assisted HfO2 polymorphism, where first-principles
  VC-NEB paths are used to analyze kinetic promotion of the polar phase.
- P2-type layered sodium-manganese oxides, where USPEX VC-NEB is used to study
  shear/gliding and tetrahedral transition mechanisms.
- LiV3O8, Fe2O3, carbonates, WTe2, superhydrides, and other high-pressure or
  reconstructive transitions, generally using USPEX VC-NEB or TSASE/VTST
  solid-state NEB variants.

These papers establish that the HfO2 T -> Pca2_1 fixture in this repository is a
scientifically meaningful validation target. They do not, by themselves, give a
portable open implementation or a reproducible multi-engine workflow.

## Software and licensing implications

The USPEX website currently says that recent versions are obtained by
registration, that users may not distribute the code, and that reverse
engineering is prohibited. We should therefore derive the implementation from
the published equations and our own tests rather than copy USPEX internals.

TSASE and VTST are useful references, but their code and license obligations
must be checked before any code is ported. FD-NEB's public repository is useful
for comparison, but its repository should not be treated as reusable under a
permissive license unless an explicit license is present. ASE and ABACUS have
their own open-source licenses and can be used through clean interfaces.

## Position of the current repository

The current prototype already occupies a useful gap:

- `vcneb/core.py` implements a calculator-agnostic optimizer target with
  fractional coordinates plus a deformation gradient.
- `vcneb/vasp.py` and `vcneb/abacus.py` separate engine adapters from path
  mechanics.
- The force/stress transform, cell mask, trajectory resume, snapshot append,
  image ownership, and optimizer API have regression coverage.
- The HfO2 fixture has 12 atoms, mapped tetragonal T and polar orthorhombic PO
  endpoints, and a variable-cell geometry path.

The current synthetic tests are not a DFT validation of the HfO2 barrier. A
real production barrier still needs to be run on the cluster with converged
ABACUS or VASP settings and independently compared with published USPEX results.

## Recommended publication-grade roadmap

### Phase 1: stabilize the scientific contract

1. Define row/column conventions for cells, fractional coordinates, forces, and
   stress in one formal document.
2. Add explicit modes: Qian-style hydrostatic VC-NEB, G-SSNEB/general stress,
   and FD-NEB work-conjugate stress handling.
3. Make rotation gauge, atom mapping, fractional MIC, cell interpolation, cell
   mask, pressure sign, and stress sign explicit in the input and output.
4. Keep FIRE as the robust default and retain ASE BFGS/LBFGS as optional
   optimizers, with warnings for methods that are poorly conditioned for CI-NEB.

### Phase 2: backend and HPC reproducibility

1. Define one backend contract: `energy`, `forces`, `stress`, units, tensor
   ordering, and a structured failure result.
2. Finish VASP and ABACUS subprocess adapters with per-image directories,
   restart-safe output, SCF/convergence provenance, and scheduler templates.
3. Add QE and LAMMPS/ML-potential adapters after the two primary engines are
   stable.
4. Support serial, image-parallel, and scheduler-array execution without
   allowing two images to share a mutable working directory.

### Phase 3: validation ladder

1. Analytic toy surface: exact barrier and finite-difference force checks.
2. Ar fcc -> hcp with GULP or an equivalent smooth potential, matching the
   USPEX example.
3. Si diamond -> beta-tin under compression, comparing VC/GSS/FD stress modes
   and the published FD-NEB trend.
4. HfO2 T -> PO with 12 atoms for the current unit-cell fixture, then a larger
   supercell to test nucleation and the cell-size dependence of the barrier.
5. A perovskite or layered oxide transition with a constrained epitaxial cell to
   demonstrate cell masks and non-hydrostatic boundary conditions.

Every benchmark should publish endpoints, mapping, image count, cell metric,
stress convention, calculator input, convergence history, final path, barrier,
and raw per-image energy/force/stress records.

### Phase 4: publication claim

A defensible paper title would emphasize an open, reproducible, multi-calculator
framework for variable-cell minimum-energy paths, not a first-principles claim
that variable-cell NEB itself is new. The main claims should be:

- a tested implementation of clearly separated VC-NEB/G-SSNEB/FD-NEB modes;
- interoperable VASP and ABACUS execution through a common contract;
- robust structure preparation and restartable HPC orchestration;
- benchmark agreement with established implementations where available;
- improved reproducibility and diagnostics for crystal phase-transition barriers.

## URLs checked

- ASE NEB: https://ase-lib.org/ase/neb.html
- ASE variable-cell discussion: https://matsci.org/t/variable-cell-for-neb/51605
- Qian VC-NEB paper: https://doi.org/10.1016/j.cpc.2013.04.004
- Qian paper PDF: https://uspex-team.org/static/file/Qian-vcNEB-2013.pdf
- USPEX VCNEB manual: https://uspex-team.org/online_utilities/uspex_manual_release/EnglishVersion/uspex_manual_english/vcneb.html
- USPEX downloads and terms: https://uspex-team.org/en/uspex/downloads
- Sheppard G-SSNEB paper: https://doi.org/10.1063/1.3684549
- TSASE SSNEB: https://www.henkelmanlab.org/tsase/ssneb.html
- TSASE source tree: https://theory.cm.utexas.edu/svn/tsase/tsase/neb/
- VTST solid-state NEB: https://vtstools.readthedocs.io/en/latest/neb.html
- FLOS repository/docs: https://github.com/siesta-project/flos and https://flos.readthedocs.io/
- FD-NEB repository: https://github.com/Gao-Group/FD-NEB
- FD-NEB paper: https://doi.org/10.1063/1.5113716
- ABINIT cell-aware NEB: https://docs.abinit.org/variables/rlx/
- QuantumATK variable-cell NEB: https://docs.quantumatk.com/manual/Types/OptimizeNudgedElasticBand/OptimizeNudgedElasticBand.html
- VASP native NEB limitation: https://vasp.at/wiki/Nudged_elastic_bands
- ABACUS documentation: https://abacus.deepmodeling.com/
- ATST-Tools: https://github.com/QuantumMisaka/ATST-Tools
- HfO2 thin-film VC-NEB study: https://arxiv.org/abs/1812.09180
- HfO2 vacancy VC-NEB study: https://arxiv.org/abs/2204.09374
- HfO2 polarization-switching VC-NEB study: https://arxiv.org/abs/2301.06248
- P2 layered oxide VC-NEB study: https://pubs.acs.org/doi/abs/10.1021/acsenergylett.4c03335

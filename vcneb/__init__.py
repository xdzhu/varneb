"""Variable-cell nudged elastic band helpers built on ASE."""

from .version import __version__
from .core import (
    VCNEB,
    VCNEBState,
    apply_chain_state,
    infer_atom_mapping,
    interpolate_vcneb,
    path_geometry_diagnostics,
    read_chain_trajectory,
    run_vcneb,
    validate_atom_mapping,
    validate_path_geometry,
)
from .calculator import (
    CalculatorCapabilities,
    CalculatorCapabilityError,
    classify_calculator_failure,
    inspect_calculator,
    validate_calculator,
    validate_image_calculators,
)
from .modes import (
    Mode,
    build_direction_basis,
    build_mode_basis,
    direction_basis_conflicts,
    mode_guided_path,
    project_path_onto_modes,
)
from .executor import ImageEvaluation, ThreadedCalculatorExecutor
from .qe import attach_qe_calculators, make_ase_espresso_factory, static_qe_input_data
from .phonons import (
    GammaModes,
    diagonalize_gamma_modes,
    load_gamma_force_constants,
    project_displacements_onto_gamma_modes,
    tangent_mode_overlaps,
)

__all__ = [
    "VCNEB",
    "__version__",
    "VCNEBState",
    "apply_chain_state",
    "infer_atom_mapping",
    "interpolate_vcneb",
    "path_geometry_diagnostics",
    "read_chain_trajectory",
    "run_vcneb",
    "validate_atom_mapping",
    "validate_path_geometry",
    "CalculatorCapabilities",
    "CalculatorCapabilityError",
    "classify_calculator_failure",
    "inspect_calculator",
    "validate_calculator",
    "validate_image_calculators",
    "Mode",
    "build_mode_basis",
    "build_direction_basis",
    "direction_basis_conflicts",
    "mode_guided_path",
    "project_path_onto_modes",
    "ImageEvaluation",
    "ThreadedCalculatorExecutor",
    "attach_qe_calculators",
    "make_ase_espresso_factory",
    "static_qe_input_data",
    "GammaModes",
    "diagonalize_gamma_modes",
    "load_gamma_force_constants",
    "project_displacements_onto_gamma_modes",
    "tangent_mode_overlaps",
]

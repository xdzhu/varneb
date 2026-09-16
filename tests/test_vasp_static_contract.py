"""Regression tests for VASP's static VCNEB image contract."""

from __future__ import annotations

import pytest

from vcneb.vasp import REQUIRED_VCNEB_STATIC_PARAMETERS, validate_vasp_static_parameters


def test_vasp_static_contract_accepts_required_parameters_case_insensitively() -> None:
    parameters = {key.upper(): value for key, value in REQUIRED_VCNEB_STATIC_PARAMETERS.items()}
    validated = validate_vasp_static_parameters(parameters)
    assert validated == REQUIRED_VCNEB_STATIC_PARAMETERS


@pytest.mark.parametrize("key, value", [("ibrion", 2), ("nsw", 5), ("isif", 3), ("isym", 2)])
def test_vasp_static_contract_rejects_internal_relaxation(key, value) -> None:
    parameters = dict(REQUIRED_VCNEB_STATIC_PARAMETERS)
    parameters[key] = value
    with pytest.raises(ValueError, match=key.upper()):
        validate_vasp_static_parameters(parameters)

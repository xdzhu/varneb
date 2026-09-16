"""Calculator-independent endpoint-identity regression tests."""

from ase import Atoms
import numpy as np

from vcneb.provenance import compare_endpoint_records, endpoint_structure_record


def test_endpoint_structure_record_is_translation_wrapped_but_order_sensitive() -> None:
    reference = Atoms("HHe", scaled_positions=[[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]], cell=[4, 5, 6], pbc=True)
    wrapped = reference.copy()
    wrapped.set_scaled_positions(wrapped.get_scaled_positions(wrap=False) + [1.0, -2.0, 3.0])
    reordered = reference[[1, 0]]
    record = endpoint_structure_record(reference)
    assert record["sha256"] == endpoint_structure_record(wrapped)["sha256"]
    assert record["sha256"] != endpoint_structure_record(reordered)["sha256"]
    assert record["composition"] == {"H": 1, "He": 1}
    assert np.isclose(record["volume_A3"], 120.0)


def test_endpoint_record_comparison_requires_each_endpoint_to_match() -> None:
    reference = {"initial": {"sha256": "a"}, "final": {"sha256": "b"}}
    assert compare_endpoint_records(reference, {"initial": {"sha256": "a"}, "final": {"sha256": "b"}})["matches"]
    mismatch = compare_endpoint_records(reference, {"initial": {"sha256": "a"}, "final": {"sha256": "c"}})
    assert not mismatch["matches"]
    assert not mismatch["endpoints"]["final"]["matches"]

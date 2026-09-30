"""The four manuscript landscapes must regenerate from archived audit data."""

from scripts.rebuild_cpc_landscapes import rebuild


def test_four_landscapes_rebuild_with_identical_source_tables(tmp_path):
    results = rebuild(tmp_path / "figures")
    assert {item["figure"]: item["DFT_nodes"] for item in results} == {
        "bto_restricted27": 27,
        "bto_frozen289": 289,
        "gan_joint289": 289,
        "gan_transverse162": 162,
    }
    assert all(item["source_csv_byte_identical"] for item in results)

import csv

import pytest

from scripts.plot_gan_gamma_mode_path import validate_source_rows


def _rows():
    return [
        {"reference_phase": phase, "image_index": image, "value": image / 28}
        for phase in ("B4", "B1") for image in range(29)
    ]


def _write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def test_gan_mode_figure_rejects_changed_source_data(tmp_path):
    rows = _rows()
    source = tmp_path / "data.csv"
    _write_csv(source, rows)
    validate_source_rows(source, rows)
    changed = _rows()
    changed[16]["value"] += 0.01
    with pytest.raises(ValueError, match="changed"):
        validate_source_rows(source, changed)


def test_gan_mode_figure_rejects_incomplete_endpoint_basis(tmp_path):
    rows = _rows()
    source = tmp_path / "data.csv"
    _write_csv(source, rows[:-1])
    with pytest.raises(ValueError, match="29 images"):
        validate_source_rows(source, rows)

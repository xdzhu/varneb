import importlib
import pytest


flow = importlib.import_module("benchmarks.hfo2_channels.20261008.switching_path_polarization_20261009.check_export")


def test_native_runtime_does_not_read_only_leading_hour_or_invent_precision():
    assert flow.native_elapsed_seconds(" Total  Time  : 0 h 1 mins 49 secs \n") == 109
    assert flow.native_elapsed_seconds(" Total  Time  : 1 h 2 mins 3 secs \n") == 3723
    for bad in ("Total Time : 0", "Total Time : 0 h 70 mins 1 secs",
                "Total Time : 0 h 1 mins 99 secs", "Total Time : 0 h 1 mins 1 secs\n"*2):
        with pytest.raises(ValueError):
            flow.native_elapsed_seconds(bad)

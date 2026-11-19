import numpy as np

from isoquant_next import __version__
from isoquant_next.core._cy_smoke import fast_interval_sum


def test_package_import_and_cython_smoke() -> None:
    assert __version__ == "0.1.0"
    starts = np.array([10, 100, 50], dtype=np.int32)
    ends = np.array([25, 140, 50], dtype=np.int32)
    assert fast_interval_sum(starts, ends) == 55

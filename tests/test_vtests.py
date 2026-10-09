import numpy as np
from dlrverify import vtests as vt


def test_kupiec_regions_table3():
    # Non-rejection regions at the 5 % level reported in Table 3 of the paper
    assert vt.kupiec_region(8760, 0.01) == (70, 106)
    assert vt.kupiec_region(365, 0.01) == (1, 7)
    assert vt.kupiec_region(2628, 0.01) == (17, 36)


def test_christoffersen_detects_clustering():
    I = np.zeros(5000, int)
    I[::100] = 1                                   # isolated hits
    assert vt.christoffersen(I, 0.01)["p_ind"] > 0.05
    J = np.zeros(5000, int)
    for s in range(0, 5000, 400):
        J[s:s + 4] = 1                             # runs of 4
    assert vt.christoffersen(J, 0.01)["p_ind"] < 0.001


def test_clopper_pearson_contains_rate():
    lo, hi = vt.clopper_pearson(14, 2038)
    assert lo < 14 / 2038 < hi


def test_runs():
    assert list(vt.runs(np.array([0, 1, 1, 0, 1, 0, 0, 1, 1, 1]))) == [2, 1, 3]

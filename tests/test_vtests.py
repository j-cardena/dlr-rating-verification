import numpy as np
from dlrverify import vtests as vt


def test_kupiec_regions_table4():
    # Non-rejection regions at the 5 % level reported in Table 4 of the paper
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


def test_hac_tests_reduce_to_classical_at_lag0():
    rng = np.random.default_rng(1)
    I = (rng.random(5000) < 0.012).astype(int)
    X = np.column_stack([np.ones(5000), np.r_[0, I[:-1]], (rng.random(5000) < 0.3).astype(float)])
    assert abs(vt.dq_hac(I, 0.01, X, 0)[1] - vt.dq_test(I, 0.01, X)[1]) < 1e-9
    # HAC coverage test: z-test with variance a(1-a) at lag 0
    p, pz, ucb = vt.uc_hac(I, 0.01, 0)
    assert 0 < pz < 1 and ucb > p


def test_hac_widens_for_overlapping_hits():
    # runs of 24 hits: the HAC coverage p-value must be larger than the iid one
    I = np.zeros(24000, int)
    for s in range(0, 24000, 1600):
        I[s:s + 24] = 1                      # rate 1.5 %
    p_iid = vt.uc_hac(I, 0.01, 0)[1]
    p_hac = vt.uc_hac(I, 0.01, 23)[1]
    assert p_hac > p_iid

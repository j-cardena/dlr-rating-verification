"""Coverage tests for a binary exceedance (hit) sequence I_t = 1{R_t > A_t}."""
import numpy as np
from scipy import stats, optimize


def _xlogy(n, p):
    return 0.0 if n == 0 else n * np.log(p)


def kupiec(I, a):
    I = np.asarray(I, int); T = I.size; x = I.sum(); p = x / T
    ll0 = _xlogy(T - x, 1 - a) + _xlogy(x, a)
    ll1 = _xlogy(T - x, 1 - p) + _xlogy(x, p) if 0 < x < T else 0.0
    lr = -2 * (ll0 - ll1)
    return lr, stats.chi2.sf(lr, 1)


def christoffersen(I, a):
    I = np.asarray(I, int)
    i0, i1 = I[:-1], I[1:]
    n00 = np.sum((i0 == 0) & (i1 == 0)); n01 = np.sum((i0 == 0) & (i1 == 1))
    n10 = np.sum((i0 == 1) & (i1 == 0)); n11 = np.sum((i0 == 1) & (i1 == 1))
    p01 = n01 / max(n00 + n01, 1); p11 = n11 / max(n10 + n11, 1)
    p2 = (n01 + n11) / (n00 + n01 + n10 + n11)
    l0 = _xlogy(n00 + n10, 1 - p2) + _xlogy(n01 + n11, p2)
    l1 = (_xlogy(n00, 1 - p01) + _xlogy(n01, p01)
          + _xlogy(n10, 1 - p11) + _xlogy(n11, p11))
    lr_ind = -2 * (l0 - l1)
    lr_uc, _ = kupiec(I, a)
    lr_cc = lr_uc + lr_ind
    return dict(lr_ind=lr_ind, p_ind=stats.chi2.sf(lr_ind, 1),
                lr_cc=lr_cc, p_cc=stats.chi2.sf(lr_cc, 2), pi01=p01, pi11=p11)


def duration_test(I):
    """Christoffersen-Pelletier (2004) Weibull duration test, H0: b = 1.
    Uses durations between hits with censored first/last spells."""
    I = np.asarray(I, int); idx = np.flatnonzero(I)
    if idx.size < 3:
        return np.nan, np.nan
    d = np.diff(idx).astype(float)
    d_first = float(idx[0] + 1); d_last = float(I.size - idx[-1])
    def nll(theta):
        a_, b_ = np.exp(theta)
        ll = np.sum(np.log(b_) + b_ * np.log(a_) + (b_ - 1) * np.log(d) - (a_ * d) ** b_)
        ll += -(a_ * d_first) ** b_ - (a_ * d_last) ** b_     # censored spells: log S(d)
        return -ll
    r1 = optimize.minimize(nll, [np.log(1 / d.mean()), 0.0], method="Nelder-Mead",
                           options=dict(xatol=1e-8, fatol=1e-10, maxiter=5000))
    r0 = optimize.minimize_scalar(lambda la: nll([la, 0.0]), bounds=(-15, 2), method="bounded")
    lr = 2 * (r0.fun - r1.fun)
    return np.exp(r1.x[1]), stats.chi2.sf(max(lr, 0.0), 1)


def dq_test(I, a, X):
    """Engle-Manganelli dynamic quantile test; X includes a constant column."""
    hit = np.asarray(I, float) - a
    X = np.asarray(X, float)
    beta, *_ = np.linalg.lstsq(X, hit, rcond=None)
    dq = beta @ X.T @ X @ beta / (a * (1 - a))
    return dq, stats.chi2.sf(dq, X.shape[1])


def clopper_pearson(x, n, conf=0.95):
    lo = stats.beta.ppf((1 - conf) / 2, x, n - x + 1) if x > 0 else 0.0
    hi = stats.beta.ppf(1 - (1 - conf) / 2, x + 1, n - x) if x < n else 1.0
    return lo, hi


def runs(I):
    I = np.asarray(I, int)
    d = np.diff(np.r_[0, I, 0]); starts = np.flatnonzero(d == 1); ends = np.flatnonzero(d == -1)
    return ends - starts


def kupiec_region(T, a, level=0.05):
    xs = np.arange(0, int(3 * a * T) + 20)
    ok = [x for x in xs if kupiec(np.r_[np.ones(x), np.zeros(T - x)], a)[1] >= level]
    return min(ok), max(ok)


# ---------------------------------------------------------------- dependence-robust versions
def newey_west(g, L):
    """Bartlett-kernel long-run covariance of the rows of g (T x k), bandwidth L."""
    g = np.atleast_2d(np.asarray(g, float).T).T
    g = g - g.mean(0)
    T = g.shape[0]; S = g.T @ g / T
    for k in range(1, L + 1):
        G = g[k:].T @ g[:-k] / T
        S = S + (1 - k / (L + 1)) * (G + G.T)
    return S


def uc_hac(I, a, L):
    """Unconditional coverage with a HAC variance for overlapping h-step forecasts
    (bandwidth L >= h-1; lag-0 variance a(1-a) imposed under H0, autocovariances at lags
    1..L from the data with Bartlett weights). Returns (rate, two-sided p for H0: rate = a,
    one-sided 95% upper bound on the rate, using the data-based long-run variance)."""
    I = np.asarray(I, float); T = I.size; p = I.mean(); H = I - a
    v0 = a * (1 - a)
    for k in range(1, L + 1):
        v0 += 2 * (1 - k / (L + 1)) * (H[k:] @ H[:-k]) / T
    z = (p - a) / np.sqrt(max(v0, 1e-12) / T)
    v1 = float(newey_west(I[:, None], L)[0, 0])
    return p, 2 * stats.norm.sf(abs(z)), p + 1.645 * np.sqrt(max(v1, p * (1 - p)) / T)


def dq_hac(I, a, X, L):
    """Dynamic quantile test as a moment test E[X_t (I_t - a)] = 0. The lag-0 covariance is
    imposed under H0, a(1-a) X'X/T, as in Engle and Manganelli (2004); for overlapping
    h-step forecasts the autocovariances at lags 1..L (L >= h-1) are added with Bartlett
    weights (Newey-West). With L = 0 this is the Engle-Manganelli statistic."""
    X = np.asarray(X, float); H = np.asarray(I, float) - a
    g = X * H[:, None]; T = g.shape[0]; m = g.mean(0)
    S = a * (1 - a) * X.T @ X / T
    for k in range(1, L + 1):
        G = g[k:].T @ g[:-k] / T
        S = S + (1 - k / (L + 1)) * (G + G.T)
    W = T * m @ np.linalg.pinv(S) @ m
    return W, stats.chi2.sf(W, g.shape[1])

"""Monte Carlo size of the coverage tests at the sample sizes of the worked example (v1.2).

Under H0 the hit sequence has rate alpha and, for h-step forecasts issued hourly, the
dependence of overlapping forecast errors: I_t = 1{z_t < Phi^-1(alpha)}, z_t a normalised
MA(h-1) of iid shocks (iid for h = 1). Reported: rejection rate of the unconditional
coverage test (Kupiec at h = 1, HAC at h = 24) and of the HAC dynamic quantile test at the
nominal 5 % level, and the coverage of the one-sided 95 % upper bound. A second design adds
a persistent AR(1) weather component (rho = 0.95 per hour, half the variance), which keeps the
marginal rate at alpha but makes the hits dependent beyond the bandwidth; this checks the
UC test and the bound only (the DQ null is then false by construction).
DQ regressors as in the paper: constant, hits lagged h and 2h, a persistent binary regime
indicator observed at issue time (frequency 0.2), and a daylight indicator.
Usage: size_sim.py [reps]"""
import sys
import numpy as np
from scipy import stats
from dlrverify import vtests as vt

REPS = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
rng = np.random.default_rng(20261013)
CASES = [(1, 5933), (1, 7344), (24, 5910), (24, 7321)]


def latent(T, h, persist):
    e = rng.standard_normal(T + h)
    z = np.convolve(e, np.ones(h), "valid")[:T] / np.sqrt(h)
    if persist:
        u = np.empty(T); u[0] = rng.standard_normal(); r = 0.95
        s = rng.standard_normal(T) * np.sqrt(1 - r ** 2)
        for t in range(1, T):
            u[t] = r * u[t - 1] + s[t]
        z = np.sqrt(0.5) * z + np.sqrt(0.5) * u
    return z


def regime(T):
    x = np.empty(T, int); x[0] = rng.random() < 0.2
    for t in range(1, T):                     # Markov chain, stationary frequency 0.2
        x[t] = (rng.random() < 0.9) if x[t - 1] else (rng.random() < 0.025)
    return x


def one(T, h, a, persist):
    I = (latent(T, h, persist) < stats.norm.ppf(a)).astype(int)
    L = h - 1
    p_uc = vt.kupiec(I, a)[1] if h == 1 else vt.uc_hac(I, a, L)[1]
    ub, _ = vt.upper_bound_eff(I)          # Bartlett bandwidth 95 h, as in the paper
    lag = lambda k: np.r_[np.zeros(k), I[:-k]]
    day = (np.arange(T) % 24 >= 6) & (np.arange(T) % 24 < 18)
    X = np.column_stack([np.ones(T), lag(h), lag(2 * h), regime(T), day])
    p_dq = vt.dq_hac(I, a, X, L)[1]
    return p_uc < 0.05, p_dq < 0.05, ub >= a


print(f"Monte Carlo size of coverage tests, {REPS} replications, nominal level 5 %")
print(f"{'design':10s} {'h':>3s} {'T':>5s} {'alpha':>5s}  {'UC size':>7s}  {'DQ size':>7s}  {'UB cover':>8s}")
for persist in [False, True]:
    for h, T in CASES:
        for a in [0.01, 0.05]:
            r = np.array([one(T, h, a, persist) for _ in range(REPS)])
            se = np.sqrt(0.05 * 0.95 / REPS) * 100
            dq = f"{r[:, 1].mean()*100:6.1f}%" if not persist else "   n/a "
            print(f"{'persist' if persist else 'overlap':10s} {h:3d} {T:5d} {a:5.2f}  {r[:, 0].mean()*100:6.1f}%  "
                  f"{dq}  {r[:, 2].mean()*100:7.1f}%", flush=True)
print(f"(Monte Carlo standard error of a 5 % rejection rate: {se:.1f} points)")

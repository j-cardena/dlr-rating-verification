"""How much evidence a one-sided safety test needs (v1.2.1).

For a nominal level alpha and a true exceedance probability pi < alpha, find the number n of
independent forecasts for which a one-sided 95 % Clopper-Pearson upper bound on the
exceedance frequency falls below alpha with probability 0.8 (power), and the corresponding
number of years of hourly ratings when the effective sample is T/D, with D the design effect
(1 for independent hourly hits; about 10 at h = 24 and alpha = 0.01 for equally weighted
overlap, Section 4.4). Exact binomial calculation.
Usage: power_analysis.py"""
import numpy as np
from scipy import stats

CONF, POWER = 0.95, 0.80


def ub(x, n):
    return stats.beta.ppf(CONF, x + 1, n - x) if x < n else 1.0


def power(n, pi, a):
    xs = np.arange(0, int(n * a) + 2)                  # bound < a needs x well below n*a
    ok = xs[[ub(x, n) < a for x in xs]]
    return 0.0 if ok.size == 0 else stats.binom.cdf(ok.max(), n, pi)


def n_required(pi, a, nmax=5_000_000):
    lo, hi = 10, 10
    while power(hi, pi, a) < POWER:
        hi *= 2
        if hi > nmax:
            return np.nan
    while hi - lo > max(1, hi // 200):                 # 0.5 % resolution; power is saw-toothed
        mid = (lo + hi) // 2
        lo, hi = (lo, mid) if power(mid, pi, a) >= POWER else (mid, hi)
    return hi


print(f"One-sided {CONF:.0%} Clopper-Pearson bound below alpha with power {POWER:.0%}")
print(f"{'alpha':>6} {'pi/alpha':>8} {'n_eff':>9} {'years, D=1':>11} {'years, D=10':>12}")
for a in (0.01, 0.02, 0.05):
    for r in (0.5, 0.8, 0.9):
        n = n_required(r * a, a)
        print(f"{a:6.2f} {r:8.1f} {n:9.0f} {n / 8760:11.2f} {10 * n / 8760:12.1f}", flush=True)

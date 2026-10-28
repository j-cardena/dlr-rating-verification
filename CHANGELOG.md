# Changelog

## 1.3.0 (2026-10-28)
- Version archived with the revised manuscript (Zenodo 10.5281/zenodo.23732605).
- Causal split fixed at h > 1: each monthly refit now uses only pairs whose outcome was known at
  the issue time of the month's first forecast (target time <= that issue time). v1.2.1 trained
  on all pairs with targets in earlier months, so the first h - 1 forecasts of each month used
  outcomes realised after issue: 230 of 7344 forecasts at Jena and 207 of 5910 at M2, all at
  h = 24 (none at h = 1). Counts in results.json (`v121_lookahead_forecasts`).
- Dependence calculations on elapsed hours: `vtests.newey_west`, `uc_hac`, `dq_hac`,
  `upper_bound_eff` and `christoffersen` take optional grid positions `pos`; pairs that straddle
  a data gap (M2, February-May 2021) are skipped. Without `pos` the results are unchanged.
- `upper_bound_eff` with no hits: effective sample T / (1 + L), the largest design effect of the
  Bartlett estimator at bandwidth L, instead of the independent bound with the nominal T.
  No reported rule has zero hits.
- Effect on the paper: h = 24 rows of Tables 10-11 at Jena change by up to 0.2 percentage points
  for the wind-informed rules and by 0.9 for the probabilistic AAR at alpha = 0.05
  (10.09 -> 9.23 %); at M2 frequencies change by at most 0.03 points and bounds in the second decimal.
  Table 13 selects the same rule in all twelve cases. Tables 6, 9, 12 and S2 unchanged
  (manuscript numbering; files table7/8/10 and table4/6/9).
- expected/ refreshed.

## 1.2.1 (2026-10)
- Quantile (pinball) score at alpha for every rule (`qscore` in results.json).
- Table 10: rule selected by lowest quantile score, by coverage closest to alpha, and by a
  safety-first criterion (largest gain with one-sided 95 % bound <= alpha). Pre-specified.
- New `scripts/power_analysis.py`: forecasts needed for a one-sided bound below alpha.
- Quantile-tracking step-size sensitivity (eta = 5 and 20 A) in results.json.
- `literature/screening_merged.csv`: new "Evidence basis" column (full text / abstract).
- No change to any number in Tables 4 and 6-9.

## 1.2.0 (2026-10)
- Seasonal ratings without look-ahead: design temperatures from before the test period (M2: 2020;
  Jena: MPI-BGC roof station 2017-2019, files `data/raw/mpi_roof_201[789]*.csv`). Jena seasonal
  exceedance 0.16 % -> 0.11 %; M2 unchanged (2020 maxima equal the full-record maxima).
- expected/ refreshed for v1.2.
- Monte Carlo p-values (Dufour 2006; 999 null samples, overlap of h-1 hours) for the h = 24
  coverage and dynamic quantile tests; asymptotic HAC values kept in results.json (p_uc_hac, p_dq_hac).
- One-sided 95 % upper bound on the exceedance rate: Clopper-Pearson with a long-run-variance
  effective sample (Bartlett, 95 h); new `vtests.upper_bound_eff`.
- New `size_sim.py`: size of the coverage tests and coverage of the bound (Supplementary Table S2).
- `diag_subhourly.py` covers both sites and writes subhourly.json; new Table 9 (reference definition).
- Tables 7-8 report the bound; Table 8 adds the dynamic quantile test; AAR rows labelled "realised T".

## 1.1.0 (2026-10)
- Causal evaluation (models refitted monthly on all earlier data) replaces leave-one-month-out.
- New reference ratings: seasonal (four seasons) and AAR with a 10 F margin in 5 F steps; wind-assumption sensitivity.
- New rules: input quantiles with fixed 45 degree wind angle; probabilistic AAR without wind data.
- "Online correction" renamed quantile tracking; hits applied when known (h hours after issue, on the time axis).
- HAC (Newey-West) unconditional coverage with one-sided upper bound, and HAC dynamic quantile test, for overlapping forecasts; lag-0 variance imposed under the null.
- Dynamic quantile regressors: hits lagged by h and 2h, low wind at issue, daylight at target.
- Gain in the hottest 10 % of hours; regime bins in results.json.
- New scripts: hac_bandwidth.py, diag_subhourly.py; diag_input.py removed (diagnostic now in run_example.py as InputRealDir).
- Tables 6-8 (paper v1.1) replace Tables 6-7.

## 1.0.0 (2026-10-12)
- Version archived with the submitted manuscript draft (Zenodo 10.5281/zenodo.23324915).

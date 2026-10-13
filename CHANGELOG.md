# Changelog

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

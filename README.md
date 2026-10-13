# dlr-rating-verification

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23324914.svg)](https://doi.org/10.5281/zenodo.23324914)

Reproduction code and data provenance for

> J. L. Cárdenas-Barrera, "From Predictive Distributions to Operational Line Ratings: A Review of
> Uncertainty Quantification, Verification, and Risk-Based Use of Ambient-Adjusted and Dynamic
> Ratings," *Forecasting* (MDPI), 2026 (under review).

The repository reproduces every computed result in the paper: the Kupiec non-rejection regions
(Table 4), the clustering simulation (Section 4.4), and the two-site worked example (Section 5:
Tables 6–8, Figure 1, the bandwidth and sub-hourly sensitivity checks). It also contains the
literature search strings, the screening of 176 Scopus records, and the full-text extraction sheet
(Supplementary Table S1).

**Version 1.1** replaces the leave-one-month-out evaluation of v1.0 with a causal design, adds
seasonal and margined ambient-adjusted reference ratings, a probabilistic AAR without wind data and
an input-quantile variant with a fixed wind angle, renames the online correction to quantile
tracking, and uses HAC (Newey–West) coverage and dynamic quantile tests for the overlapping
24-hour forecasts. See `CHANGELOG.md`.

## What is computed

For each site, hourly reference ampacity $A_t$ of a single ACSR Drake span is computed from
measured weather with the IEEE Std 738-2012 steady-state model. Rating rules are back-tested at
$\alpha = 0.01$ and $0.05$ for 1 h and 24 h lead times in a causal design: each month is forecast
with models fitted on all earlier data (Jena: test March–December 2020; M2: test 2021).

| Rule | Paper section | Description |
|---|---|---|
| QR | 3.1 | Linear quantile regression of $A_{t+h}$; rating $= \hat q_\alpha$ |
| QR + tracking | 4.6 | Online quantile tracking: offset $\theta \leftarrow \theta - \eta (I-\alpha)$, each hit applied $h$ hours after issue |
| Median − margin | 3.2 | QR median less the $\alpha$-quantile of training residuals |
| Input quantiles | 6.1 | Wind speed at $\alpha$, temperature and irradiance at $1-\alpha$; direction persisted or fixed at 45° |
| Probabilistic AAR | 5.1 | QR of $A_{t+h}$ on target-hour temperature, daylight and calendar harmonics (no wind data) |

Reference ratings: annual static, seasonal (four seasons), AAR with measured temperature, and AAR
with measured temperature + 10 °F rounded up to 5 °F steps. Each rating is evaluated with the tests
of Section 4: Kupiec coverage (h = 1), HAC coverage with a one-sided upper bound (h = 24),
Christoffersen independence, dynamic quantile test (hits lagged by h and 2h; HAC at h = 24),
regime bins, runs, magnitudes, mean gain over the static rating overall and in the hottest 10 % of
hours.

## Quick start

```bash
git clone <this repository> && cd dlr-rating-verification
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt && pip install -e .
make data        # downloads Jena; prints instructions for the M2 export; checks SHA-256
make all         # tests, simulations, both sites, Figure 1, Tables 6-8, comparison with expected/
```

Outputs go to `out/`:

| Output | Paper item |
|---|---|
| `out/table4_kupiec.txt` | Table 4 |
| `out/cluster_simulation.txt` | Section 4.4 (false-rejection rates 23/37/54 %) |
| `out/{jena,m2_ew,m2_ns}/results.json`, `log.txt` | Tables 6–8 and all numbers in Section 5 |
| `out/hac_bandwidth.txt` | Bandwidth sensitivity of the h = 24 tests (Section 5.2) |
| `out/jena/diag_subhourly.txt` | Hourly-mean vs 10-min reference ampacity (Section 5.3) |
| `out/fig_regime.pdf`, `.csv` | Figure 1 and its data |
| `out/table{6,7,8}_*_rows.tex` | LaTeX body rows of Tables 6–8 |

`make verify` compares the regenerated results with `expected/`. With the pinned versions in
`requirements.txt` the results are identical; other BLAS or statsmodels builds may move an
exceedance count by one or two, which the check tolerates while flagging any change in a test
decision.

## Data

Raw data are not redistributed; see [`data/README.md`](data/README.md).

- **Jena, Germany, 2020** — Max Planck Institute for Biogeochemistry weather station, 10-min,
  downloaded automatically (Darts library copy).
- **NREL NWTC M2 tower, Colorado, 2020–2021** — 1-min, exported manually from the NREL MIDC
  portal (instructions in `data/README.md`). Irradiance is missing February–May 2021; those hours
  are excluded.

## Layout

```
src/dlrverify/      ieee738.py (steady-state ampacity), vtests.py (coverage tests)
scripts/            run_example.py, fig_regime2.py, make_tables.py, hac_bandwidth.py,
                    diag_subhourly.py, kupiec_table.py, cluster_sim.py, verify_results.py,
                    get_data.sh
expected/           published results.json per site/orientation
literature/         search_queries.md, screening_merged.csv, table_S1_extraction.csv
tests/              unit tests (IEEE 738 annex example, Table 4 regions, HAC tests, test behaviour)
```

## Scope and limitations

The example demonstrates the verification protocol; it is not a forecasting benchmark. The
forecasts use no numerical weather prediction, the solar term uses measured horizontal
irradiance on the projected diameter instead of the standard's solar geometry, ratings are
steady state for a single span, and the reference ampacity is computed from hourly-mean station
weather, not measured on a conductor; `diag_subhourly.py` shows how much the definition of the
reference matters.

## Licence and citation

Code: MIT (see `LICENSE`). Data remain under the terms of their providers. Please cite the
article (see `CITATION.cff`).

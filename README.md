# dlr-rating-verification

Reproduction code and data provenance for

> J. L. Cárdenas-Barrera, "From Predictive Distributions to Operational Line Ratings: A Review of
> Uncertainty Quantification, Verification, and Risk-Based Use of Ambient-Adjusted and Dynamic
> Ratings," *Forecasting* (MDPI), 2026 (under review).

The repository reproduces every computed result in the paper: the Kupiec non-rejection regions
(Table 3), the clustering simulation (Section 4.4), and the two-site worked example (Section 5:
Tables 6–7, Figure 1 and the input-quantile diagnostic). It also contains the literature search
strings, the screening of 176 Scopus records, and the full-text extraction sheet (Supplementary
Table S1).

## What is computed

For each site, hourly reference ampacity $A_t$ of a single ACSR Drake span is computed from
measured weather with the IEEE Std 738-2012 steady-state model. Four rating rules are then
back-tested at $\alpha = 0.01$ and $0.05$ for 1 h and 24 h lead times with leave-one-month-out
validation:

| Rule | Paper section | Description |
|---|---|---|
| QR | 3.1 | Linear quantile regression of $A_{t+h}$; rating $= \hat q_\alpha$ |
| QR + online correction | 4.6 | Additive offset updated by $\theta \leftarrow \theta - \eta (I-\alpha)$, using hits known at issue time |
| Median − margin | 3.2 | QR median less the $\alpha$-quantile of training residuals |
| Input quantiles | 6.1 | Wind speed at $\alpha$, temperature and irradiance at $1-\alpha$, direction persisted |

plus a static rating and an ambient-adjusted rating (AAR, Order No. 881 structure, measured
temperature). Each rule is evaluated with the tests of Section 4: Kupiec unconditional coverage,
Christoffersen independence, Christoffersen–Pelletier duration, Engle–Manganelli dynamic quantile
(hit lagged by $h$), Clopper–Pearson intervals by regime, exceedance runs and magnitudes, and
capacity gain over the static rating.

## Quick start

```bash
git clone <this repository> && cd dlr-rating-verification
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt && pip install -e .
make data        # downloads Jena; prints instructions for the M2 export; checks SHA-256
make all         # tests, simulations, both sites, Figure 1, Tables 6-7, comparison with expected/
```

Run time is about 5 minutes on a laptop (most of it in the quantile regressions). Outputs go to `out/`:

| Output | Paper item |
|---|---|
| `out/table3_kupiec.txt` | Table 3 |
| `out/cluster_simulation.txt` | Section 4.4 (false-rejection rates 23/37/54 %) |
| `out/{jena,m2_ew,m2_ns}/results.json`, `log.txt` | Tables 6–7 and all numbers in Section 5 |
| `out/{site}/diag_input.txt` | Input-quantile rule with realised wind direction (Section 5.2) |
| `out/fig_regime.pdf`, `.csv` | Figure 1 and its data |
| `out/table6_jena_rows.tex`, `out/table7_m2_rows.tex` | LaTeX body rows of Tables 6–7 |

`make verify` compares the regenerated results with `expected/`, the files that produced the
published numbers. With the pinned versions in `requirements.txt` the results are identical;
other BLAS or statsmodels builds may move an exceedance count by one or two, which the check
tolerates while flagging any change in a test decision.

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
scripts/            run_example.py, diag_input.py, fig_regime2.py, make_tables.py,
                    kupiec_table.py, cluster_sim.py, verify_results.py, get_data.sh
expected/           published results.json per site/orientation
literature/         search_queries.md, screening_merged.csv, table_S1_extraction.csv
tests/              unit tests (IEEE 738 annex example, Table 3 regions, test behaviour)
```

## Scope and limitations

The example demonstrates the verification protocol; it is not a forecasting benchmark. The
forecasts use no numerical weather prediction, the solar term uses measured horizontal
irradiance on the projected diameter instead of the standard's solar geometry, ratings are
steady state for a single span, and the reference ampacity is computed from station weather,
not measured on a conductor.

## Licence and citation

Code: MIT (see `LICENSE`). Data remain under the terms of their providers. Please cite the
article (see `CITATION.cff`).

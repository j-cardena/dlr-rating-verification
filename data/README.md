# Input data

Raw data are **not** redistributed here. `scripts/get_data.sh` downloads the Jena file and
checks both files against `checksums.sha256` (SHA-256 of the exact files used for the paper).

## Jena, Germany (2020)
10-min measurements from the weather station of the Max Planck Institute for Biogeochemistry
(https://www.bgc-jena.mpg.de/wetter/), as redistributed in the Darts library repository
(`datasets/weather.csv`). Variables used: `T (degC)`, `wv (m/s)`, `wd (deg)`, `SWDR (W/m²)`.

## Jena, Germany (2017-2019): pre-test climatology for the seasonal ratings
Half-year archives of the same station (WS Beutenberg, roof, 10-min), `mpi_roof_2017a` to
`mpi_roof_2019b`, from https://www.bgc-jena.mpg.de/wetter/ (Data Download, left column; CC BY 4.0).
The archived January-June 2020 file is identical to the Darts copy in T, wind and SWDR. Unzip the six CSV files into `data/raw/` keeping their names;
the Makefile passes them to `run_example.py` automatically. Only air temperature is used (the
highest hourly mean of each season). Without them, seasons with no pre-test data fall back to the
2020 record and are flagged in `results.json` (`seasonal_insample`).

## NREL NWTC M2 tower, Boulder, Colorado (2020-2021)
Measurement and Instrumentation Data Center (MIDC), site NWTC:
https://midcdmz.nlr.gov/apps/daily.pl?site=NWTC

1. Start date 2020-01-01, end date 2021-12-31, "Entire Day".
2. Output type: "Selected 1-Min Data (ASCII Text)".
3. Select: Global Horizontal, Temperature (2m), Wind Speed (2m), Wind Speed (10m),
   Wind Direction (10m), Zenith Angle, Azimuth Angle.
4. Save as `data/raw/nwtc_m2_2020_2021.csv`.

The loader matches columns by header text, so column order does not matter. Values below -5000
are MIDC missing-value flags. Global horizontal irradiance is missing from February to May 2021;
those hours are excluded (14,740 complete hours remain).

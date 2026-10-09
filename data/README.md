# Input data

Raw data are **not** redistributed here. `scripts/get_data.sh` downloads the Jena file and
checks both files against `checksums.sha256` (SHA-256 of the exact files used for the paper).

## Jena, Germany (2020)
10-min measurements from the weather station of the Max Planck Institute for Biogeochemistry
(https://www.bgc-jena.mpg.de/wetter/), as redistributed in the Darts library repository
(`datasets/weather.csv`). Variables used: `T (degC)`, `wv (m/s)`, `wd (deg)`, `SWDR (W/m²)`.

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

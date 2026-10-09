#!/usr/bin/env bash
# Fetch the public input data. Raw data are not redistributed in this repository.
set -euo pipefail
mkdir -p data/raw
# 1) Jena 2020 (Max Planck Institute for Biogeochemistry), 10-min, as distributed with Darts
if [ ! -f data/raw/jena_2020_weather.csv ]; then
  curl -fsSL -o data/raw/jena_2020_weather.csv \
    https://raw.githubusercontent.com/unit8co/darts/master/datasets/weather.csv
fi
# 2) NREL NWTC M2 tower, 2020-01-01 to 2021-12-31, 1-min: manual export (see data/README.md)
if [ ! -f data/raw/nwtc_m2_2020_2021.csv ]; then
  echo "M2 data missing: export it from MIDC as described in data/README.md and save as"
  echo "  data/raw/nwtc_m2_2020_2021.csv"
fi
# Verify against the files used for the paper (a mismatch means upstream data changed)
sha256sum -c data/checksums.sha256 || echo "WARNING: checksum mismatch; results may differ from the paper"

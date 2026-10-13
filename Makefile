# Reproduce all computed results in the paper.  `make all` after `make data`.
PY ?= python3
export PYTHONPATH := src:scripts
OUT := out
# Pre-test climatology for the Jena seasonal ratings (MPI-BGC roof station, 2017-2019),
# comma-separated; see data/README.md. Empty: seasonal maxima from the pre-test record only.
JENA_CLIM ?= $(shell ls data/raw/mpi_roof_201[789]*.csv 2>/dev/null | paste -sd, -)

.PHONY: all data test sims jena m2 figures tables verify verify-jena clean

all: test sims jena m2 figures tables verify

data:
	bash scripts/get_data.sh

test:
	$(PY) -m pytest -q tests

sims: | $(OUT)
	$(PY) scripts/kupiec_table.py   > $(OUT)/table4_kupiec.txt     # Table 4
	$(PY) scripts/cluster_sim.py    > $(OUT)/cluster_simulation.txt # Section 4.4 simulation
	$(PY) scripts/size_sim.py 2000  > $(OUT)/size_simulation.txt    # size of the coverage tests (Table S1)

jena: | $(OUT)
	mkdir -p $(OUT)/jena
	$(PY) scripts/run_example.py jena data/raw/jena_2020_weather.csv $(OUT)/jena 90 $(JENA_CLIM) > $(OUT)/jena/log.txt
	$(PY) scripts/diag_subhourly.py jena data/raw/jena_2020_weather.csv $(OUT)/jena/series_h1_a01.csv 90 > $(OUT)/jena/diag_subhourly.txt

m2: jena | $(OUT)
	mkdir -p $(OUT)/m2_ew $(OUT)/m2_ns
	$(PY) scripts/run_example.py m2 data/raw/nwtc_m2_2020_2021.csv $(OUT)/m2_ew 90 > $(OUT)/m2_ew/log.txt
	$(PY) scripts/run_example.py m2 data/raw/nwtc_m2_2020_2021.csv $(OUT)/m2_ns 0  > $(OUT)/m2_ns/log.txt
	$(PY) scripts/diag_subhourly.py m2 data/raw/nwtc_m2_2020_2021.csv $(OUT)/m2_ew/series_h1_a01.csv 90 > $(OUT)/m2_ew/diag_subhourly.txt
	$(PY) scripts/diag_subhourly.py m2 data/raw/nwtc_m2_2020_2021.csv $(OUT)/m2_ns/series_h1_a01.csv 0  > $(OUT)/m2_ns/diag_subhourly.txt
	$(PY) scripts/hac_bandwidth.py $(OUT)/jena/series_h24_a01.csv $(OUT)/m2_ew/series_h24_a01.csv \
	      $(OUT)/m2_ns/series_h24_a01.csv > $(OUT)/hac_bandwidth.txt

figures:
	$(PY) scripts/fig_regime2.py $(OUT)/jena/series_h1_a01.csv "Jena" \
	      $(OUT)/m2_ns/series_h1_a01.csv "M2, N–S line" $(OUT)/fig_regime.pdf   # Figure 1

tables:
	$(PY) scripts/make_tables.py $(OUT)/jena $(OUT)/m2_ew $(OUT)/m2_ns $(OUT)   # Tables 6-9, S1

verify:
	$(PY) scripts/verify_results.py $(OUT)/jena  expected/jena
	$(PY) scripts/verify_results.py $(OUT)/m2_ew expected/m2_ew
	$(PY) scripts/verify_results.py $(OUT)/m2_ns expected/m2_ns

verify-jena:
	$(PY) scripts/verify_results.py $(OUT)/jena  expected/jena

$(OUT):
	mkdir -p $(OUT)

clean:
	rm -rf $(OUT)

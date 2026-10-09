# Reproduce all computed results in the paper.  `make all` after `make data`.
PY ?= python3
export PYTHONPATH := src:scripts
OUT := out

.PHONY: all data test sims jena m2 figures tables verify verify-jena clean

all: test sims jena m2 figures tables verify

data:
	bash scripts/get_data.sh

test:
	$(PY) -m pytest -q tests

sims: | $(OUT)
	$(PY) scripts/kupiec_table.py   > $(OUT)/table3_kupiec.txt     # Table 3
	$(PY) scripts/cluster_sim.py    > $(OUT)/cluster_simulation.txt # Section 4.4 simulation

jena: | $(OUT)
	mkdir -p $(OUT)/jena
	$(PY) scripts/run_example.py jena data/raw/jena_2020_weather.csv $(OUT)/jena 90 > $(OUT)/jena/log.txt
	$(PY) scripts/diag_input.py  jena data/raw/jena_2020_weather.csv $(OUT)/jena 90 > $(OUT)/jena/diag_input.txt

m2: | $(OUT)
	mkdir -p $(OUT)/m2_ew $(OUT)/m2_ns
	$(PY) scripts/run_example.py m2 data/raw/nwtc_m2_2020_2021.csv $(OUT)/m2_ew 90 > $(OUT)/m2_ew/log.txt
	$(PY) scripts/run_example.py m2 data/raw/nwtc_m2_2020_2021.csv $(OUT)/m2_ns 0  > $(OUT)/m2_ns/log.txt
	$(PY) scripts/diag_input.py  m2 data/raw/nwtc_m2_2020_2021.csv $(OUT)/m2_ew 90 > $(OUT)/m2_ew/diag_input.txt
	$(PY) scripts/diag_input.py  m2 data/raw/nwtc_m2_2020_2021.csv $(OUT)/m2_ns 0  > $(OUT)/m2_ns/diag_input.txt

figures:
	$(PY) scripts/fig_regime2.py $(OUT)/jena/series_h1_a01.csv "Jena" \
	      $(OUT)/m2_ns/series_h1_a01.csv "M2, N–S line" $(OUT)/fig_regime.pdf   # Figure 1

tables:
	$(PY) scripts/make_tables.py $(OUT)/jena $(OUT)/m2_ew $(OUT)/m2_ns $(OUT)   # Tables 6-7

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

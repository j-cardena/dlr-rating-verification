# Literature search (Section 2 of the paper)

Database: Scopus, title-abstract-keywords. Both searches run on 8 October 2026.

**Q1 — risk, decision rules and verification (67 records)**
```
TITLE-ABS-KEY(("dynamic line rating" OR "dynamic thermal rating" OR "ambient adjusted rating*") AND (probabilistic OR quantile OR "chance constrain*" OR risk) AND (verification OR reliability OR calibration OR backtest* OR "decision rule*" OR "safe rating*")) AND PUBYEAR > 2014
```

**Q2 — probabilistic rating forecasting (129 records)**
```
TITLE-ABS-KEY(("dynamic line rating" OR "dynamic thermal rating" OR "ampacity" OR "ambient adjusted rating*") AND (forecast* OR predict*) AND (probabilistic OR quantile OR interval OR ensemble OR conformal OR uncertainty)) AND PUBYEAR > 2009
```

After removing 20 duplicates: 176 unique records (`screening_merged.csv`), tiered on title and abstract:

| Tier | Q1 | Q2 | Total |
|---|---|---|---|
| 1A Core: rating rule, verification or operations | 8 | 20 | 28 |
| 1B Core: probabilistic rating forecasting | 5 | 19 | 24 |
| 2 Support | 24 | 49 | 73 |
| 3 Peripheral | 18 | 0 | 18 |
| 4 Excluded | 12 | 21 | 33 |

A record found by both searches is listed once, under Q1. Two Q2 records are proceedings
front matter without titles ("Conference review") and are excluded.

`table_S1_extraction.csv` is Supplementary Table S1: full-text extraction for the 16 studies read in full.
Coverage was cross-checked with LeapSpace (37 references) and SciSpace (23 entries); every reference
taken from either tool was verified against the publisher's record.

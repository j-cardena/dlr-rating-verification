"""Sensitivity of the h = 24 coverage and dynamic quantile p-values to the Newey-West bandwidth
(Section 5.2). Usage: hac_bandwidth.py <series_h24_a01.csv> [<more series> ...]"""
import sys, numpy as np, pandas as pd
from dlrverify import vtests as vt
a, h = 0.01, 24
for f in sys.argv[1:]:
    d = pd.read_csv(f, index_col=0, parse_dates=True)
    tt = d.index + pd.Timedelta(hours=h)
    for k in ["QR", "QR+tracking", "Margin", "Input", "Input45", "ProbAAR"]:
        I = (d[k] > d["y"]).astype(int).values
        Is = pd.Series(I, index=tt)
        lag = lambda q: Is.reindex(tt - pd.Timedelta(hours=q)).fillna(0).values
        X = np.column_stack([np.ones(len(I)), lag(h), lag(2 * h), (d["V0"].values < 1).astype(float), d["y_day"].values])
        cells = [f"L={L}: pUC {vt.uc_hac(I, a, L)[1]:.3f} pDQ {vt.dq_hac(I, a, X, L)[1]:.3f}" for L in (23, 47, 95)]
        print(f"{f} {k:12s} " + " | ".join(cells))

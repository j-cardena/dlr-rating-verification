"""Diagnostic: input-quantile rule with persisted vs realised wind direction."""
import sys, numpy as np, pandas as pd
sys.argv = [sys.argv[0]] + sys.argv[1:]
import run_example as rx
from dlrverify import vtests as vt
for h in [1, 24]:
    df = rx.design(h); a = 0.01
    res = {"persist": [], "realised": [], "realisedV": []}
    for mth in sorted(df.month.unique()):
        tr, te = df[df.month != mth], df[df.month == mth]
        Vq, _ = rx.qfit(tr, te, "y_V", rx.XIN["V"], a)
        Tq, _ = rx.qfit(tr, te, "y_T", rx.XIN["T"], 1 - a)
        Gq, _ = rx.qfit(tr, te, "y_G", rx.XIN["G"], 1 - a)
        for k, phi, V in [("persist", te.phi0, Vq), ("realised", te.y_phi, Vq), ("realisedV", te.phi0, te.y_V)]:
            R = rx.ampacity(Tq, np.maximum(V, 0), phi, np.clip(Gq, 0, None), Tmax=rx.TMAX, elev=rx.ELEV)
            res[k].append((R > te.y.values).astype(int))
    for k, v in res.items():
        I = np.concatenate(v); print(h, k, round(I.mean() * 100, 2), round(vt.kupiec(I, a)[1], 3))

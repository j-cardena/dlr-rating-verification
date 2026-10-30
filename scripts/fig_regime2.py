"""Figure: exceedance frequency by wind regime at two sites (h=1, alpha=0.01).
Usage: fig_regime2.py <series site A> <label A> <series site B> <label B> <out.pdf>"""
import sys, pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from dlrverify import vtests as vt
fa, la, fb, lb, out = sys.argv[1:6]
bins = [0, 1, 2, 3, 99]; labels = ["<1", "1–2", "2–3", ">3"]
meths = [("QR", "QR", "o"), ("QR+tracking", "QR + quantile tracking", "s"),
         ("Margin", "Median − margin", "^"), ("Input", "Input quantiles", "D")]
fig, axs = plt.subplots(2, 2, figsize=(7.2, 5.6), sharey=True)
rows = []
for r, (f, site) in enumerate([(fa, la), (fb, lb)]):
    s = pd.read_csv(f, index_col=0, parse_dates=True)
    for c, (col, title) in enumerate([("V0", "regime at issue time"), ("y_V", "realized regime at target")]):
        ax = axs[r, c]; cat = pd.cut(s[col], bins, labels=labels, right=False)
        for j, (m, lab, mk) in enumerate(meths):
            I = (s[m] > s["y"]).astype(int).values
            rr, lo, hi = [], [], []
            for b in labels:
                msk = (cat == b).values; n = msk.sum(); k = I[msk].sum()
                a_, b_ = vt.clopper_pearson(k, n); rr.append(100 * k / n); lo.append(100 * a_); hi.append(100 * b_)
                rows.append((site, col, m, b, int(n), int(k), 100 * k / n, 100 * a_, 100 * b_))
            x = np.arange(4) + (j - 1.5) * 0.15
            ax.errorbar(x, rr, yerr=[np.array(rr) - lo, np.array(hi) - rr], fmt=mk, ms=4, capsize=2, label=lab, lw=0.8)
        ax.axhline(1.0, color="k", lw=0.8, ls="--")
        ax.set_xticks(range(4)); ax.set_xticklabels(labels)
        if r == 1: ax.set_xlabel("Hourly mean wind speed (m/s)")
        ax.set_title(f"({'abcd'[2*r+c]}) {site}: {title}", fontsize=8.5)
        ax.set_yscale("symlog", linthresh=1); ax.set_ylim(-0.08, 30)
        ax.set_yticks([0, 0.5, 1, 2, 5, 10, 20]); ax.set_yticklabels(["0", "0.5", "1", "2", "5", "10", "20"])
    axs[r, 0].set_ylabel("Exceedance frequency (%)")
axs[0, 0].legend(fontsize=7, loc="upper right", frameon=False)
fig.tight_layout(); fig.savefig(out, bbox_inches="tight")
pd.DataFrame(rows, columns=["site", "regime", "method", "bin", "n", "k", "rate", "lo", "hi"]).to_csv(out.replace(".pdf", ".csv"), index=False)

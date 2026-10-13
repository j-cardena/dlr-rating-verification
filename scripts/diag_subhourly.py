"""Sensitivity of the reference ampacity to hourly averaging (Jena, 10-min data).
Compares A computed from hourly-mean weather (the reference used in the paper) with the
mean and the minimum of the six 10-min ampacities in each hour, and re-scores the h = 1
ratings at alpha = 0.01 against those alternatives.
Usage: diag_subhourly.py <jena weather.csv> <series_h1_a01.csv>"""
import sys, numpy as np, pandas as pd
from dlrverify.ieee738 import ampacity
src, ser = sys.argv[1:3]
raw = pd.read_csv(src)
raw["t"] = pd.to_datetime(raw["Date Time"], format="%d.%m.%Y %H:%M:%S"); raw = raw.set_index("t")
V = raw["wv (m/s)"].where(raw["wv (m/s)"] >= 0)
d = np.abs(((raw["wd (deg)"] - 90.0) + 180) % 360 - 180); phi = np.where(d > 90, 180 - d, d)
A10 = pd.Series(ampacity(raw["T (degC)"], V, phi, raw["SWDR (W/m²)"], Tmax=75.0, elev=155.0), index=raw.index)
g = A10.resample("1h", label="right", closed="right")
alt = pd.DataFrame({"mean10": g.mean(), "min10": g.min()})
s = pd.read_csv(ser, index_col=0, parse_dates=True)
tgt = s.index + pd.Timedelta(hours=1)
m = alt.reindex(tgt).values
ok = ~np.isnan(m).any(1)
y = s["y"].values[ok]
print(f"hours compared: {ok.sum()}")
print(f"A(hourly mean weather) / mean of 10-min A: median {np.median(y / m[ok, 0]):.3f}")
print(f"A(hourly mean weather) / min of 10-min A:  median {np.median(y / m[ok, 1]):.3f}, 95th pct {np.percentile(y / m[ok, 1], 95):.3f}")
for k in ["QR", "QR+tracking", "Margin", "AAR881"]:
    R = s[k].values[ok]
    print(f"{k:12s} exceedance vs hourly-mean reference {100*np.mean(R > y):5.2f}%  "
          f"vs mean of 10-min A {100*np.mean(R > m[ok, 0]):5.2f}%  vs min of 10-min A {100*np.mean(R > m[ok, 1]):5.2f}%")

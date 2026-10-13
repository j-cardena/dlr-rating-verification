"""Sensitivity of the reference ampacity to its definition within the hour (v1.2).
Compares A computed from hourly-mean weather (the reference used in the paper) with the
mean and the minimum of the 10-min ampacities in each hour, and re-scores the h = 1 ratings
at alpha = 0.01 against those alternatives. Jena: 10-min station data. M2: 1-min tower data
averaged to 10 min (wind as a vector mean).
Usage: diag_subhourly.py <jena|m2> <raw data file> <series_h1_a01.csv> [azimuth, default 90]
Writes a JSON summary next to the series file (subhourly.json) and prints a table."""
import sys, os, json, numpy as np, pandas as pd
from dlrverify.ieee738 import ampacity
site, src, ser = sys.argv[1:4]
AZ = float(sys.argv[4]) if len(sys.argv) > 4 else 90.0
ELEV = {"jena": 155.0, "m2": 1855.0}[site]
raw = pd.read_csv(src)
if site == "jena":
    raw["t"] = pd.to_datetime(raw["Date Time"], format="%d.%m.%Y %H:%M:%S"); raw = raw.set_index("t")
    T, V, WD, G = raw["T (degC)"], raw["wv (m/s)"].where(raw["wv (m/s)"] >= 0), raw["wd (deg)"], raw["SWDR (W/m²)"]
else:
    pick = lambda key: next(c for c in raw.columns if key.lower() in c.lower())
    raw["t"] = pd.to_datetime(raw[pick("DATE")] + " " + raw[pick("MST")], format="%m/%d/%Y %H:%M")
    raw = raw.set_index("t")
    f = lambda c: raw[c].where(raw[c] > -5000)
    ws, wd = f(pick("Wind Speed @ 10m")), np.radians(f(pick("Wind Direction @ 10m")))
    m = pd.DataFrame({"T": f(pick("Temperature @ 2m")), "V": ws.where(ws >= 0),
                      "u": ws * np.sin(wd), "v": ws * np.cos(wd),
                      "G": f(pick("Global Horizontal")).clip(lower=0)}).resample(
        "10min", label="right", closed="right").mean()
    T, V, G = m["T"], m["V"], m["G"]
    WD = (np.degrees(np.arctan2(m["u"], m["v"])) + 360) % 360
d = np.abs(((WD - AZ) + 180) % 360 - 180); phi = np.where(d > 90, 180 - d, d)
A10 = pd.Series(ampacity(T, V, phi, G, Tmax=75.0, elev=ELEV), index=T.index)
g = A10.resample("1h", label="right", closed="right")
alt = pd.DataFrame({"mean10": g.mean(), "min10": g.min(), "n": g.count()})
s = pd.read_csv(ser, index_col=0, parse_dates=True)
tgt = s.index + pd.Timedelta(hours=1)
a = alt.reindex(tgt)
ok = (a["n"].values == 6) & ~a[["mean10", "min10"]].isna().any(axis=1).values
y = s["y"].values[ok]; m10, n10 = a["mean10"].values[ok], a["min10"].values[ok]
out = dict(hours=int(ok.sum()), ratio_mean_median=float(np.median(y / m10)),
           ratio_min_median=float(np.median(y / n10)), ratio_min_p95=float(np.percentile(y / n10, 95)), rules={})
print(f"hours compared: {ok.sum()}")
print(f"A(hourly mean weather) / mean of 10-min A: median {out['ratio_mean_median']:.3f}")
print(f"A(hourly mean weather) / min of 10-min A:  median {out['ratio_min_median']:.3f}, 95th pct {out['ratio_min_p95']:.3f}")
for k in ["QR", "QR+tracking", "Margin", "ProbAAR", "AAR881"]:
    R = s[k].values[ok]
    r = [float(np.mean(R > ref)) for ref in (y, m10, n10)]
    out["rules"][k] = dict(hourly=r[0], mean10=r[1], min10=r[2])
    print(f"{k:12s} exceedance vs hourly-mean reference {100*r[0]:5.2f}%  "
          f"vs mean of 10-min A {100*r[1]:5.2f}%  vs min of 10-min A {100*r[2]:5.2f}%")
json.dump(out, open(os.path.join(os.path.dirname(ser), "subhourly.json"), "w"), indent=1)

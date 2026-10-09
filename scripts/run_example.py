"""Worked example: coverage back-test of probabilistic line ratings.

Data: 10-min weather at the Max Planck Institute for Biogeochemistry station,
Jena (Germany), 2020, as distributed with the Darts library (weather.csv).
Line: hypothetical ACSR Drake, east-west axis, T_max = 75 C, elevation 155 m.
Reference ampacity A_t: IEEE 738-2012 steady state from hourly-mean measured weather.

Usage: python3 -I run_example.py <jena|m2> <data file> <outdir> [azimuth]
"""
import sys, json
import numpy as np, pandas as pd
import statsmodels.api as sm
from statsmodels.regression.quantile_regression import QuantReg
from dlrverify.ieee738 import ampacity
from dlrverify import vtests as vt

# Usage: run_example.py <site: jena|m2> <data file> <outdir> [line azimuth deg, default 90]
SITE, SRC, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
AZ = float(sys.argv[4]) if len(sys.argv) > 4 else 90.0
TMAX = 75.0
ELEV = {"jena": 155.0, "m2": 1855.0}[SITE]
ALPHAS = [0.01, 0.05]
HORIZONS = [1, 24]
ETA = 10.0  # A, step of the online correction (fixed a priori)


# ---------------------------------------------------------------- data

def _m2_columns(cols):
    """Map MIDC export headers to short names, independent of column order."""
    want = {"date": ["DATE"], "time": ["MST"], "G": ["Global Horizontal"],
            "T": ["Temperature @ 2m"], "V10": ["Wind Speed @ 10m"], "WD": ["Wind Direction @ 10m"]}
    out = {}
    for short, keys in want.items():
        hit = [c for c in cols if any(k.lower() in c.lower() for k in keys)]
        if not hit:
            raise ValueError(f"M2 file lacks a column for {short}: expected header containing {keys}")
        out[hit[0]] = short
    return out

def load(site, src):
    """Hour-ending means of T, V, wind vector (u, v) and irradiance G."""
    if site == "jena":       # MPI-BGC Jena 2020, 10-min (Darts weather.csv)
        raw = pd.read_csv(src)
        raw["t"] = pd.to_datetime(raw["Date Time"], format="%d.%m.%Y %H:%M:%S")
        raw = raw.set_index("t")
        T, ws, wdeg, G = raw["T (degC)"], raw["wv (m/s)"], raw["wd (deg)"], raw["SWDR (W/m²)"]
        span = ("2020-01-01 01:00", "2020-12-31 23:00")
    else:                    # NREL NWTC M2 tower, 1-min, MIDC export (MST)
        raw = pd.read_csv(src)
        raw = raw.rename(columns=_m2_columns(raw.columns))
        raw["t"] = pd.to_datetime(raw["date"] + " " + raw["time"], format="%m/%d/%Y %H:%M")
        raw = raw.set_index("t")
        for c in ["G", "T", "V10", "WD"]:
            raw[c] = raw[c].where(raw[c] > -5000)          # MIDC missing-value flags
        T, ws, wdeg, G = raw["T"], raw["V10"], raw["WD"], raw["G"].clip(lower=0)
        span = (str(raw.index.min().floor("h") + pd.Timedelta(hours=1)), str(raw.index.max().floor("h")))
    ws = ws.where(ws >= 0)
    wd = np.radians(wdeg)
    u, v = ws * np.sin(wd), ws * np.cos(wd)
    agg = pd.DataFrame({"T": T, "V": ws, "u": u, "v": v, "G": G}).resample(
        "1h", label="right", closed="right").mean()
    agg = agg.loc[span[0]:span[1]].interpolate(limit=3, limit_area="inside")
    # keep a complete hourly index: lags and leads must be time-aligned across data gaps;
    # rows with missing inputs are dropped later, after lags are formed
    return agg


agg = load(SITE, SRC)
agg["dir"] = (np.degrees(np.arctan2(agg.u, agg.v)) + 360) % 360
d = np.abs(((agg["dir"] - AZ) + 180) % 360 - 180)        # 0..180
agg["phi"] = np.where(d > 90, 180 - d, d)                    # angle to conductor axis
agg["A"] = ampacity(agg["T"], agg["V"], agg["phi"], agg["G"], Tmax=TMAX, elev=ELEV)
agg["day"] = (agg["G"] > 5).astype(float).where(agg["G"].notna())
hours = agg.index

# static and ambient-adjusted ratings (FERC 881-style: forecast temperature, day/night sun)
STATIC = float(ampacity(30.0, 0.61, 90.0, 1000.0, Tmax=TMAX, elev=ELEV))
agg["AAR"] = ampacity(agg["T"], 0.61, 90.0, np.where(agg["day"] == 1, 1000.0, 0.0),
                      Tmax=TMAX, elev=ELEV)   # perfect temperature foresight: lower bound on AAR risk


# ---------------------------------------------------------------- features
def design(h):
    df = pd.DataFrame(index=hours)
    s = agg
    df["y"] = s["A"].shift(-h)                       # target A_{t+h}
    df["A0"] = s["A"]; df["A1"] = s["A"].shift(1); df["A2"] = s["A"].shift(2)
    df["A24"] = s["A"].shift(24 - h) if h < 24 else s["A"]
    df["Am24"] = s["A"].rolling(24).mean()
    df["V0"] = s["V"]; df["T0"] = s["T"]; df["G0"] = s["G"]
    tgt = hours + pd.Timedelta(hours=h)
    df["sh"] = np.sin(2 * np.pi * tgt.hour / 24); df["ch"] = np.cos(2 * np.pi * tgt.hour / 24)
    for var in ["V", "T", "G"]:
        df["y_" + var] = s[var].shift(-h)
        df[var + "1"] = s[var].shift(1)
    df["y_phi"] = s["phi"].shift(-h); df["phi0"] = s["phi"]
    df["y_day"] = s["day"].shift(-h); df["y_AAR"] = s["AAR"].shift(-h)
    df["month"] = tgt.year * 100 + tgt.month
    return df.dropna()


XA = ["A0", "A1", "A2", "A24", "Am24", "V0", "T0", "G0", "sh", "ch"]
XIN = {"V": ["V0", "V1", "T0", "G0", "sh", "ch"],
       "T": ["T0", "T1", "G0", "sh", "ch"],
       "G": ["G0", "G1", "T0", "sh", "ch"]}


def qfit(df_tr, df_te, ycol, xcols, q):
    Xtr = sm.add_constant(df_tr[xcols].values, has_constant="add")
    Xte = sm.add_constant(df_te[xcols].values, has_constant="add")
    m = QuantReg(df_tr[ycol].values, Xtr).fit(q=q, max_iter=5000)
    return Xte @ m.params, Xtr @ m.params


def lomo(df, a):
    """Leave-one-month-out forecasts for every method at level a."""
    out = []
    for mth in sorted(df.month.unique()):
        tr, te = df[df.month != mth], df[df.month == mth]
        r = pd.DataFrame(index=te.index)
        # M1 direct quantile regression of ampacity
        r["QR"], _ = qfit(tr, te, "y", XA, a)
        # M3 point forecast (median) minus constant margin from training residuals
        med_te, med_tr = qfit(tr, te, "y", XA, 0.5)
        m = np.quantile(tr["y"].values - med_tr, a)
        r["Margin"] = med_te + m
        # M2 input-level confidence: V at a, T at 1-a, G at 1-a, direction persisted
        Vq, _ = qfit(tr, te, "y_V", XIN["V"], a)
        Tq, _ = qfit(tr, te, "y_T", XIN["T"], 1 - a)
        Gq, _ = qfit(tr, te, "y_G", XIN["G"], 1 - a)
        r["Input"] = ampacity(Tq, np.maximum(Vq, 0), te["phi0"], np.clip(Gq, 0, None),
                              Tmax=TMAX, elev=ELEV)
        out.append(r)
    return pd.concat(out).sort_index()


def online_correct(q, y, a, h, eta=ETA):
    """Quantile tracking on an additive offset; hits become known h steps later."""
    theta = np.zeros(len(q)); R = np.empty(len(q)); th = 0.0
    hits = np.zeros(len(q), int)
    for i in range(len(q)):
        if i - h >= 0:                      # information available at issue time
            th -= eta * (hits[i - h] - a)
        R[i] = q[i] + th; hits[i] = int(R[i] > y[i]); theta[i] = th
    return R


# ---------------------------------------------------------------- evaluate
def evaluate(R, df, a, h):
    y = df["y"].values; I = (R > y).astype(int); T = len(I); x = int(I.sum())
    lr_uc, p_uc = vt.kupiec(I, a)
    ch = vt.christoffersen(I, a)
    b, p_dur = vt.duration_test(I)
    lowV = (df["V0"].values < 1.0).astype(float)          # regime known at issue time
    lag = np.r_[np.zeros(h), I[:-h]]                        # hit known at issue time
    X = np.column_stack([np.ones(T), lag, lowV, df["y_day"].values])
    dq, p_dq = vt.dq_test(I, a, X)
    rl = vt.runs(I)
    reg = {}
    lowT = df["y_V"].values < 1.0
    for name, mask in [("lowwind", lowT), ("otherwind", ~lowT),
                       ("day", df["y_day"].values == 1), ("night", df["y_day"].values == 0)]:
        n = int(mask.sum()); k = int(I[mask].sum())
        lo, hi = vt.clopper_pearson(k, n)
        reg[name] = dict(n=n, k=k, rate=k / n, lo=lo, hi=hi)
    mag = (R[I == 1] - y[I == 1]) / y[I == 1] if x else np.array([np.nan])
    return dict(T=T, x=x, rate=x / T, p_uc=p_uc, p_ind=ch["p_ind"], p_cc=ch["p_cc"],
                pi01=ch["pi01"], pi11=ch["pi11"], weib_b=b, p_dur=p_dur, p_dq=p_dq,
                episodes=int(rl.size), mean_run=float(rl.mean()) if rl.size else 0.0,
                max_run=int(rl.max()) if rl.size else 0,
                mag_med=float(np.nanmedian(mag)), mag_p95=float(np.nanpercentile(mag, 95)),
                ratio_med=float(np.median(R / y)), mean_R=float(R.mean()),
                gain_vs_static=float(R.mean() / STATIC - 1), regimes=reg)


def main():
    results = []
    for h in HORIZONS:
        df = design(h)
        # benchmarks (independent of alpha)
        bench = {"Static": np.full(len(df), STATIC), "AAR": df["y_AAR"].values}
        for a in ALPHAS:
            f = lomo(df, a)
            f["QR+online"] = online_correct(f["QR"].values, df["y"].values, a, h)
            meths = dict(bench); meths.update({k: f[k].values for k in ["QR", "Input", "Margin", "QR+online"]})
            for name, R in meths.items():
                r = evaluate(np.asarray(R, float), df, a, h); r.update(method=name, h=h, alpha=a)
                results.append(r)
                print(f"h={h:2d} a={a:.2f} {name:10s} x={r['x']:4d}/{r['T']} rate={r['rate']*100:5.2f}% "
                      f"pUC={r['p_uc']:.3f} pIND={r['p_ind']:.3f} pDur={r['p_dur']:.3f} pDQ={r['p_dq']:.3f} "
                      f"ep={r['episodes']} mrun={r['mean_run']:.2f} max={r['max_run']} "
                      f"low={r['regimes']['lowwind']['rate']*100:.1f}% oth={r['regimes']['otherwind']['rate']*100:.2f}% "
                      f"mag={r['mag_med']*100:.1f}/{r['mag_p95']*100:.1f}% R/A={r['ratio_med']:.2f} "
                      f"gain={r['gain_vs_static']*100:+.1f}%", flush=True)
            if h == 1 and a == ALPHAS[0]:
                sav = df[["y", "V0", "y_V", "y_day"]].copy()
                for k in ["QR", "Input", "Margin", "QR+online"]:
                    sav[k] = f[k].values
                sav["AAR"] = bench["AAR"]; sav.to_csv(f"{OUT}/series_h1_a01.csv")
            if h == 24 and a == ALPHAS[0]:
                sav = df[["y", "V0", "y_V", "y_day"]].copy()
                for k in ["QR", "Input", "Margin", "QR+online"]:
                    sav[k] = f[k].values
                sav.to_csv(f"{OUT}/series_h24_a01.csv")

    ok = agg.dropna(subset=["A"])
    summary = dict(static=STATIC, n_hours=int(len(ok)), n_index=int(len(agg)),
                   A_quantiles={str(q): float(ok["A"].quantile(q)) for q in [0.01, 0.05, 0.5, 0.95]},
                   frac_lowwind=float((ok["V"] < 1).mean()),
                   static_exceed=float((STATIC > ok["A"]).mean()),
                   kupiec_region={f"{T}_{a}": vt.kupiec_region(T, a) for T, a in
                                  [(len(ok), 0.01), (len(ok), 0.05)]})
    json.dump(dict(summary=summary, results=results), open(f"{OUT}/results.json", "w"),
              indent=1, default=float)
    print(json.dumps(summary, indent=1, default=float))


if __name__ == "__main__":
    main()

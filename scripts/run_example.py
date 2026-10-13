"""Worked example: coverage back-test of probabilistic line ratings (v1.1: causal design).

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

# ---------------------------------------------------------------- reference ratings
# Annual static rating (0.61 m/s perpendicular, 30 C, 1000 W/m2), the seasonal ratings
# that Order No. 881 keeps beyond 10 days (four seasons, design temperature = highest hourly
# temperature of that season in the record), and two ambient-adjusted ratings with the same
# wind assumption: with the measured temperature (no forecast error, no margin) and with an
# Order No. 881-style margin (measured temperature + 10 F, rounded up to the next 5 F step).
STATIC = float(ampacity(30.0, 0.61, 90.0, 1000.0, Tmax=TMAX, elev=ELEV))
SUN = lambda day: np.where(day == 1, 1000.0, 0.0)
agg["season"] = (agg.index.month % 12) // 3
TSEAS = agg.groupby("season")["T"].max()
agg["SEAS"] = ampacity(agg["season"].map(TSEAS), 0.61, 90.0, 1000.0, Tmax=TMAX, elev=ELEV)
agg["AAR"] = ampacity(agg["T"], 0.61, 90.0, SUN(agg["day"]), Tmax=TMAX, elev=ELEV)
TF = agg["T"] * 9 / 5 + 32
agg["T881"] = (np.ceil((TF + 10) / 5) * 5 - 32) * 5 / 9
agg["AAR881"] = ampacity(agg["T881"], 0.61, 90.0, SUN(agg["day"]), Tmax=TMAX, elev=ELEV)
FIRST_TEST = {"jena": 202003, "m2": 202101}[SITE]   # models refitted monthly on all earlier data


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
    df["sd"] = np.sin(2 * np.pi * tgt.dayofyear / 365.25); df["cd"] = np.cos(2 * np.pi * tgt.dayofyear / 365.25)
    for var in ["V", "T", "G"]:
        df["y_" + var] = s[var].shift(-h)
        df[var + "1"] = s[var].shift(1)
    df["y_phi"] = s["phi"].shift(-h); df["phi0"] = s["phi"]
    df["y_day"] = s["day"].shift(-h)
    for k in ["AAR", "AAR881", "SEAS"]:
        df["y_" + k] = s[k].shift(-h)
    df["month"] = tgt.year * 100 + tgt.month
    df["t_target"] = tgt
    return df.dropna()


XA = ["A0", "A1", "A2", "A24", "Am24", "V0", "T0", "G0", "sh", "ch"]
XIN = {"V": ["V0", "V1", "T0", "G0", "sh", "ch"],
       "T": ["T0", "T1", "G0", "sh", "ch"],
       "G": ["G0", "G1", "T0", "sh", "ch"]}
XPA = ["y_T", "sh", "ch", "sd", "cd", "y_day"]       # probabilistic AAR: no wind or ampacity telemetry


def qfit(df_tr, df_te, ycol, xcols, q):
    Xtr = sm.add_constant(df_tr[xcols].values, has_constant="add")
    Xte = sm.add_constant(df_te[xcols].values, has_constant="add")
    m = QuantReg(df_tr[ycol].values, Xtr).fit(q=q, max_iter=5000)
    return Xte @ m.params, Xtr @ m.params


def causal(df, a):
    """Forecasts for every test month from models fitted on all earlier months."""
    out = []
    for mth in sorted(m for m in df.month.unique() if m >= FIRST_TEST):
        tr, te = df[df.month < mth], df[df.month == mth]
        r = pd.DataFrame(index=te.index)
        r["QR"], _ = qfit(tr, te, "y", XA, a)
        med_te, med_tr = qfit(tr, te, "y", XA, 0.5)
        r["Margin"] = med_te + np.quantile(tr["y"].values - med_tr, a)
        Vq, _ = qfit(tr, te, "y_V", XIN["V"], a)
        Tq, _ = qfit(tr, te, "y_T", XIN["T"], 1 - a)
        Gq, _ = qfit(tr, te, "y_G", XIN["G"], 1 - a)
        args = dict(Tmax=TMAX, elev=ELEV)
        r["Input"] = ampacity(Tq, np.maximum(Vq, 0), te["phi0"], np.clip(Gq, 0, None), **args)
        r["Input45"] = ampacity(Tq, np.maximum(Vq, 0), 45.0, np.clip(Gq, 0, None), **args)
        r["InputRealDir"] = ampacity(Tq, np.maximum(Vq, 0), te["y_phi"], np.clip(Gq, 0, None), **args)
        r["ProbAAR"], _ = qfit(tr, te, "y", XPA, a)
        out.append(r)
    return pd.concat(out).sort_index()


def quantile_tracking(q, y, issue, a, h, eta=ETA):
    """Additive online quantile tracking: R_t = q_t + theta_t, theta updated with each hit
    once it is observed (h hours after issue). Long-run coverage converges to a whenever
    theta stays bounded (Angelopoulos et al.); not adaptive conformal inference."""
    R = np.empty(len(q)); hits = np.zeros(len(q), int); th = 0.0; j = 0
    known = issue + pd.Timedelta(hours=h)                     # when each hit becomes known
    for i in range(len(q)):
        while j < i and known[j] <= issue[i]:
            th -= eta * (hits[j] - a); j += 1
        R[i] = q[i] + th; hits[i] = int(R[i] > y[i])
    return R


# ---------------------------------------------------------------- evaluate
def hourly(I, t):
    """Hit series on the full hourly grid of target times (gaps filled with 0)."""
    s = pd.Series(np.asarray(I, int), index=pd.DatetimeIndex(t))
    return s.reindex(pd.date_range(s.index.min(), s.index.max(), freq="h"), fill_value=0)


def evaluate(R, df, a, h):
    y = df["y"].values; I = (R > y).astype(int); T = len(I); x = int(I.sum())
    L = max(h - 1, 0)
    rate, p_uc_hac, ucb = vt.uc_hac(I, a, L)
    lr_uc, p_uc = vt.kupiec(I, a)
    ch = vt.christoffersen(I, a)
    tt = pd.DatetimeIndex(df["t_target"].values)
    Is = pd.Series(I, index=tt)
    lag = lambda k: Is.reindex(tt - pd.Timedelta(hours=k)).fillna(0).values
    lowV = (df["V0"].values < 1.0).astype(float)          # regime known at issue time
    X = np.column_stack([np.ones(T), lag(h), lag(2 * h), lowV, df["y_day"].values])
    _, p_dq = vt.dq_hac(I, a, X, L)
    rl = vt.runs(hourly(I, tt).values)
    hot = df["y_T"].values >= np.quantile(df["y_T"].values, 0.9)
    mag = (R[I == 1] - y[I == 1]) / y[I == 1] if x else np.array([np.nan])
    reg = {}
    for name, col in [("issue", "V0"), ("realised", "y_V")]:
        for lo, hi in [(0, 1), (1, 2), (2, 3), (3, 99)]:
            m = (df[col].values >= lo) & (df[col].values < hi); n = int(m.sum()); k = int(I[m].sum())
            reg[f"{name}_{lo}_{hi}"] = dict(n=n, k=k, rate=k / max(n, 1))
    return dict(T=T, x=x, rate=x / T, p_uc=p_uc if h == 1 else p_uc_hac, p_uc_kupiec=p_uc,
                p_uc_hac=p_uc_hac, ucb95=ucb, p_ind=ch["p_ind"], p_dq=p_dq,
                episodes=int(rl.size), max_run=int(rl.max()) if rl.size else 0,
                mag_med=float(np.nanmedian(mag)), mag_p95=float(np.nanpercentile(mag, 95)),
                mean_R=float(R.mean()), gain_vs_static=float(R.mean() / STATIC - 1),
                gain_hot=float(R[hot].mean() / STATIC - 1), regimes=reg)


RULES = ["QR", "QR+tracking", "Margin", "Input", "Input45", "ProbAAR"]


def main():
    results = []
    for h in HORIZONS:
        dfall = design(h)
        for a in ALPHAS:
            f = causal(dfall, a)
            df = dfall.loc[f.index]
            f["QR+tracking"] = quantile_tracking(f["QR"].values, df["y"].values, df.index, a, h)
            meths = {"Static": np.full(len(df), STATIC), "Seasonal": df["y_SEAS"].values,
                     "AAR": df["y_AAR"].values, "AAR881": df["y_AAR881"].values}
            meths.update({k: f[k].values for k in RULES + ["InputRealDir"]})
            for name, R in meths.items():
                r = evaluate(np.asarray(R, float), df, a, h); r.update(method=name, h=h, alpha=a)
                results.append(r)
                print(f"h={h:2d} a={a:.2f} {name:12s} x={r['x']:4d}/{r['T']} rate={r['rate']*100:5.2f}% "
                      f"pUC={r['p_uc']:.3f} ucb={r['ucb95']*100:5.2f}% pIND={r['p_ind']:.3f} pDQ={r['p_dq']:.3f} "
                      f"max={r['max_run']} mag95={r['mag_p95']*100:.1f}% gain={r['gain_vs_static']*100:+.1f}% "
                      f"hot={r['gain_hot']*100:+.1f}%", flush=True)
            if a == ALPHAS[0]:
                sav = df[["y", "V0", "y_V", "y_day", "y_T"]].copy()
                for k in RULES:
                    sav[k] = f[k].values
                for k in ["AAR", "AAR881", "SEAS"]:
                    sav[k] = df["y_" + k].values
                sav.to_csv(f"{OUT}/series_h{h}_a01.csv")

    ok = agg.dropna(subset=["A"])
    test = ok[(ok.index.year * 100 + ok.index.month) >= FIRST_TEST]
    def exc(col, frame):
        f2 = frame.dropna(subset=[col]); return float((f2[col] > f2["A"]).mean())
    sens = {}
    for V, phi in [(0.61, 90.0), (0.61, 45.0), (1.2, 90.0), (1.2, 45.0)]:
        R = ampacity(test["T881"], V, phi, SUN(test["day"]), Tmax=TMAX, elev=ELEV)
        sens[f"V{V}_phi{phi:.0f}"] = dict(rate=float((R > test["A"]).mean()),
                                          gain=float(R.mean() / STATIC - 1))
    summary = dict(static=STATIC, seasonal_T=TSEAS.to_dict(), n_hours=int(len(ok)),
                   n_test_hours=int(len(test)), first_test_month=FIRST_TEST,
                   mean_wind=float(ok["V"].mean()), frac_lowwind=float((ok["V"] < 1).mean()),
                   static_exceed=float((STATIC > ok["A"]).mean()),
                   aar881_wind_sensitivity=sens)
    json.dump(dict(summary=summary, results=results), open(f"{OUT}/results.json", "w"),
              indent=1, default=float)
    print(json.dumps(summary, indent=1, default=float))


if __name__ == "__main__":
    main()

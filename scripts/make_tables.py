"""Write the LaTeX body rows of Tables 6 (Jena) and 7 (M2) from results.json files.
Usage: make_tables.py <jena dir> <m2_ew dir> <m2_ns dir> <outdir>"""
import sys, json
jd, ed, nd, out = sys.argv[1:5]
NAMES = {"Static": "Static rating", "AAR": "AAR (measured $T_a$)", "QR": "QR",
         "QR+online": "QR + online correction", "Margin": "Median $-$ margin", "Input": "Input quantiles"}
RULES = ["QR", "QR+online", "Margin", "Input"]


def load(d):
    return json.load(open(f"{d}/results.json"))["results"]


def get(R, m, h, a=0.01):
    return next(r for r in R if r["method"] == m and r["h"] == h and r["alpha"] == a)


def p(x):
    return "$<$0.001" if x < 0.001 else (f"{x:.3f}" if x < 0.01 else f"{x:.2f}")


def gain(r):
    g = r["gain_vs_static"] * 100
    return "0" if abs(g) < 0.05 else f"${g:+.1f}$"


J, E, N = load(jd), load(ed), load(nd)
rows = []
for m in ["Static", "AAR"]:
    r = get(J, m, 1)
    rows.append(f"{NAMES[m]} & -- & {r['rate']*100:.2f} & -- & -- & -- & {r['max_run']} & "
                f"{r['mag_med']*100:.1f} / {r['mag_p95']*100:.1f} & {gain(r)} \\\\")
for h in [1, 24]:
    rows.append("\\midrule")
    for m in RULES:
        r = get(J, m, h)
        rows.append(f"{NAMES[m]} & {h} & {r['rate']*100:.2f} & {p(r['p_uc'])} & {p(r['p_ind'])} & "
                    f"{p(r['p_dq'])} & {r['max_run']} & {r['mag_med']*100:.1f} / {r['mag_p95']*100:.1f} & {gain(r)} \\\\")
open(f"{out}/table6_jena_rows.tex", "w").write("\n".join(rows) + "\n")

rows = []
for m in ["Static", "AAR"]:
    e, n = get(E, m, 1), get(N, m, 1)
    rows.append(f"{NAMES[m]} & -- & {e['rate']*100:.2f} & {n['rate']*100:.2f} & -- & -- & -- & -- & "
                f"{e['mag_p95']*100:.1f} & {n['mag_p95']*100:.1f} & {gain(e)} & {gain(n)} \\\\")
for h in [1, 24]:
    rows.append("\\midrule")
    for m in RULES:
        e, n = get(E, m, h), get(N, m, h)
        rows.append(f"{NAMES[m]} & {h} & {e['rate']*100:.2f} & {n['rate']*100:.2f} & {p(e['p_uc'])} & {p(n['p_uc'])} & "
                    f"{p(e['p_dq'])} & {p(n['p_dq'])} & {e['mag_p95']*100:.1f} & {n['mag_p95']*100:.1f} & {gain(e)} & {gain(n)} \\\\")
open(f"{out}/table7_m2_rows.tex", "w").write("\n".join(rows) + "\n")
print(f"wrote {out}/table6_jena_rows.tex and {out}/table7_m2_rows.tex")

"""Write the LaTeX body rows of Tables 6-8 from results.json files (v1.1).
Usage: make_tables.py <jena dir> <m2_ew dir> <m2_ns dir> <outdir>"""
import sys, json
dirs, out = sys.argv[1:4], sys.argv[4]
CASES = [json.load(open(f"{d}/results.json"))["results"] for d in dirs]
REF = [("Static", "Static (annual)"), ("Seasonal", "Seasonal (4 seasons)"),
       ("AAR", "AAR, no margin"), ("AAR881", "AAR, 10\\,$^\\circ$F margin")]
RULES = [("QR", "QR"), ("QR+tracking", "QR + tracking"), ("Margin", "Median $-$ margin"),
         ("Input", "Input q., persisted"), ("Input45", "Input q., 45$^\\circ$"),
         ("ProbAAR", "Prob.\\ AAR")]


def get(R, m, h, a):
    return next(r for r in R if r["method"] == m and r["h"] == h and r["alpha"] == a)


def p(x):
    return "$<$0.001" if x < 0.001 else (f"{x:.3f}" if x < 0.01 else f"{x:.2f}")


def g(x):
    v = x * 100
    return "0" if abs(v) < 0.05 else f"${v:+.1f}$"


rows = []
for m, lab in REF:
    rs = [get(R, m, 1, 0.01) for R in CASES]
    rows.append(lab + " & " + " & ".join(f"{r['rate']*100:.2f}" for r in rs) + " & "
                + " & ".join(g(r["gain_vs_static"]) for r in rs) + " & "
                + " & ".join(g(r["gain_hot"]) for r in rs) + " \\\\")
open(f"{out}/table6_reference_rows.tex", "w").write("\n".join(rows) + "\n")

for a, fn, cols in [(0.01, "table7_a01_rows.tex", ("rate", "p_uc", "p_dq", "gain")),
                    (0.05, "table8_a05_rows.tex", ("rate", "p_uc", "gain"))]:
    rows = []
    for h in [1, 24]:
        rows.append("\\midrule")
        for m, lab in RULES:
            cells = []
            for R in CASES:
                r = get(R, m, h, a)
                for c in cols:
                    cells.append(f"{r['rate']*100:.2f}" if c == "rate" else
                                 (g(r["gain_vs_static"]) if c == "gain" else p(r[c])))
            rows.append(f"{lab} & {h} & " + " & ".join(cells) + " \\\\")
    open(f"{out}/{fn}", "w").write("\n".join(rows[1:]) + "\n")
print("wrote table6_reference_rows.tex, table7_a01_rows.tex, table8_a05_rows.tex")

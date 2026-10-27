"""Write the LaTeX body rows of Tables 6-10 from results.json files (v1.2.1).
Usage: make_tables.py <jena dir> <m2_ew dir> <m2_ns dir> <outdir>"""
import sys, json
dirs, out = sys.argv[1:4], sys.argv[4]
CASES = [json.load(open(f"{d}/results.json"))["results"] for d in dirs]
REF = [("Static", "Static (annual)"), ("Seasonal", "Seasonal (4 seasons)"),
       ("AAR", "AAR (realised $T$), no margin"), ("AAR881", "AAR (realised $T$), $+10$\\,$^\\circ$F")]
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
                    (0.05, "table8_a05_rows.tex", ("rate", "p_uc", "p_dq", "gain"))]:
    rows = []
    for h in [1, 24]:
        rows.append("\\midrule")
        for m, lab in RULES:
            cells = []
            for R in CASES:
                r = get(R, m, h, a)
                for c in cols:
                    cells.append(f"{r['rate']*100:.2f} ({r['ucb95']*100:.2f})" if c == "rate" else
                                 (g(r["gain_vs_static"]) if c == "gain" else p(r[c])))
            rows.append(f"{lab} & {h} & " + " & ".join(cells) + " \\\\")
    open(f"{out}/{fn}", "w").write("\n".join(rows[1:]) + "\n")
# Table 9: exceedance frequency at h = 1, alpha = 0.01 against three reference definitions
SUB = [json.load(open(f"{d}/subhourly.json")) for d in dirs]
rows = []
for k, lab in [("QR", "QR"), ("QR+tracking", "QR + tracking"), ("Margin", "Median $-$ margin"),
               ("ProbAAR", "Prob.\\ AAR"), ("AAR881", "AAR (realised $T$), $+10$\\,$^\\circ$F")]:
    cells = [f"{S['rules'][k][c]*100:.2f}" for S in SUB for c in ("hourly", "mean10", "min10")]
    rows.append(lab + " & " + " & ".join(cells) + " \\\\")
open(f"{out}/table9_reference_rows.tex", "w").write("\n".join(rows) + "\n")
# Table 10: which rule each evaluation criterion selects (pre-specified criteria)
#   QS: lowest mean quantile score at alpha; coverage: frequency closest to alpha;
#   safety-first: largest mean gain among rules whose one-sided 95 % bound is <= alpha
#   (static and seasonal ratings included as candidates, so a fallback always exists).
SHORT = {"Static": "Static", "Seasonal": "Seasonal", "AAR881": "AAR $+10\\,^\\circ$F",
         "QR": "QR", "QR+tracking": "QR + tr.", "Margin": "Margin", "Input": "Input q.",
         "Input45": "Input q.\\ 45$^\\circ$", "ProbAAR": "Prob.\\ AAR"}
CAND = list(SHORT)
SITES = ["Jena", "M2 E--W", "M2 N--S"]
rows, summary = [], dict(cases=0, qs_unsafe=0, cov_unsafe=0, agree_qs=0)
for a in (0.01, 0.05):
    for h in (1, 24):
        for site, R in zip(SITES, CASES):
            rr = {m: get(R, m, h, a) for m in CAND}
            qs = min(CAND, key=lambda m: rr[m]["qscore"])
            cov = min(CAND, key=lambda m: abs(rr[m]["rate"] - a))
            safe = [m for m in CAND if rr[m]["ucb95"] <= a]
            sf = max(safe, key=lambda m: rr[m]["gain_vs_static"])
            summary["cases"] += 1
            summary["qs_unsafe"] += rr[qs]["ucb95"] > a
            summary["cov_unsafe"] += rr[cov]["ucb95"] > a
            summary["agree_qs"] += qs == sf
            cell = lambda m: f"{SHORT[m]} ({rr[m]['ucb95']*100:.2f})"
            rows.append(f"{a:.2f} & {h} & {site} & {cell(qs)} & {cell(cov)} & {SHORT[sf]} "
                        f"({rr[sf]['ucb95']*100:.2f}; {g(rr[sf]['gain_vs_static'])}) \\\\")
open(f"{out}/table10_selection_rows.tex", "w").write("\n".join(rows) + "\n")
json.dump(summary, open(f"{out}/table10_summary.json", "w"), indent=1, default=int)
print("wrote table6_reference_rows.tex, table7_a01_rows.tex, table8_a05_rows.tex, "
      "table9_reference_rows.tex, table10_selection_rows.tex")

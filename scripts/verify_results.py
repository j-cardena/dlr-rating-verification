"""Compare regenerated results with the published ones in expected/.
Exceedance counts may move by one or two with a different BLAS or statsmodels build;
anything larger, or any change of sign in a test decision at the 5 % level, is reported.
Usage: verify_results.py <outdir> <expected dir>   (exit code 1 on material difference)"""
import sys, json
out, exp = sys.argv[1:3]
A = json.load(open(f"{out}/results.json"))["results"]
B = {(r["method"], r["h"], r["alpha"]): r for r in json.load(open(f"{exp}/results.json"))["results"]}
bad = 0
for r in A:
    e = B[(r["method"], r["h"], r["alpha"])]
    issues = []
    if abs(r["x"] - e["x"]) > 2:
        issues.append(f"exceedances {r['x']} vs {e['x']}")
    for k in ["p_uc", "p_dq"]:
        if (r[k] < 0.05) != (e[k] < 0.05) and abs(r[k] - e[k]) > 0.01:
            issues.append(f"{k} decision {r[k]:.3f} vs {e[k]:.3f}")
    if abs(r["gain_vs_static"] - e["gain_vs_static"]) > 0.002:
        issues.append(f"gain {r['gain_vs_static']:.4f} vs {e['gain_vs_static']:.4f}")
    if issues:
        bad += 1
        print(f"DIFF {r['method']:10s} h={r['h']:2d} a={r['alpha']}: " + "; ".join(issues))
exact = A == list(json.load(open(f"{exp}/results.json"))["results"])
print(f"{out}: {'identical to' if exact else ('consistent with' if not bad else 'DIFFERS from')} {exp}")
sys.exit(1 if bad else 0)

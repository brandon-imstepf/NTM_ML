"""Supplementary Table S2 (R1.m1): error and R^2 of each sNTM against the numeric
NTM for each of the 24 simulations in the surrogates' training convention
(results/network_checks_l2f_consistent.csv, from network_checks_l2f("consistent")).
Writes results/table_sy_consistent.tex (a LaTeX table body for the supplement).

    python make_table_sy.py
"""
import csv
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
METHODS = ["mlp", "gp", "lr", "sr"]
SEEDS = ["ECl", "CA1", "RSPd"]


def sci(x):
    if not math.isfinite(x):
        return "n.f."
    m, e = f"{x:.1e}".split("e")
    return rf"${m}\mathrm{{e}}{{{int(e)}}}$"


def r2(x):
    return "n.f." if not math.isfinite(x) else f"${x:.3f}$"


rows = list(csv.DictReader(open(HERE / "results" / "network_checks_l2f_consistent.csv")))
truths = [t for t in csv.DictReader(open(HERE / "manifests" / "recovery_truths.csv")) if t["convention"] == "l2f"]
val = {(r["seed"], r["truth_id"], r["method"]): (float(r["nmse"]), float(r["r2"])) for r in rows}

lines = []
for seed in SEEDS:
    for t in truths:
        g = float(t["gamma1"])
        mant, ex = f"{g:.0e}".split("e")
        cells = [seed, rf"${mant}\times10^{{{int(ex)}}}$", t["lambda1"], f"{t['delta']}/{t['epsilon']}"]
        for m in METHODS:
            e, r = val[(seed, t["truth_id"], m)]
            if not math.isfinite(r):
                e = float("nan")
            cells += [sci(e), r2(r)]
        lines.append(" & ".join(cells) + r" \\")
    lines.append(r"\midrule")
lines[-1] = r"\bottomrule"
(HERE / "results" / "table_sy_consistent.tex").write_text("\n".join(lines) + "\n")
print(f"results/table_sy_consistent.tex: {sum(1 for l in lines if l.endswith(chr(92) * 2))} rows")

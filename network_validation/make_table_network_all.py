"""Supplementary tables of network error (R1.m1): NMSE of each sNTM against the
numeric NTM over all regions and time points, per simulation, from
results/network_error_all.csv (gp_manuscript/network_error_all.m).
Writes results/table_network_interior.tex (24 interior simulations) and
results/table_network_edge_ood.tex (8 edge and 8 out-of-domain settings): table rows
only, for the supplement. R^2 is not listed: in this metric R^2 = 1 - NMSE.

    python make_table_network_all.py
"""
import csv
import math
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
METHODS = ["lr", "sr", "interp", "gp", "mlp"]   # paper order (11/17)
SEEDS = ["ECl", "CA1", "RSPd"]
PRIMARY = ("CA1", "g0.005_l0.02_d20_e80")     # the Primary Parameter Set (Figs. 9-10 of the main text)


def fmt(v):
    if v is None:
        return "n.f."
    e = math.floor(math.log10(abs(v)))
    m = v / 10 ** e
    if round(m, 1) >= 10:
        m, e = m / 10, e + 1
    return f"${m:.1f}\\mathrm{{e}}{{{e}}}$"


def sci(x):
    e = math.floor(math.log10(x))
    m = x / 10 ** e
    return f"{m:g}\\times10^{{{e}}}" if abs(m - 1) > 1e-9 else f"10^{{{e}}}"


def par(row):
    g, l, d, e = (float(row[k]) for k in ("gamma", "lambda", "delta", "epsilon"))
    return f"${sci(g)}$ & {l:g} & {d:g}/{e:g}"


rows = list(csv.DictReader(open(HERE / "results" / "network_error_all.csv")))
val = defaultdict(dict)
meta = {}
for r in rows:
    key = (r["group"], r["seed"], r["setting"])
    val[key][r["method"]] = float(r["nmse"]) if r["finite"] == "1" else None
    meta[key] = r

order = lambda k: (float(meta[k]["gamma"]), float(meta[k]["lambda"]), -float(meta[k]["delta"]))
out = []
for s in SEEDS:
    keys = sorted([k for k in val if k[0] == "interior" and k[1] == s], key=order)
    for k in keys:
        cells = [s] + par(meta[k]).split(" & ") + [fmt(val[k].get(m)) for m in METHODS]
        if k[1] == PRIMARY[0] and k[2] == PRIMARY[1]:          # Primary Parameter Set in bold
            cells = [r"{\bfseries\boldmath " + c + "}" for c in cells]
        out.append(" & ".join(cells) + r" \\")
    out.append(r"\midrule")
(HERE / "results" / "table_network_interior.tex").write_text("\n".join(out[:-1]) + "\n", encoding="utf8")

out = []
for g, name in [("edge", "Edge"), ("ood", "OOD")]:
    keys = sorted([k for k in val if k[0] == g], key=order)
    for k in keys:
        out.append(f"{name} & {par(meta[k])} & " + " & ".join(fmt(val[k].get(m)) for m in METHODS) + r" \\")
    out.append(r"\midrule")
(HERE / "results" / "table_network_edge_ood.tex").write_text("\n".join(out[:-1]) + "\n", encoding="utf8")
print("wrote results/table_network_interior.tex, results/table_network_edge_ood.tex")

"""Supplementary Table S3: single-edge NMSE of each surrogate by position in the
training domain (results/domain_edge_error.csv and results/domain_edge_tails.csv,
from domain_edge_error.m). Writes results/table_domain_edge.tex (a LaTeX table
body for the supplement).

    python make_table_domain_edge.py
"""
import csv
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
METHODS = ["mlp", "gp", "lr", "sr"]
BINS = [("interior d>=0.15", r"$d \geq 0.15$ (interior)"), ("0.05-0.15", r"$0.05 \leq d < 0.15$"),
        ("0.01-0.05", r"$0.01 \leq d < 0.05$"), ("d<0.01", r"$d < 0.01$ (edge)")]
TAILS = [("N1", "low", r"$N_1$, lower 1\%"), ("N1", "high", r"$N_1$, upper 1\%")]


def sci(x):
    if not math.isfinite(x):
        return "n.f."
    m, e = f"{x:.1e}".split("e")
    return rf"${m}\mathrm{{e}}{{{int(e)}}}$"


depth = list(csv.DictReader(open(HERE / "results" / "domain_edge_error.csv")))
tails = list(csv.DictReader(open(HERE / "results" / "domain_edge_tails.csv")))
dv = {(r["bin"], r["method"], r["target"]): r for r in depth}
tv = {(r["input"], r["tail"], r["method"], r["target"]): r for r in tails}

lines = []
for key, label in BINS:
    n = dv[(key, "mlp", "flux")]["n"]
    cells = [label, f"{int(n):,}".replace(",", "{,}")]
    cells += [sci(float(dv[(key, m, t)]["nmse_mean"])) for t in ("flux", "w1") for m in METHODS]
    lines.append(" & ".join(cells) + r" \\")
lines.append(r"\midrule")
for inp, side, label in TAILS:
    n = tv[(inp, side, "mlp", "flux")]["n"]
    cells = [label, f"{int(n):,}".replace(",", "{,}")]
    cells += [sci(float(tv[(inp, side, m, t)]["nmse_mean"])) for t in ("flux", "w1") for m in METHODS]
    lines.append(" & ".join(cells) + r" \\")
lines.append(r"\bottomrule")
(HERE / "results" / "table_domain_edge.tex").write_text("\n".join(lines) + "\n")
print(f"results/table_domain_edge.tex: {len(BINS) + len(TAILS)} rows")

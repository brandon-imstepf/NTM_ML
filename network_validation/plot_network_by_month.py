"""Network error by month (R1.m1): NMSE (Eq. 7) and R^2 of the least-squares fit
over the 30 regions at months 1-5, for each sNTM in the 40 network settings of
results/network_error_by_month.csv (gp_manuscript/network_error_all.m).
Columns: interior (24 simulations), edge (8) and out-of-domain (8) settings (main text Table 3a). Lines: median over the simulations with finite concentrations; band:
interquartile range. Also writes results/table_network_by_month.tex (median
NMSE and R^2 per month, rows of a LaTeX table).

Needs numpy, pandas, matplotlib:
    python plot_network_by_month.py      # outputs_nm/network_by_month.png/.pdf

Font sizes 25% larger than before (11/16).
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
M = pd.read_csv(HERE / "results" / "network_error_by_month.csv")
# paper order (11/17): LR, SR, 1-NN, GP, MLP (legend left to right)
METHODS = [("lr", "LR"), ("sr", "SR"), ("interp", "1-NN"), ("gp", "GP"), ("mlp", "MLP")]
COLORS = {"mlp": "#1f77b4", "gp": "#ff7f0e", "interp": "#2ca02c", "lr": "#d62728", "sr": "#9467bd"}
GROUPS = [("interior", "Interior (24 simulations)"), ("edge", "Edge (8)"), ("ood", "Out of Domain (8)")]
# (11/17) R^2 = squared Pearson correlation, as in the network figure of the main text
# (column r2_fit of network_error_by_month.csv: R(1,2)^2 of corrcoef, i.e. the squared correlation)
R2_COL = "r2_fit"
R2_MIN = -0.02         # axis floor (R^2 = squared correlation is never negative)
plt.rcParams.update({"font.family": "serif", "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
                     "mathtext.fontset": "stix", "font.size": 15})

fig, axs = plt.subplots(2, 3, figsize=(14, 7.5), sharex=True)
for j, (g, glab) in enumerate(GROUPS):
    d = M[(M.group == g) & (M.finite == 1)]
    for i, (col, ylab, log) in enumerate([("nmse", "NMSE by month", True), (R2_COL, r"$R^2$ by month", False)]):
        ax = axs[i, j]
        for m, lab in METHODS:
            q = d[d.method == m].groupby("month")[col].quantile([0.25, 0.5, 0.75]).unstack()
            if not log:
                q = q.clip(lower=R2_MIN)
            ax.plot(q.index, q[0.5], "-o", color=COLORS[m], lw=1.8, ms=4, label=lab)
            ax.fill_between(q.index, q[0.25], q[0.75], color=COLORS[m], alpha=0.12, lw=0)
        if log:
            ax.set_yscale("log")
        else:
            ax.set_ylim(R2_MIN - 0.03, 1.02)
        ax.grid(alpha=0.25)
        if i == 0:
            ax.set_title(glab, fontsize=17.5)
        if j == 0:
            ax.set_ylabel(ylab, fontsize=16.25)
        if i == 1:
            ax.set_xlabel("Month", fontsize=16.25)
        ax.set_xticks(range(1, 6))
# legend in one row below the panels, left to right in the paper order
h, l = axs[0, 0].get_legend_handles_labels()
fig.legend(h, l, loc="lower center", ncol=len(METHODS), frameon=False, fontsize=15)
fig.tight_layout(rect=(0, 0.06, 1, 1))
out = HERE / "outputs_nm" / "network_by_month"
out.parent.mkdir(exist_ok=True)
fig.savefig(str(out) + ".png", dpi=250, facecolor="white")
fig.savefig(str(out) + ".pdf", facecolor="white")
print(str(out) + ".png")


def fmt(v):
    if not np.isfinite(v):
        return "--"
    if abs(v) >= 0.1:
        return f"{v:.2f}" if abs(v) < 10 else f"{v:.0f}"
    e = int(np.floor(np.log10(abs(v))))
    return f"${v / 10 ** e:.1f}\\mathrm{{e}}{{{e}}}$"


lines = []
for g, glab in GROUPS:
    d = M[(M.group == g) & (M.finite == 1)]
    lines.append(r"\multicolumn{11}{@{}l}{\textit{" + glab + r"}} \\")
    for m, lab in METHODS:
        e = d[d.method == m]
        nm = [fmt(e[e.month == k].nmse.median()) for k in range(1, 6)]
        r2 = [f"{e[e.month == k][R2_COL].median():.3f}" for k in range(1, 6)]
        lines.append(f"\\quad {lab} & " + " & ".join(nm) + " & " + " & ".join(r2) + r" \\")
    lines.append(r"\midrule")
(HERE / "results" / "table_network_by_month.tex").write_text("\n".join(lines[:-1]) + "\n", encoding="utf8")
print("results/table_network_by_month.tex")

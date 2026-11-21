"""Computational cost of griddatan 'linear' against 'nearest' (results/
griddatan_timing.csv, gp_manuscript/griddatan_timing.m; one CPU thread, the same
queries and training samples). Measured points, and for 'linear' a power-law fit
(least squares in log-log) extrapolated to 10^6 training samples (dashed).
(a) time per query; (b) set-up time (Delaunay triangulation; k-d tree);
(c) time of the surrogate evaluations of one network simulation, 4.2e4 queries
(276 edges x 151 time points), with the run time of the numeric NTM for
reference (14,254 s, Table 4 of the main text, one thread of a cluster node;
the timings here are on a different machine).

    python plot_griddatan_timing.py      # outputs_nm/griddatan_timing.png/.pdf

Font sizes 25% larger than before (11/16), and again 25% larger with one legend
under the panels (11/20, user).
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
T = pd.read_csv(HERE / "results" / "griddatan_timing.csv")
N_QUERIES = 276 * 151
T_NUMERIC = 14254.0
plt.rcParams.update({"font.family": "serif", "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
                     "mathtext.fontset": "stix", "font.size": 18.75})
COL = {"nearest": "#1f77b4", "linear": "#d62728"}
NAME = {"nearest": "1-NN ('nearest')", "linear": "Linear (Delaunay, 'linear')"}
n_ext = np.logspace(np.log10(500), 6, 60)

fig, axs = plt.subplots(1, 3, figsize=(16, 6.4), layout="constrained")
for j, (col, ylab) in enumerate([("query_s", "Time per query (s)"), ("setup_s", "Set-up time (s)"),
                                 ("query_s", f"Surrogate time per network run (s)\n({N_QUERIES:,} queries)")]):
    ax = axs[j]
    scale = N_QUERIES if j == 2 else 1
    for m in ["nearest", "linear"]:
        d = T[T.method == m].sort_values("n_train")
        ax.loglog(d.n_train, d[col] * scale, "o-", color=COL[m], lw=1.8, ms=5, label=f"{NAME[m]}, measured")
        if m == "linear":
            b, a = np.polyfit(np.log10(d.n_train), np.log10(d[col]), 1)
            ax.loglog(n_ext[n_ext >= d.n_train.max()], 10 ** (a + b * np.log10(n_ext[n_ext >= d.n_train.max()])) * scale,
                      "--", color=COL[m], lw=1.8, label=fr"{NAME[m]}, extrapolated ($\propto n^{{{b:.2f}}}$)")
    if j == 2:
        ax.axhline(T_NUMERIC, color="0.3", ls=":", lw=1.5, label="Numeric NTM run (Table 4)")
    ax.set_xlabel("Number of training samples, $n$", fontsize=20.3)
    ax.set_ylabel(ylab, fontsize=20.3)
    ax.grid(alpha=0.25, which="both")
    ax.set_title("abc"[j], loc="left", fontweight="bold", fontsize=23.4)
# one legend for the figure, under the panels (panel c has all four entries)
h, lab = axs[2].get_legend_handles_labels()
fig.legend(h, lab, loc="outside lower center", ncol=2, fontsize=15.6, frameon=False)
out = HERE / "outputs_nm" / "griddatan_timing"
fig.savefig(str(out) + ".png", dpi=250, facecolor="white")
fig.savefig(str(out) + ".pdf", facecolor="white")
print(str(out) + ".png")
for m in ["linear"]:
    d = T[T.method == m].sort_values("n_train")
    for col in ["query_s", "setup_s"]:
        b, a = np.polyfit(np.log10(d.n_train), np.log10(d[col]), 1)
        print(f"{m} {col}: exponent {b:.2f}; extrapolated at 1e4 {10 ** (a + b * 4):.3g} s, 1e5 {10 ** (a + b * 5):.3g} s, "
              f"1e6 {10 ** (a + b * 6):.3g} s")

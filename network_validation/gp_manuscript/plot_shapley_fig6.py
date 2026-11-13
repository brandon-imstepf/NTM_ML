"""Shapley figures in the style of the MLP Shapley figure of the paper
(images/Figure10.png): (a) shap.summary_plot beeswarms for phi_tau and
Delta_tau side by side with one colour bar, (b) a table of the mean |SHAP| of
each input.

Shapley values from results/shapley_<model>_<target>.csv (shapley_surrogates.m),
divided by the standard deviation of the target over the model's training
samples (results/shapley_target_std.csv), so that they are in units of the
Z-scored target like the MLP's.

Needs numpy, pandas, matplotlib and shap (tested with shap 0.52):
    python plot_shapley_fig6.py            # figures/shapley/shapley_fig6_<model>.png/.pdf
    python plot_shapley_fig6.py lr sr      # some models only
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import shap

HERE = Path(__file__).resolve().parent
RES = HERE.parent / "results"
OUT = HERE / "figures" / "shapley"
INPUTS = ["gamma1", "lambda1", "delta", "epsilon", "N1", "N2"]
LABELS = [r"$\gamma$", r"$\lambda$", r"$\delta$", r"$\epsilon$", r"$N_1$", r"$N_2$"]
TABLE_LABELS = ["γ", "λ", "δ", "ε", "N₁", "N₂"]
TARGETS = [("flux", r"$\phi_\tau$", "φτ"), ("emc", r"$\Delta_\tau$", "Δτ")]
MODELS = {"gp": ("Gaussian process", "GP"), "interp": ("Nearest-neighbor interpolation", "NN"),
          "lr": ("Linear regression", "LR"), "sr": ("Symbolic regression", "SR")}

plt.rcParams.update({"font.family": "serif", "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
                     "mathtext.fontset": "stix"})


def load(model, target):
    T = pd.read_csv(RES / f"shapley_{model}_{target}.csv")
    s = pd.read_csv(RES / "shapley_target_std.csv")
    sd = s[(s.model == model) & (s.target == target)].train_std.iloc[0]
    return T[["phi_" + n for n in INPUTS]].to_numpy() / sd, T[INPUTS].to_numpy()


def figure(model):
    name, short = MODELS[model]
    fig = plt.figure(figsize=(13, 8.6))
    fig.suptitle(f"{name} Shapley analysis feature importance", fontsize=22, fontweight="bold", y=0.985)
    fig.text(0.01, 0.905, "a", fontsize=20, fontweight="bold")
    fig.text(0.5, 0.875, "Distribution of feature importance SHAP values", ha="center", fontsize=17)
    means = []
    for k, (target, tlab, _) in enumerate(TARGETS):
        phi, X = load(model, target)
        means.append(np.abs(phi).mean(axis=0))
        ax = fig.add_axes([0.10 + 0.47 * k, 0.40, 0.37, 0.42])
        plt.sca(ax)
        shap.summary_plot(phi, X, feature_names=LABELS, show=False, plot_size=None, color_bar=(k == 1),
                          sort=True, max_display=6)
        ax.set_title(tlab, fontsize=18)
        ax.set_xlabel("")
        ax.tick_params(axis="y", labelsize=17)
        ax.tick_params(axis="x", labelsize=12)
        if k == 0:
            ax.set_ylabel("Model Feature", fontsize=17)
    fig.text(0.5, 0.345, "SHAP Value (Influence on Model Prediction)", ha="center", fontsize=17)
    # shap's colour bar: label and size as in the MLP figure
    cb = fig.axes[-1]
    if cb not in fig.axes[:2]:
        cb.set_ylabel("Feature Value", fontsize=17)
        cb.tick_params(labelsize=12)

    fig.text(0.01, 0.285, "b", fontsize=20, fontweight="bold")
    fig.text(0.5, 0.265, "Mean magnitude of SHAP values", ha="center", fontsize=17)
    tax = fig.add_axes([0.08, 0.02, 0.84, 0.22])
    tax.axis("off")
    cells = [[f"{v:.3f}" for v in m] for m in means]
    rows = [f"{t[2]} {short}\nModel" for t in TARGETS]
    tab = tax.table(cellText=cells, rowLabels=rows, colLabels=TABLE_LABELS, loc="center", cellLoc="center",
                    rowLoc="center", bbox=[0.12, 0, 0.88, 1])
    tab.auto_set_font_size(False)
    for (r, c), cell in tab.get_celld().items():
        cell.set_edgecolor("white")
        cell.set_linewidth(1.5)
        cell.set_fontsize(15)
        if r == 0:
            cell.set_facecolor("#1b5e80")
            cell.get_text().set_color("white")
            cell.get_text().set_fontweight("bold")
        else:
            cell.set_facecolor("#cdd1d6" if r == 1 else "#e3e6ea")
            if c == -1:
                cell.get_text().set_fontweight("bold")
    hdr = tab.add_cell(0, -1, width=tab[0, 0].get_width(), height=tab[0, 0].get_height(), text="Feature",
                       loc="center", facecolor="#1b5e80", edgecolor="white")
    hdr.get_text().set_color("white")
    hdr.get_text().set_fontweight("bold")
    hdr.set_fontsize(15)
    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"shapley_fig6_{model}.{ext}", dpi=300, facecolor="white")
    plt.close(fig)
    print(OUT / f"shapley_fig6_{model}.png")


if __name__ == "__main__":
    import sys
    for m in (sys.argv[1:] or MODELS):
        figure(m)

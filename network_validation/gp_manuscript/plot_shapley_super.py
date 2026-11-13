"""Shapley superfigure: the MLP, the Gaussian process and nearest-neighbor
interpolation in one figure, in the style of the MLP Shapley figure of the paper
(images/Figure10.png). (a) shap.summary_plot beeswarms, one row per model, with
phi_tau and Delta_tau side by side and a common x range per column; (b) a table
of the mean |SHAP| of each input for each model and target.

All three sets of Shapley values come from shapley_surrogates.m (exact
interventional Shapley values, 2,000 held-out samples, 50 background samples
from each model's training data), divided by the standard deviation of the
target over the model's training samples (results/shapley_target_std.csv), so
they are in units of the Z-scored target. The MLP values therefore differ from
those of Figure10.png, which were computed with the shap package.

    python plot_shapley_super.py           # figures/shapley/shapley_super.png/.pdf
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
MODELS = [("mlp", "MLP"), ("gp", "Gaussian process"), ("interp", "Nearest-neighbor\ninterpolation")]
SHORT = {"mlp": "MLP", "gp": "GP", "interp": "NN"}

plt.rcParams.update({"font.family": "serif", "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
                     "mathtext.fontset": "stix"})


def load(model, target):
    T = pd.read_csv(RES / f"shapley_{model}_{target}.csv")
    s = pd.read_csv(RES / "shapley_target_std.csv")
    sd = s[(s.model == model) & (s.target == target)].train_std.iloc[0]
    return T[["phi_" + n for n in INPUTS]].to_numpy() / sd, T[INPUTS].to_numpy()


def main():
    data = {(m, t): load(m, t) for m, _ in MODELS for t, _, _ in TARGETS}
    xlim = {}
    for t, _, _ in TARGETS:
        lo = min(data[m, t][0].min() for m, _ in MODELS)
        hi = max(data[m, t][0].max() for m, _ in MODELS)
        pad = 0.04 * (hi - lo)
        xlim[t] = (lo - pad, hi + pad)

    fig = plt.figure(figsize=(13, 17))
    fig.suptitle("Shapley analysis feature importance", fontsize=22, fontweight="bold", y=0.993)
    fig.text(0.01, 0.955, "a", fontsize=20, fontweight="bold")
    fig.text(0.555, 0.955, "Distribution of feature importance SHAP values", ha="center", fontsize=17)
    top, row_h, gap = 0.915, 0.185, 0.035
    for r, (m, mlab) in enumerate(MODELS):
        y0 = top - (r + 1) * row_h - r * gap
        for k, (t, tlab, _) in enumerate(TARGETS):
            phi, X = data[m, t]
            ax = fig.add_axes([0.20 + 0.40 * k, y0, 0.31, row_h])
            plt.sca(ax)
            shap.summary_plot(phi, X, feature_names=LABELS, show=False, plot_size=None,
                              color_bar=False, sort=True, max_display=6)
            ax.set_xlim(*xlim[t])
            ax.set_xlabel("")
            ax.tick_params(axis="y", labelsize=16)
            ax.tick_params(axis="x", labelsize=11)
            if r == 0:
                ax.set_title(tlab, fontsize=19)
        fig.text(0.04, y0 + row_h / 2, mlab, rotation=90, ha="center", va="center", fontsize=16,
                 fontweight="bold", multialignment="center")
    fig.text(0.555, top - 3 * row_h - 2 * gap - 0.035, "SHAP Value (Influence on Model Prediction)",
             ha="center", fontsize=17)
    # one colour bar for all panels, with shap's colour map
    cax = fig.add_axes([0.925, top - 3 * row_h - 2 * gap + 0.15, 0.012, 3 * row_h + 2 * gap - 0.3])
    sm = matplotlib.cm.ScalarMappable(cmap=shap.plots.colors.red_blue, norm=matplotlib.colors.Normalize(0, 1))
    cb = fig.colorbar(sm, cax=cax, ticks=[0, 1])
    cb.set_ticklabels(["Low", "High"])
    cb.outline.set_visible(False)
    cb.set_label("Feature Value", fontsize=17, labelpad=-5)
    cax.tick_params(labelsize=12, length=0)

    fig.text(0.01, 0.215, "b", fontsize=20, fontweight="bold")
    fig.text(0.555, 0.215, "Mean magnitude of SHAP values", ha="center", fontsize=17)
    tax = fig.add_axes([0.06, 0.01, 0.9, 0.19])
    tax.axis("off")
    cells, rows = [], []
    for t, _, tt in TARGETS:
        for m, _ in MODELS:
            cells.append([f"{v:.3f}" for v in np.abs(data[m, t][0]).mean(axis=0)])
            rows.append(f"{tt} {SHORT[m]}")
    tab = tax.table(cellText=cells, rowLabels=rows, colLabels=TABLE_LABELS, loc="center", cellLoc="center",
                    rowLoc="center", bbox=[0.13, 0, 0.87, 1])
    tab.auto_set_font_size(False)
    for (r, c), cell in tab.get_celld().items():
        cell.set_edgecolor("white")
        cell.set_linewidth(1.5)
        cell.set_fontsize(14)
        if r == 0:
            cell.set_facecolor("#1b5e80")
            cell.get_text().set_color("white")
            cell.get_text().set_fontweight("bold")
        else:
            cell.set_facecolor("#cdd1d6" if r <= 3 else "#e3e6ea")
            if c == -1:
                cell.get_text().set_fontweight("bold")
    hdr = tab.add_cell(0, -1, width=tab[0, 0].get_width(), height=tab[0, 0].get_height(), text="Feature",
                       loc="center", facecolor="#1b5e80", edgecolor="white")
    hdr.get_text().set_color("white")
    hdr.get_text().set_fontweight("bold")
    hdr.set_fontsize(14)
    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"shapley_super.{ext}", dpi=300, facecolor="white")
    plt.close(fig)
    print(OUT / "shapley_super.png")


if __name__ == "__main__":
    main()

"""Point-by-point comparison of Shapley values between surrogates: does each input
get the same attribution at the same input from the MLP and the GP (and from
nearest-neighbor interpolation)?

Data: results/shapley_<model>_<target>_commonbg.csv (shapley_surrogates.m with
background "common": the same 2,000 held-out inputs and the same 50 background
samples of data_bias_e5_1.csv for every model). All Shapley values are divided
by the same scale, the s.d. of the target over data_bias_e5_1.csv (the MLP row
of shapley_target_std.csv), so that the axes are directly comparable.

Figures (figures/shapley/):
  shapley_compare_mlp_<other>.png  one panel per input (columns) and target (rows):
      SHAP of <other> against SHAP of the MLP at each of the 2,000 inputs; dashed
      line y = x. R^2 about y = x: 1 - sum (y - x)^2 / sum (x - mean x)^2. The 5% of
      inputs with the largest total disagreement (sum over the six inputs of
      |SHAP_other - SHAP_MLP|) colours each point (log scale; blue = best
      agreement, red = worst), and the points are drawn in order of it.
  shapley_disagreement_mlp_<other>.png  total disagreement at each input against
      each of the six inputs (top: phi_tau, bottom: Delta_tau), and against the
      prediction errors of the two models on the numeric truth (held-out targets).
Writes results/shapley_compare.csv (R^2 and correlation per model pair, target
and input) and prints where the disagreement is concentrated.

    python plot_shapley_compare.py          # gp and interp against the MLP
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
RES = HERE.parent / "results"
OUT = HERE / "figures" / "shapley"
ROOT = HERE.parent.parent
INPUTS = ["gamma1", "lambda1", "delta", "epsilon", "N1", "N2"]
LABELS = [r"$\gamma$", r"$\lambda$", r"$\delta$", r"$\epsilon$", r"$N_1$", r"$N_2$"]
LOGX = [True, True, False, False, False, False]
TARGETS = [("flux", r"$\phi_\tau$", "FValue"), ("emc", r"$\Delta_\tau$", "W1")]
NAMES = {"mlp": "MLP", "gp": "GP", "interp": "1-NN"}
CMAP = "RdYlBu_r"     # blue = small disagreement (good), red = large (bad)

plt.rcParams.update({"font.family": "serif", "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
                     "mathtext.fontset": "stix", "font.size": 12})


def load(model, target, scale):
    T = pd.read_csv(RES / f"shapley_{model}_{target}_commonbg.csv")
    return T[["phi_" + n for n in INPUTS]].to_numpy() / scale, T


def r2_identity(x, y):
    return 1 - np.sum((y - x) ** 2) / np.sum((x - x.mean()) ** 2)


def truth(X):
    """Numeric phi_tau, Delta_tau at the query inputs (rows of the held-out set)."""
    H = pd.read_csv(ROOT / "data" / "data_neural_network_model" / "training_data_e6" / "data_bias_e5_11.csv")
    key = lambda A: [tuple(float(f"{v:.10g}") for v in r) for r in A]     # 10 significant digits
    idx = {k: i for i, k in enumerate(key(H[INPUTS[:4] + ["NRow", "NCol"]].to_numpy()))}
    rows = [idx[k] for k in key(X)]
    return H.iloc[rows][["FValue", "W1"]].to_numpy()


def compare(other):
    sd = pd.read_csv(RES / "shapley_target_std.csv")
    rows = []
    fig, axs = plt.subplots(2, 6, figsize=(20, 7.2), layout="constrained")
    fig2, axs2 = plt.subplots(2, 8, figsize=(25, 6.6), layout="constrained")
    for r, (t, tlab, tcol) in enumerate(TARGETS):
        scale = sd[(sd.model == "mlp") & (sd.target == t)].train_std.iloc[0]
        pm, Tm = load("mlp", t, scale)
        po, To = load(other, t, scale)
        X = Tm[INPUTS].to_numpy()
        assert np.array_equal(X, To[INPUTS].to_numpy())
        dis = np.abs(po - pm).sum(axis=1)
        top = dis >= np.quantile(dis, 0.95)
        # colour = total disagreement (log scale): blue = best agreement, red = worst;
        # worst points drawn last so they are not hidden
        order = np.argsort(dis)
        norm = LogNorm(vmin=np.quantile(dis, 0.01), vmax=dis.max())
        for i in range(6):
            ax = axs[r, i]
            x, y = pm[:, i], po[:, i]
            lim = [min(x.min(), y.min()), max(x.max(), y.max())]
            lim = [lim[0] - 0.05 * (lim[1] - lim[0]), lim[1] + 0.05 * (lim[1] - lim[0])]
            ax.plot(lim, lim, "--", color="0.45", lw=1)
            sc = ax.scatter(x[order], y[order], s=7, c=dis[order], cmap=CMAP, norm=norm, alpha=0.85, lw=0)
            ax.set_xlim(lim); ax.set_ylim(lim); ax.set_aspect("equal")
            r2, rho = r2_identity(x, y), np.corrcoef(x, y)[0, 1]
            ax.set_title(f"{LABELS[i]}   $R^2$ = {r2:.3f}", fontsize=13)
            if r == 1: ax.set_xlabel("MLP SHAP")
            if i == 0: ax.set_ylabel(f"{tlab}\n{NAMES[other]} SHAP", fontsize=14)
            rows.append(dict(other=other, target=t, input=INPUTS[i], r2_identity=r2, corr=rho,
                             mean_abs_mlp=np.abs(x).mean(), mean_abs_other=np.abs(y).mean(),
                             mean_abs_diff=np.abs(y - x).mean()))
        fig.colorbar(sc, ax=axs[r, :], shrink=0.9, pad=0.01,
                     label=r"Disagreement $\sum_i |\mathrm{SHAP}_i^{\mathrm{" + NAMES[other].split("-")[0] +
                     r"}} - \mathrm{SHAP}_i^{\mathrm{MLP}}|$")
        # where they disagree
        f_truth = truth(X)[:, r]
        err_o = (To.f_x.to_numpy() - f_truth) / scale
        err_m = (Tm.f_x.to_numpy() - f_truth) / scale
        for i in range(6):
            ax = axs2[r, i]
            sc2 = ax.scatter(X[order, i], dis[order], s=6, c=dis[order], cmap=CMAP, norm=norm, alpha=0.85, lw=0)
            if LOGX[i]: ax.set_xscale("log")
            ax.set_yscale("log")
            if r == 1: ax.set_xlabel(LABELS[i], fontsize=14)
            if i == 0: ax.set_ylabel(f"{tlab}\n" + r"$\sum_i |\Delta$SHAP$_i|$", fontsize=13)
        for j, (e, lab) in enumerate([(err_o, NAMES[other]), (err_m, "MLP")]):
            ax = axs2[r, 6 + j]
            ax.scatter(np.abs(e)[order], dis[order], s=6, c=dis[order], cmap=CMAP, norm=norm, alpha=0.85, lw=0)
            ax.set_xscale("log"); ax.set_yscale("log")
            c = np.corrcoef(np.log(np.abs(e) + 1e-12), np.log(dis))[0, 1]
            ax.set_title(f"corr (log) = {c:.2f}", fontsize=12)
            if r == 1: ax.set_xlabel(f"|{lab} prediction error|")
        fig2.colorbar(sc2, ax=axs2[r, :], shrink=0.9, pad=0.01, label="Disagreement")
        # summary of where the top 5% lie: median percentile of each input
        pct = {n: 100 * (pd.Series(X[:, i]).rank(pct=True)[top].median()) for i, n in enumerate(INPUTS)}
        print(f"{other} vs mlp, {t}: R2 per input " +
              ", ".join(f"{n} {row['r2_identity']:.3f}" for n, row in zip(INPUTS, rows[-6:])))
        print(f"   top-5% disagreement: median percentile of each input among all points: " +
              ", ".join(f"{n} {v:.0f}" for n, v in pct.items()))
        print(f"   |{other} error| median {np.median(np.abs(err_o)):.3g}, top-5% {np.median(np.abs(err_o[top])):.3g}; "
              f"|MLP error| median {np.median(np.abs(err_m)):.3g}, top-5% {np.median(np.abs(err_m[top])):.3g}")
    # (11/17) no suptitle on the comparison figure: it is a Supplementary figure and its caption says this
    fig2.suptitle(f"Where {NAMES[other]} and MLP Shapley values disagree", fontsize=16, fontweight="bold")
    OUT.mkdir(parents=True, exist_ok=True)
    for f, name in [(fig, f"shapley_compare_mlp_{other}"), (fig2, f"shapley_disagreement_mlp_{other}")]:
        f.savefig(OUT / f"{name}.png", dpi=200, facecolor="white")
        f.savefig(OUT / f"{name}.pdf", facecolor="white")
        plt.close(f)
    return rows


if __name__ == "__main__":
    rows = []
    for o in ["gp", "interp"]:
        rows += compare(o)
    pd.DataFrame(rows).to_csv(RES / "shapley_compare.csv", index=False)
    print(RES / "shapley_compare.csv")

"""Supplementary table of all 22 symbolic regression expressions (11 PySR runs for
phi_tau, 11 for Delta_tau): for each run, the expression PySR's "best" rule
selected (best.txt), its complexity (model_equations.csv of that run), PySR's own
training and test NMSE (best.txt), and the NMSE (Eq. 7) on the 10^5-sample
held-out set data_bias_e5_11.csv (never used for fitting or selection), computed
with that run's input/target scaling (scaling_params.mat). The run used in the
paper is marked; for Delta_tau that run's selected expression is the one used, for
phi_tau the paper uses the complexity-41 expression of the 2025-10-21 run (Fig. 4).

Writes results/sr_runs.csv and results/table_sr_runs.tex.
Run with the UCM_AD_ML venv (numpy, scipy, pandas):
  AD_ML_Venv/Scripts/python.exe make_sr_runs_table.py
"""
import os
import re
import numpy as np
import pandas as pd
from scipy.io import loadmat

HERE = os.path.dirname(os.path.abspath(__file__))
SR_ROOT = r"C:\Users\USER\Documents\Code\UCM_AD_ML\datasets\sim_csv_v11\train"
HELDOUT = os.path.join(HERE, "..", "data", "data_neural_network_model", "training_data_e6", "data_bias_e5_11.csv")
TEST = os.path.join(HERE, "..", "data", "data_baseline_models", "e3_data", "bias_val_e3.csv")   # Table 1 test set
RUNS = {"phi": "flux_only_bias_train_e3", "Delta": "w1_only_bias_train_e3"}
USED = {"phi": ("2025-10-21_09-36_259200_Info", 41), "Delta": ("2025-10-21_09-14_259200_Info", None)}
SYM = ["\\gamma", "\\lambda", "\\delta", "\\epsilon", "N_1", "N_2"]   # x0..x5


def nmse(p, y):
    return np.mean((p - y) ** 2) / np.var(y)


def evaluate(expr, X):
    env = {f"x{k}": X[:, k] for k in range(6)}
    env.update(exp=np.exp, log=np.log, sqrt=np.sqrt, abs=np.abs)
    with np.errstate(all="ignore"):
        v = eval(expr, {"__builtins__": {}}, env)
    return np.broadcast_to(v, X.shape[:1]).astype(float)


def to_latex(expr):
    """sympy-style string -> LaTeX via sympy, inputs renamed to their symbols."""
    import sympy as sp
    xs = sp.symbols("x0:6")
    e = sp.sympify(expr, locals={f"x{k}": xs[k] for k in range(6)})
    s = sp.latex(e)
    s = re.sub(r"(\d) \{", r"\1 \\cdot {", s)       # number times a power: 0.509 \cdot {0.888}^{...}
    s = re.sub(r"x_\{(\d)\}", lambda m: SYM[int(m.group(1))], s)
    s = re.sub(r"x_(\d)", lambda m: SYM[int(m.group(1))], s)
    s = re.sub(r"(\d+\.\d{3})\d+", r"\1", s)       # 3 decimals at most
    return s


H = pd.read_csv(HELDOUT).to_numpy(float)
if False:
    H = H[1:]
X, Y = H[:, :6], H[:, 6:8]
Tt = pd.read_csv(TEST).to_numpy(float)
Xt, Yt = Tt[:, :6], Tt[:, 6:8]
rows = []
for ti, (tgt, sub) in enumerate(RUNS.items()):
    for d in sorted(os.listdir(os.path.join(SR_ROOT, sub))):
        if not d.endswith("_Info"):
            continue
        p = os.path.join(SR_ROOT, sub, d)
        best = open(os.path.join(p, "best.txt")).read()
        expr = re.search(r"Selected Best Function:\s*\d+,\s*(.*)", best).group(1).strip()
        tr = float(re.search(r"Train NMSE:\s*([\d.eE+-]+)", best).group(1))
        te = float(re.search(r"Test NMSE:\s*([\d.eE+-]+)", best).group(1))
        E = pd.read_csv(os.path.join(p, "model_equations.csv"))
        cx = E.loc[E["sympy_format"].str.strip() == expr, "complexity"]
        used_run, used_cx = USED[tgt]
        if d == used_run and used_cx is not None:      # phi_tau: the complexity-41 expression
            r = E.loc[E["complexity"] == used_cx].iloc[0]
            expr, cx = r["sympy_format"].strip(), pd.Series([used_cx])
        S = loadmat(os.path.join(p, "scaling_params.mat"))
        fm, fs = S["feature_means"].ravel(), S["feature_stds"].ravel()
        tm, tsd = float(np.squeeze(S["target_mean"])), float(np.squeeze(S["target_std"]))
        pred = evaluate(expr, (X - fm) / fs) * tsd + tm
        pt = evaluate(expr, (Xt - fm) / fs) * tsd + tm
        rows.append(dict(target=tgt, run=d[:16], complexity=int(cx.iloc[0]) if len(cx) else np.nan,
                         expression=expr, nmse_train_pysr=tr, nmse_test_pysr=te,
                         nmse_test=nmse(pt, Yt[:, ti]), nmse_heldout=nmse(pred, Y[:, ti]), used=(d == used_run)))
T = pd.DataFrame(rows)
T.to_csv(os.path.join(HERE, "results", "sr_runs.csv"), index=False)
print(T[["target", "run", "complexity", "nmse_test_pysr", "nmse_test", "nmse_heldout", "used"]].to_string())

lines = []
for tgt, name in [("phi", "\\phi_\\tau"), ("Delta", "\\Delta_\\tau")]:
    lines.append(f"\\multicolumn{{5}}{{l}}{{${name}$}} \\\\")
    for k, r in enumerate(T[T.target == tgt].itertuples(), 1):
        mark = "*" if r.used else ""
        lines.append(f"{k}{mark} & {r.complexity} & ${to_latex(r.expression)}$ & {r.nmse_test:.2f} & {r.nmse_heldout:.2f} \\\\")
    lines.append("\\midrule")
open(os.path.join(HERE, "results", "table_sr_runs.tex"), "w").write("\n".join(lines[:-1]) + "\n")
print("wrote results/sr_runs.csv and results/table_sr_runs.tex")

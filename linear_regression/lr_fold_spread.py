"""Spread of the 10-fold cross-validated NMSE of the selected linear regressions
(R2.M3: variability across independent training runs).

Uses the ten features of each selected model as listed in Table 2 of the paper
(TABLE2), on data/data_baseline_models/e3_data/bias_train_e3.csv, with the folds of
feature_selection_and_importance.ipynb (KFold(10, shuffle=True, random_state=42);
NMSE = MSE / Var as in Eq. 7), and reports the validation NMSE of each of the 10
folds. Rerunning the notebook's forward selection here (forward_stepwise) gives the
same features for Delta_tau, and for phi_tau eps*N2 in place of eps (7th of 10).
As a check against Table 1, the test NMSE on bias_val_e3.csv (1000 bootstrap
resamples of 1000, random_state 42 + i, as classical_regressions.ipynb) is printed
for the model of the best fold and for the OLS fit on all training data (the
notebook's last cell).

Writes linear_regression/results/lr_fold_spread.csv.

    python lr_fold_spread.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import KFold
from sklearn.utils import resample

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data" / "data_baseline_models" / "e3_data"
ELBOW = 10
K = 10
FEATS_LIN = ["x1", "x2", "x3", "x4", "x5", "x6"]
# x1..x6 = gamma, lambda, delta, epsilon, N1, N2; order of Table 2
TABLE2 = {"FValue": ["x6", "x1", "x3x5", "x1x1", "x1x5", "x1x3", "x4", "x2x5", "x1x6", "x2x6"],
          "W1": ["x1", "x1x1", "x2x3", "x4x5", "x4x6", "x1x2", "x2", "x2x2", "x2x5", "x2x4"]}


def build_features(data):
    df = pd.DataFrame()
    for name, col in zip(FEATS_LIN, ["gamma1", "lambda1", "delta", "epsilon", "NRow", "NCol"]):
        df[name] = data[col]
    for i, fi in enumerate(FEATS_LIN):
        for j in range(i + 1, len(FEATS_LIN)):
            df[fi + FEATS_LIN[j]] = df[fi] * df[FEATS_LIN[j]]
    for f in FEATS_LIN:
        df[f + f] = df[f] ** 2
    return df


def nmse(y, p):
    return np.mean((y - p) ** 2) / np.var(y)


def fold_scores(X, y, kf):
    s = np.zeros(K)
    models = []
    for k, (tr, te) in enumerate(kf.split(X)):
        m = LinearRegression().fit(X[tr], y[tr])
        s[k] = nmse(y[te], m.predict(X[te]))
        models.append(m)
    return s, models


def forward_stepwise(X, y, kf):
    selected, remaining = [], list(X.columns)
    while remaining:
        best = min(remaining, key=lambda f: fold_scores(np.array(X[selected + [f]]), y, kf)[0].mean())
        selected.append(best)
        remaining.remove(best)
    return selected


train = pd.read_csv(DATA / "bias_train_e3.csv")
val = pd.read_csv(DATA / "bias_val_e3.csv")
Xtr_all, Xva_all = build_features(train), build_features(val)
rows = []
for target, name in [("FValue", "phi_tau"), ("W1", "Delta_tau")]:
    kf = KFold(n_splits=K, shuffle=True, random_state=42)
    y, yv = np.array(train[target]), np.array(val[target])
    feats = TABLE2[target]
    s, models = fold_scores(np.array(Xtr_all[feats]), y, kf)
    best = int(np.argmin(s))
    Xv = np.array(Xva_all[feats])
    def boot(m):
        b = [nmse(yb, m.predict(Xb)) for Xb, yb in
             (resample(Xv, yv, n_samples=1000, random_state=42 + i) for i in range(1000))]
        return np.mean(b), np.std(b)
    full = LinearRegression().fit(np.array(Xtr_all[feats]), y)
    print(f"{name}: features {feats}")
    print(f"  fold validation NMSE: {np.array2string(s, precision=3)}")
    print(f"  mean {s.mean():.3f} +/- {s.std(ddof=1):.3f} (sd), median {np.median(s):.3f}, "
          f"range {s.min():.3f}-{s.max():.3f}; best fold {best + 1}")
    print("  test NMSE (1000 x 1000 bootstrap): best-fold model %.3f +/- %.3f; full-data fit %.3f +/- %.3f"
          % (*boot(models[best]), *boot(full)))
    for k in range(K):
        rows.append({"target": name, "fold": k + 1, "val_nmse": s[k], "selected": k == best,
                     "features": " ".join(feats)})
out = HERE / "results"
out.mkdir(exist_ok=True)
pd.DataFrame(rows).to_csv(out / "lr_fold_spread.csv", index=False)
print(f"wrote {out / 'lr_fold_spread.csv'}")

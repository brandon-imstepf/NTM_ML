"""Refit the feature-selected linear regression surrogates and export them for MATLAB.

Reproduces linear_regression/feature_selection_and_importance.ipynb and the
model-selection rule in the paper's Methods:
  1. forward stepwise selection over the 27 'All' features (10-fold CV, seed 42)
  2. keep the first ELBOW features
  3. refit those features in each of the same 10 folds and keep the fold
     with the lowest validation NMSE
Then reports bootstrapped train/test NMSE as in classical_regressions.ipynb
(1000 resamples of 1000 rows, seeds 42+i) for comparison with Table 1.

Usage:
    python fit_lr_surrogate.py [out.mat]
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import KFold
from sklearn.utils import resample

ROOT = Path(__file__).resolve().parent.parent
TRAIN = ROOT / "data" / "data_baseline_models" / "e3_data" / "bias_train_e3.csv"
TEST = ROOT / "data" / "data_baseline_models" / "e3_data" / "bias_val_e3.csv"
ELBOW = 10
K_FOLDS = 10
BOOTS = 1000

FEATS_LIN = ["x1", "x2", "x3", "x4", "x5", "x6"]
SYMBOL = {"x1": "gamma", "x2": "lambda", "x3": "delta", "x4": "epsilon", "x5": "N1", "x6": "N2"}
# Feature order in the paper's tab:linreg-ttest
PAPER_ORDER = {
    "FValue": ["x6", "x1", "x3x5", "x1x1", "x1x5", "x1x3", "x4", "x2x5", "x1x6", "x2x6"],
    "W1": ["x1", "x1x1", "x2x3", "x4x5", "x4x6", "x1x2", "x2", "x2x2", "x2x5", "x2x4"],
}


def build_features(data):
    df = pd.DataFrame({f: data[c] for f, c in zip(
        FEATS_LIN, ["gamma1", "lambda1", "delta", "epsilon", "NRow", "NCol"])})
    for i, fi in enumerate(FEATS_LIN):
        for fj in FEATS_LIN[i + 1:]:
            df[fi + fj] = df[fi] * df[fj]
    for f in FEATS_LIN:
        df[f + f] = df[f] ** 2
    return df


def nmse(y, p):
    # the notebooks' definition (MSE / variance), as used for Table 1
    return mean_squared_error(y, p) / np.var(y)


def forward_stepwise(X, y):
    kf = KFold(n_splits=K_FOLDS, shuffle=True, random_state=42)
    selected, remaining = [], list(X.columns)
    while remaining:
        scores = {}
        for f in remaining:
            Xs = X[selected + [f]].values
            scores[f] = np.mean([
                nmse(y[te], LinearRegression().fit(Xs[tr], y[tr]).predict(Xs[te]))
                for tr, te in kf.split(Xs)])
        best = min(remaining, key=lambda f: scores[f])  # first minimum, as in the notebook
        selected.append(best)
        remaining.remove(best)
    return selected


def pretty(f):
    a, b = f[:2], f[2:]
    if not b:
        return SYMBOL[a]
    return f"{SYMBOL[a]}^2" if a == b else f"{SYMBOL[a]}*{SYMBOL[b]}"


def fit_target(target, train, test):
    Xtr, Xte = build_features(train), build_features(test)
    ytr, yte = train[target].values, test[target].values

    order = forward_stepwise(Xtr, ytr)
    feats = order[:ELBOW]
    match = feats == PAPER_ORDER[target]
    print(f"\n{target}: selected {[pretty(f) for f in feats]}")
    print(f"  matches paper feature order: {match}")

    kf = KFold(n_splits=K_FOLDS, shuffle=True, random_state=42)
    best, best_score = None, np.inf
    for tr, te in kf.split(Xtr):
        m = LinearRegression().fit(Xtr[feats].values[tr], ytr[tr])
        s = nmse(ytr[te], m.predict(Xtr[feats].values[te]))
        if s < best_score:
            best, best_score = m, s

    def boot(X, y):
        return np.array([nmse(yb, best.predict(Xb)) for Xb, yb in (
            resample(X[feats].values, y, n_samples=1000, random_state=42 + i) for i in range(BOOTS))])
    b_tr, b_te = boot(Xtr, ytr), boot(Xte, yte)
    print(f"  bootstrapped NMSE  train {b_tr.mean():.3f} +/- {b_tr.std():.3f}   test {b_te.mean():.3f} +/- {b_te.std():.3f}")

    return {
        "features": np.array(feats, dtype=object),
        "feature_symbols": np.array([pretty(f) for f in feats], dtype=object),
        "coef": best.coef_.astype(np.float64),
        "intercept": float(best.intercept_),
        "matches_paper_order": match,
        "boot_nmse_train": b_tr, "boot_nmse_test": b_te,
    }


def main():
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent / "models" / "lr_weights.mat"
    out.parent.mkdir(parents=True, exist_ok=True)
    train, test = pd.read_csv(TRAIN), pd.read_csv(TEST)
    scipy.io.savemat(out, {"flux": fit_target("FValue", train, test), "w1": fit_target("W1", train, test)})
    print(f"\nsaved {out}")


if __name__ == "__main__":
    main()

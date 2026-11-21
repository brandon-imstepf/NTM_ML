"""Forward feature selection curve of the linear regression surrogates (Fig. 3).

Recomputes the curve of linear_regression/feature_selection_and_importance.ipynb
(whose saved .npy outputs are lost), with the same features, folds and rule as
fit_lr_surrogate.py: at each step the feature with the lowest mean 10-fold
cross-validation NMSE is added. For the feature added at each step, records the
mean and standard deviation over the folds of the training NMSE (on the 9 fitting
folds) and of the cross-validation NMSE (on the held-out fold).

Usage:
    python lr_selection_curve.py      -> results/lr_selection_curve.csv
"""
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import KFold

from fit_lr_surrogate import ELBOW, K_FOLDS, PAPER_ORDER, TEST, TRAIN, build_features, nmse, pretty

OUT = Path(__file__).resolve().parent / "results" / "lr_selection_curve.csv"


def curve(X, y):
    # z-score the features: with an intercept the least-squares fit (and so every
    # NMSE) is unchanged in exact arithmetic, but on the raw inputs (N^2 ~ 1e-12)
    # the solve loses precision once ~22 features are in and the training NMSE
    # rises, which a nested least-squares fit cannot do (11/19)
    X = (X - X.mean()) / X.std()
    kf = KFold(n_splits=K_FOLDS, shuffle=True, random_state=42)
    folds = list(kf.split(X))
    selected, remaining, rows = [], list(X.columns), []
    while remaining:
        best = None
        for f in remaining:
            Xs = X[selected + [f]].values
            tr_s, cv_s = [], []
            for tr, te in folds:
                m = LinearRegression().fit(Xs[tr], y[tr])
                tr_s.append(nmse(y[tr], m.predict(Xs[tr])))
                cv_s.append(nmse(y[te], m.predict(Xs[te])))
            if best is None or np.mean(cv_s) < best[1]:     # first minimum, as in the notebook
                best = (f, np.mean(cv_s), np.std(cv_s), np.mean(tr_s), np.std(tr_s))
        selected.append(best[0])
        remaining.remove(best[0])
        rows.append({"n_features": len(selected), "feature": best[0], "symbol": pretty(best[0]),
                     "train_nmse_mean": best[3], "train_nmse_sd": best[4],
                     "cv_nmse_mean": best[1], "cv_nmse_sd": best[2]})
    return pd.DataFrame(rows)


def main():
    train = pd.read_csv(TRAIN)
    X = build_features(train)
    out = []
    for target, name in [("FValue", "flux"), ("W1", "emc")]:
        T = curve(X, train[target].values)
        T.insert(0, "target", name)
        ok = list(T.feature[:ELBOW]) == PAPER_ORDER[target]
        print(f"{target}: first {ELBOW} features match the paper's order: {ok}")
        print(T.head(12).to_string(index=False))
        out.append(T)
    pd.concat(out).to_csv(OUT, index=False)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()

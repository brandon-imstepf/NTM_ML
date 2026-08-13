from sklearn.metrics import mean_squared_error, r2_score
from scipy.stats import pearsonr
import numpy as np


def model_fit(y_true, y_pred, verbose=True):
    
    pr = pearsonr(y_true, y_pred)
    r2 = r2_score(y_true, y_pred)
    mse = mean_squared_error(y_true, y_pred)
    var = np.var(y_true)

    if verbose:
        print(f"R^2: {r2}")
        print(f"Pearson R: {pr[0]}")
        print(f"MSE: {mse}")
        print(f"RMSE: {mse**0.5}")
        print(f"NMSE: {mse/var}")
        print(f"NRMSE: {(mse**0.5)/var}")

    return r2, pr[0], mse, mse**0.5, mse/var, (mse**0.5)/var

"""Export the final e6 MLPs (TorchScript .pt) to a .mat file for mlp_forward.m.

Usage:
    python export_mlp_weights.py <flux_model.pt> <w1_model.pt> [out.mat]

The scalers are refit exactly as nn-hyp-tuning.py does: StandardScaler on
data_bias_e5_1.csv. The script checks that the NumPy forward pass reproduces
the TorchScript output before saving.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io
import torch

FEATURES = ["gamma1", "lambda1", "delta", "epsilon", "NRow", "NCol"]
ROOT = Path(__file__).resolve().parent.parent
SCALE_DATA = ROOT / "data" / "data_neural_network_model" / "training_data_e6" / "data_bias_e5_1.csv"


def layer_order(name):
    if "layer_first" in name:
        return -1
    if "layer_last" in name:
        return 10**6
    return int(name.split("layer_hidden_")[1].split(".")[0])


def export(pt_path, target, df):
    model = torch.jit.load(pt_path, map_location="cpu").eval()
    sd = {k: v.detach().numpy().astype(np.float64) for k, v in model.state_dict().items()}
    names = sorted({k.rsplit(".", 1)[0] for k in sd}, key=layer_order)
    W = [sd[n + ".weight"] for n in names]
    b = [sd[n + ".bias"] for n in names]

    X = df[FEATURES].values
    y = df[target].values
    net = {
        "x_mean": X.mean(0), "x_std": X.std(0),  # population std, as StandardScaler
        "y_mean": y.mean(), "y_std": y.std(),
    }

    # Self-check against TorchScript on 1000 rows
    Xs = (X[:1000] - net["x_mean"]) / net["x_std"]
    with torch.no_grad():
        ref = model(torch.tensor(Xs, dtype=torch.float32)).numpy().ravel()
    h = Xs.T
    for k in range(len(W)):
        h = W[k] @ h + b[k][:, None]
        if k < len(W) - 1:
            h = np.maximum(h, 0)
    err = np.max(np.abs(h.ravel() - ref))
    print(f"{target}: {len(W)} linear layers, widths {[w.shape[0] for w in W]}, max |numpy - torch| = {err:.2e}")
    assert err < 1e-4, "NumPy forward pass does not match TorchScript"

    net["W"] = np.empty(len(W), dtype=object)
    net["b"] = np.empty(len(b), dtype=object)
    for k in range(len(W)):
        net["W"][k] = W[k]
        net["b"][k] = b[k]
    return net


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    out = Path(sys.argv[3]) if len(sys.argv) > 3 else Path(__file__).resolve().parent / "models" / "mlp_e6_weights.mat"
    out.parent.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(SCALE_DATA)
    scipy.io.savemat(out, {
        "flux": export(sys.argv[1], "FValue", df),
        "w1": export(sys.argv[2], "W1", df),
    })
    print(f"saved {out}")


if __name__ == "__main__":
    main()

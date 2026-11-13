"""Manifest for the numeric network runs at the 8 edge and 8 out-of-domain (ood)
settings (ECl seed) in the surrogates' training convention, as for the interior
runs of make_consistent_manifests.py (outputs_consistent/): flux at lambda2 = 0.01,
W1 at lambda2 = lambda1 (SNTM_LAMBDA2_W1=lambda1), each edge's flux solved
separately (SNTM_FLUX_PER_EDGE=1).

The settings are those of manifests/numeric.csv (groups edge and ood). The
surrogate runs need no rerun: they do not depend on lambda2 (identical N in
outputs/ and outputs_consistent/ for the interior settings, checked 11/12), so
the surrogate files in outputs/ are used with these numeric runs.

    python make_edge_ood_manifest.py
    # then, in MATLAB with SNTM_OUT_DIR=outputs_consistent, SNTM_LAMBDA2_W1=lambda1,
    # SNTM_FLUX_PER_EDGE=1:
    #   run_local('manifests/consistent_edge_ood_numeric.csv', 8)
"""
import csv
from pathlib import Path

HERE = Path(__file__).resolve().parent
COLUMNS = ["case_id", "method", "setting_id", "group", "seed", "subset",
           "gamma1", "lambda1", "delta", "epsilon", "init_rescale", "lambda2"]

with open(HERE / "manifests" / "numeric.csv", newline="") as f:
    src = [r for r in csv.DictReader(f) if r["group"] in ("edge", "ood")]

rows = []
for r in src:
    sid = r["setting_id"].rsplit("__", 1)[0]          # e.g. edge__g5e-05_l0.01_d10_e100
    rows.append({
        "case_id": f"{sid}__{r['seed']}__numeric_l2f",
        "method": "numeric",
        "setting_id": f"{sid}__{r['seed']}",
        "group": f"consistent_{r['group']}", "seed": r["seed"], "subset": r["subset"],
        "gamma1": r["gamma1"], "lambda1": r["lambda1"],
        "delta": r["delta"], "epsilon": r["epsilon"],
        "init_rescale": r["init_rescale"], "lambda2": 0.01})

out = HERE / "manifests" / "consistent_edge_ood_numeric.csv"
with open(out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=COLUMNS)
    w.writeheader()
    w.writerows(rows)
print(f"{out.relative_to(HERE)}: {len(rows)} rows")

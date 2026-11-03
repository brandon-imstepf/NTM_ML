"""Manifests for the network runs in the surrogates' training convention
(solver_check/SOLVER_CHECK.md, sections 8-9; decided 11/03/26).

consistent_l2f_numeric.csv, consistent_l2f_surrogate.csv
    the 8 l2f settings of recovery_truths.csv seeded in ECl, CA1 and RSPd:
    24 numeric and 96 surrogate cases, with the same case ids as the earlier
    runs. Run with SNTM_LAMBDA2_W1=lambda1 (W1 at lambda2 = lambda1, as in the
    training data; flux keeps lambda2 = 0.01) and, for the numeric cases,
    SNTM_FLUX_PER_EDGE=1, into a separate output folder (outputs_consistent/).

    python make_consistent_manifests.py
"""
import csv
from pathlib import Path

HERE = Path(__file__).resolve().parent
COLUMNS = ["case_id", "method", "setting_id", "group", "seed", "subset",
           "gamma1", "lambda1", "delta", "epsilon", "init_rescale", "lambda2"]
SUBSET = "Hippocampus+PC+RSP"
SEEDS = ["ECl", "CA1", "RSPd"]
SURROGATES = ["mlp", "sr", "gp", "lr"]


def write(name, rows):
    with open(HERE / "manifests" / name, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS)
        w.writeheader()
        w.writerows(rows)
    print(f"manifests/{name}: {len(rows)} rows")


with open(HERE / "manifests" / "recovery_truths.csv", newline="") as f:
    truths = [t for t in csv.DictReader(f) if t["convention"] == "l2f"]

rows = {"numeric": [], "surrogate": []}
for seed in SEEDS:
    for t in truths:
        for method in ["numeric"] + SURROGATES:
            rows["numeric" if method == "numeric" else "surrogate"].append({
                "case_id": f"rec_l2f__{method}__{seed}__{t['truth_id']}",
                "method": method,
                "setting_id": f"rec_l2f__{t['truth_id']}__{seed}",
                "group": "consistent_l2f", "seed": seed, "subset": SUBSET,
                "gamma1": t["gamma1"], "lambda1": t["lambda1"],
                "delta": t["delta"], "epsilon": t["epsilon"],
                "init_rescale": 0.01, "lambda2": t["lambda2"]})
write("consistent_l2f_numeric.csv", rows["numeric"])
write("consistent_l2f_surrogate.csv", rows["surrogate"])

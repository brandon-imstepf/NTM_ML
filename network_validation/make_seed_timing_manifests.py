"""Manifests for two sets of network runs for the reviewer responses.

seeds_l2f_numeric.csv, seeds_l2f_surrogate.csv (R1.m1: several seeds)
    the 8 l2f settings of recovery_truths.csv (lambda2 = 0.01, as in the
    surrogates' training data) seeded in CA1 and RSPd; the ECl runs are the
    recovery truths in outputs_recovery/.
timing_fig8.csv (R1.m4: Table 3 with repeats)
    the Fig. 8 setting (gamma1 5e-5, lambda1 = lambda2 = 0.01, delta 10,
    epsilon 50, init_rescale 0.002), REPEATS runs of each method, case ids
    ending in __rep<k> so that each repeat is its own output file.

    python make_seed_timing_manifests.py
"""
import csv
from pathlib import Path

HERE = Path(__file__).resolve().parent
COLUMNS = ["case_id", "method", "setting_id", "group", "seed", "subset",
           "gamma1", "lambda1", "delta", "epsilon", "init_rescale", "lambda2"]
SUBSET = "Hippocampus+PC+RSP"
SEEDS = ["CA1", "RSPd"]
SURROGATES = ["mlp", "sr", "gp", "lr"]
METHODS = ["numeric"] + SURROGATES
REPEATS = 5


def write(name, rows):
    with open(HERE / "manifests" / name, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS)
        w.writeheader()
        w.writerows(rows)
    print(f"manifests/{name}: {len(rows)} rows")


with open(HERE / "manifests" / "recovery_truths.csv", newline="") as f:
    truths = [t for t in csv.DictReader(f) if t["convention"] == "l2f"]

seed_rows = {"numeric": [], "surrogate": []}
for t in truths:
    for seed in SEEDS:
        for method in METHODS:
            seed_rows["numeric" if method == "numeric" else "surrogate"].append({
                "case_id": f"rec_l2f__{method}__{seed}__{t['truth_id']}",
                "method": method,
                "setting_id": f"rec_l2f__{t['truth_id']}__{seed}",
                "group": "seed_l2f", "seed": seed, "subset": SUBSET,
                "gamma1": t["gamma1"], "lambda1": t["lambda1"],
                "delta": t["delta"], "epsilon": t["epsilon"],
                "init_rescale": 0.01, "lambda2": t["lambda2"]})
write("seeds_l2f_numeric.csv", seed_rows["numeric"])
write("seeds_l2f_surrogate.csv", seed_rows["surrogate"])

timing_rows = []
for k in range(1, REPEATS + 1):
    for method in METHODS:
        timing_rows.append({
            "case_id": f"timing__{method}__ECl__g5e-05_l0.01_d10_e50__rep{k}",
            "method": method,
            "setting_id": "timing__g5e-05_l0.01_d10_e50__ECl",
            "group": "timing_fig8", "seed": "ECl", "subset": SUBSET,
            "gamma1": 5e-05, "lambda1": 0.01, "delta": 10, "epsilon": 50,
            "init_rescale": 0.002, "lambda2": 0.01})
write("timing_fig8.csv", timing_rows)

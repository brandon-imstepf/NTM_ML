"""Build the network-validation manifests (R1.m1, R2.M1).

Writes manifests/numeric.csv and manifests/surrogate.csv. Each row is one
network simulation; surrogate rows are matched to numeric rows by setting_id.

Groups
  interior : full factorial inside the training ranges, 3 seeds (24 sims)
  edge     : corners of the training ranges, ECl seed (8 sims)
  ood      : slightly outside the training ranges, ECl seed (8 sims)
"""
import csv
import itertools
from pathlib import Path

SEEDS_ALL = ["ECl", "CA1", "RSPd"]
GROUPS = {
    # group: (gamma1 values, lambda1 values, (delta, epsilon) pairs, seeds)
    "interior": ([5e-4, 5e-3], [0.02, 0.05], [(80, 20), (20, 80)], SEEDS_ALL),
    "edge": ([5e-5, 5e-2], [0.01, 0.1], [(10, 100), (100, 10)], ["ECl"]),
    "ood": ([2.5e-5, 1e-1], [0.008, 0.12], [(8, 110), (110, 8)], ["ECl"]),
}
SURROGATES = ["mlp", "sr", "gp", "interp", "lr"]
SUBSET = "Hippocampus+PC+RSP"
COLS = ["case_id", "method", "setting_id", "group", "seed", "subset",
        "gamma1", "lambda1", "delta", "epsilon"]


def settings():
    for group, (gammas, lambdas, de_pairs, seeds) in GROUPS.items():
        for g, l, (d, e), s in itertools.product(gammas, lambdas, de_pairs, seeds):
            sid = f"{group}__g{g:g}_l{l:g}_d{d:g}_e{e:g}__{s}"
            yield dict(setting_id=sid, group=group, seed=s, subset=SUBSET,
                       gamma1=g, lambda1=l, delta=d, epsilon=e)


def main():
    rows = list(settings())
    outdir = Path(__file__).resolve().parent / "manifests"
    outdir.mkdir(exist_ok=True)
    interior = [r for r in rows if r["group"] == "interior"]
    for name, methods, subset in [("numeric", ["numeric"], rows),
                                  ("surrogate", SURROGATES, rows),
                                  ("surrogate_interior", SURROGATES, interior)]:
        with open(outdir / f"{name}.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=COLS)
            w.writeheader()
            for m in methods:
                for r in subset:
                    w.writerow({**r, "method": m, "case_id": f"{r['setting_id']}__{m}"})
        print(f"{name}.csv: {len(subset) * len(methods)} cases")


if __name__ == "__main__":
    main()

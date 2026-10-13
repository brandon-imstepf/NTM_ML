"""Build the network-validation manifests (R1.m1, R2.M1).

Writes manifests/<name>.csv. Each row is one network simulation; surrogate rows
are matched to numeric rows by setting_id.

Groups (init_rescale 0.01 unless the name says otherwise)
  interior        : full factorial inside the training ranges, 3 seeds (24 sims)
  edge            : corners of the training ranges, ECl seed (8 sims)
  ood             : slightly outside the training ranges, ECl seed (8 sims)
  interior_r0.002 : interior settings at init_rescale 0.002 (24 sims)
  fig8            : the Fig. 8 setting at init_rescale 0.002, as in the paper (1 sim)
  fig8_r0.01      : the Fig. 8 setting at init_rescale 0.01 (1 sim)
"""
import csv
import itertools
from pathlib import Path

SEEDS_ALL = ["ECl", "CA1", "RSPd"]
INTERIOR = ([5e-4, 5e-3], [0.02, 0.05], [(80, 20), (20, 80)], SEEDS_ALL)
FIG8 = ([5e-5], [0.01], [(10, 50)], ["ECl"])
GROUPS = {
    # group: ((gamma1 values, lambda1 values, (delta, epsilon) pairs, seeds), init_rescale)
    "interior": (INTERIOR, 1e-2),
    "edge": (([5e-5, 5e-2], [0.01, 0.1], [(10, 100), (100, 10)], ["ECl"]), 1e-2),
    "ood": (([2.5e-5, 1e-1], [0.008, 0.12], [(8, 110), (110, 8)], ["ECl"]), 1e-2),
    "interior_r0.002": (INTERIOR, 2e-3),
    "fig8": (FIG8, 2e-3),
    "fig8_r0.01": (FIG8, 1e-2),
}
SURROGATES = ["mlp", "sr", "gp", "interp", "lr"]
SUBSET = "Hippocampus+PC+RSP"
COLS = ["case_id", "method", "setting_id", "group", "seed", "subset",
        "gamma1", "lambda1", "delta", "epsilon", "init_rescale"]


def settings():
    for group, ((gammas, lambdas, de_pairs, seeds), rescale) in GROUPS.items():
        for g, l, (d, e), s in itertools.product(gammas, lambdas, de_pairs, seeds):
            sid = f"{group}__g{g:g}_l{l:g}_d{d:g}_e{e:g}__{s}"
            yield dict(setting_id=sid, group=group, seed=s, subset=SUBSET,
                       gamma1=g, lambda1=l, delta=d, epsilon=e, init_rescale=rescale)


def write(outdir, name, methods, rows):
    with open(outdir / f"{name}.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLS)
        w.writeheader()
        for m in methods:
            for r in rows:
                w.writerow({**r, "method": m, "case_id": f"{r['setting_id']}__{m}"})
    print(f"{name}.csv: {len(rows) * len(methods)} cases")


def main():
    rows = list(settings())
    outdir = Path(__file__).resolve().parent / "manifests"
    outdir.mkdir(exist_ok=True)
    first = [r for r in rows if r["group"] in ("interior", "edge", "ood")]
    rescale = [r for r in rows if r["group"] in ("interior_r0.002", "fig8", "fig8_r0.01")]
    write(outdir, "numeric", ["numeric"], first)
    write(outdir, "surrogate", SURROGATES, first)
    write(outdir, "surrogate_interior", SURROGATES, [r for r in first if r["group"] == "interior"])
    write(outdir, "rescale_numeric", ["numeric"], rescale)
    write(outdir, "rescale_surrogate", SURROGATES, rescale)


if __name__ == "__main__":
    main()

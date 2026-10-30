"""Build the network parameter-recovery manifests.

Each row is one surrogate network simulation at a candidate parameter value. The
loss against the numeric truth is computed afterwards (analyze_network_recovery.m).

Truths: the interior settings, ECl seed. Two lambda2 conventions, as separate manifests:
  l2eq : lambda2 = lambda1 (also while lambda1 is swept). The truths are the existing
         numeric runs in outputs/ (network validation), so no new numeric runs.
  l2f  : lambda2 = 0.01 (the flux convention of the surrogates' training data). The
         truths are new numeric runs, manifests/recovery_truth_l2f.csv.

Candidates: each parameter is swept over a fixed grid across its training range
(log-spaced for gamma1 and lambda1, linear for delta and epsilon), plus every true
value of that parameter, with the other three at their true values. case_id is
built from the parameter values, so a simulation shared by several sweeps or
truths (e.g. the true point, or two truths differing only in the swept parameter)
appears once and runs once.

Writes manifests/
  recovery_truths.csv                  truth_id, convention, values, where its numeric run is
  recovery_truth_l2f.csv               numeric truths for l2f (to run)
  recovery_truth_mlp_l2f.csv           numeric and MLP runs at each l2f truth (for the
                                       time-step test, run with SNTM_T_END / SNTM_DT_MAX)
  recovery_pilot_<conv>.csv            MLP, first truth, gamma1 sweep
  recovery_single_<conv>.csv           all single-parameter sweeps
  recovery_joint_<a>_<b>_<conv>.csv    2-D grid over one pair, the other two at truth
"""
import csv
import itertools
from pathlib import Path

import numpy as np

SEED = "ECl"
SUBSET = "Hippocampus+PC+RSP"
INIT_RESCALE = 1e-2
LAMBDA2_FIXED = 0.01
CONVENTIONS = ["l2eq", "l2f"]
MODELS = ["mlp", "sr", "gp", "lr"]

# truths: the interior settings of make_manifest.py
TRUTHS = [dict(gamma1=g, lambda1=l, delta=d, epsilon=e)
          for g, l, (d, e) in itertools.product([5e-4, 5e-3], [0.02, 0.05], [(80, 20), (20, 80)])]

# search ranges = training ranges; log-spaced or linear
PARAMS = {
    "gamma1": ([5e-5, 5e-2], True),
    "lambda1": ([0.01, 0.1], True),
    "delta": ([10, 100], False),
    "epsilon": ([10, 100], False),
}
N_SINGLE = 15  # grid points per single sweep (plus the true values)
N_JOINT = 9    # grid points per parameter in a joint grid (plus the true values)

COLS = ["case_id", "method", "setting_id", "group", "seed", "subset",
        "gamma1", "lambda1", "delta", "epsilon", "init_rescale", "lambda2"]


def rounded(x):
    return float(f"{x:.4g}")


def grid(name, n):
    (lo, hi), log = PARAMS[name]
    values = np.logspace(np.log10(lo), np.log10(hi), n) if log else np.linspace(lo, hi, n)
    true_values = [t[name] for t in TRUTHS]
    return sorted({rounded(v) for v in values} | {rounded(v) for v in true_values})


def truth_id(t):
    return f"g{t['gamma1']:g}_l{t['lambda1']:g}_d{t['delta']:g}_e{t['epsilon']:g}"


def lambda2_of(theta, conv):
    return theta["lambda1"] if conv == "l2eq" else LAMBDA2_FIXED


def row(theta, method, conv):
    return dict(case_id=f"rec_{conv}__{method}__{SEED}__{truth_id(theta)}",
                method=method, setting_id=f"rec_{conv}__{truth_id(theta)}", group=f"rec_{conv}",
                seed=SEED, subset=SUBSET, init_rescale=INIT_RESCALE,
                lambda2=lambda2_of(theta, conv), **theta)


def single_rows(conv, methods, truths, names):
    for method, t, name in itertools.product(methods, truths, names):
        for v in grid(name, N_SINGLE):
            yield row({**t, name: v}, method, conv)


def joint_rows(conv, methods, truths, a, b):
    for method, t in itertools.product(methods, truths):
        for va, vb in itertools.product(grid(a, N_JOINT), grid(b, N_JOINT)):
            yield row({**t, a: va, b: vb}, method, conv)


def write(outdir, name, rows):
    unique = list({r["case_id"]: r for r in rows}.values())
    with open(outdir / f"{name}.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLS)
        w.writeheader()
        w.writerows(unique)
    print(f"{name}.csv: {len(unique)} cases")


def main():
    outdir = Path(__file__).resolve().parent / "manifests"
    outdir.mkdir(exist_ok=True)

    # where each truth's numeric run is (folder relative to network_validation/)
    with open(outdir / "recovery_truths.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["truth_id", "convention", "gamma1", "lambda1", "delta", "epsilon",
                    "lambda2", "seed", "truth_case_id", "truth_dir"])
        for conv, t in itertools.product(CONVENTIONS, TRUTHS):
            if conv == "l2eq":
                case_id = f"interior__{truth_id(t)}__{SEED}__numeric"
                folder = "outputs"
            else:
                case_id = row(t, "numeric", conv)["case_id"]
                folder = "outputs_recovery"
            w.writerow([truth_id(t), conv, t["gamma1"], t["lambda1"], t["delta"], t["epsilon"],
                        lambda2_of(t, conv), SEED, case_id, folder])
    print(f"recovery_truths.csv: {len(CONVENTIONS) * len(TRUTHS)} truths")

    write(outdir, "recovery_truth_l2f", [row(t, "numeric", "l2f") for t in TRUTHS])
    write(outdir, "recovery_truth_mlp_l2f",
          [row(t, m, "l2f") for m in ["numeric", "mlp"] for t in TRUTHS])
    for conv in CONVENTIONS:
        write(outdir, f"recovery_pilot_{conv}", single_rows(conv, ["mlp"], TRUTHS[:1], ["gamma1"]))
        write(outdir, f"recovery_single_{conv}", single_rows(conv, MODELS, TRUTHS, list(PARAMS)))
        for a, b in itertools.combinations(PARAMS, 2):
            write(outdir, f"recovery_joint_{a}_{b}_{conv}", joint_rows(conv, MODELS, TRUTHS, a, b))


if __name__ == "__main__":
    main()

"""Manifests for the second-network test (R2.M1): the Primary Parameter Set
(gamma1 0.005, lambda 0.02, delta 20, epsilon 80) on a second mouse subnetwork of
about the same size as Hippocampus+PC+RSP (30 regions, 276 directed edges):
connectome_subset "Amygdala+CTXsp", the amygdala and cortical subplate of both
hemispheres (28 regions, 166 directed edges, no region shared with
Hippocampus+PC+RSP, all reachable from the seed), seeded in the left basolateral
amygdalar nucleus (seed code BLA, sntm_config.m). Numeric NTM and the five
surrogates, in the convention of the consistent runs (outputs_consistent/):

    subnet_amyg_numeric.csv    1 numeric case
    subnet_amyg_surrogate.csv  5 surrogate cases (lr, sr, interp, gp, mlp)

The subset needs the 'Amygdala+CTXsp' case in adjacency_cached.m and
NetworkTransportModel.m of the simulator (Tau_Transport-jtorok_dev, 11/17).
Run on ROMULUS (from network_validation/) with
    SNTM_LAMBDA2_W1=lambda1 SNTM_FLUX_PER_EDGE=1 SNTM_OUT_DIR=$PWD/outputs_subnet \
        sbatch --array=1 run_cases.slurm manifests/subnet_amyg_numeric.csv 1
    SNTM_LAMBDA2_W1=lambda1 SNTM_OUT_DIR=$PWD/outputs_subnet \
        sbatch --array=1-5 run_cases.slurm manifests/subnet_amyg_surrogate.csv 1

    python make_scale_manifest.py
"""
import csv
from pathlib import Path

HERE = Path(__file__).resolve().parent
COLUMNS = ["case_id", "method", "setting_id", "group", "seed", "subset",
           "gamma1", "lambda1", "delta", "epsilon", "init_rescale", "lambda2"]
SUBSET = "Amygdala+CTXsp"
SEED = "BLA"
P = dict(gamma1="0.005", lambda1="0.02", delta="20", epsilon="80", init_rescale="0.01", lambda2="0.01")
SETTING = "g0.005_l0.02_d20_e80"
SURROGATES = ["lr", "sr", "interp", "gp", "mlp"]


def row(method):
    return dict(case_id=f"subnet__{method}__{SEED}__{SETTING}__amyg", method=method,
                setting_id=f"subnet__{SETTING}__{SEED}__amyg", group="subnet_amyg",
                seed=SEED, subset=SUBSET, **P)


def write(name, rows):
    with open(HERE / "manifests" / name, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS)
        w.writeheader()
        w.writerows(rows)
    print(f"manifests/{name}: {len(rows)} rows")


write("subnet_amyg_numeric.csv", [row("numeric")])
write("subnet_amyg_surrogate.csv", [row(m) for m in SURROGATES])

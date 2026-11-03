# Network and single-edge checks for the reviewer responses

Dated log. Numbers are copied from the CSVs next to this file; the scripts that
write them are in `network_validation/`. Truths: the 8 l2f settings of
`manifests/recovery_truths.csv` (lambda2 = 0.01), seed ECl, 180 days, 30 regions.

## 2026-11-02

### Coverage of the training data (R2.M2) - `edge_inputs_l2f.m` -> `edge_coverage_l2f.csv`
- 44,727 single-edge inputs per numeric simulation (298 edges x 151 time points, edges with
  N_i > 0 or N_j > 0).
- Training N (NRow, NCol of data_bias_e5_*.csv, 2.2e6 values): min 6.6e-7, 99.9th pct 4.6e-3, max 8.9e-3.
- No network N above the training max; largest network N is 15-40% of it.
- 0.21% of inputs have one node empty (N = 0); 0.01-0.48% of nonzero N are below the training min.

### GP predictive s.d. (R1.M1) - `edge_inputs_l2f.m` -> `gp_variance_l2f.csv`; `gp_trajectory_l2f.m` -> `gp_trajectory_l2f.csv`
- At numeric-run inputs: median s.d. 1.2-3.7e-10, held-out training-distribution rows 2.3e-10.
- Along the GP sNTM's own trajectory (every 5th time point):
  - 5 stable runs (NMSE 5.5e-4 to 0.18): largest s.d. never above its day-0 value.
  - 3 failing runs (NMSE 3.1, 779, non-finite): largest s.d. up 1.7-7.3x by day 2.3.
  - NMSE-779 run (g0.005 l0.02 d80 e20): max N reaches 3.95x the training max by day 22,
    s.d. 14x its day-0 value.

### Mass balance (R2.M4) - `mass_balance_l2f.m` -> `mass_balance_l2f.csv`
- The flux terms of the N update cancel over the nodes for any flux matrix, so node + edge mass
  changes only through the explicit step in M. Drift = cumulative sum_i V_i [d(N+M) - (1+m) dN],
  as a fraction of initial node mass.
- Max drift over 180 d, median [min, max]: numeric 0.011 [0.0098, 0.023]; mlp 0.011 [0.0096, 0.018];
  sr 0.013 [0.004, 0.21]; gp 0.015 [0.011, 165]; lr 0.0095 [0.0067, 0.013].
- |node mass - numeric| / initial: mlp 0.042 [0.016, 0.29]; sr 0.16; gp 0.46; lr 0.16.
- Node mass at 180 d / initial (numeric): 0.73 [0.55, 1.71]. Above 1 means the linearized edge
  mass fell; not investigated.

### Surrogates on the locked single-edge sets (R2.M1) - `single_edge_error_locked.m` -> `single_edge_error_locked.csv`
NMSE over the 32 x 32 (N1, N2) grid, median [min, max]:

| group | MLP flux | MLP W1 | GP flux | SR flux | LR flux |
|---|---|---|---|---|---|
| interior (8) | 2.7e-5 [6.5e-6, 1.4e-4] | 1.1e-3 [2.2e-4, 0.016] | 0.015 | 0.12 | 0.25 |
| edge (8) | 1.7e-4 [8.4e-6, 0.015] | 0.047 [4.3e-4, 0.72] | 0.012 | 0.23 | 0.20 |
| ood (8) | 1.6e-3 [6.3e-4, 0.27] | 0.27 [5.0e-3, 9.3] | 0.26 | 0.11 | 10.5 |
| Fig. 8 params | 5.8e-3 | 3.2e-3 | 0.050 | 0.74 | 1.08 |

### Edge error at the numeric network states - `edge_inputs_l2f.m` -> `edge_error_l2f.csv`
- MLP flux NMSE against the network numeric F: 0.057 [0.014, 0.30]; W1 0.30 [0.14, 1.45].
- NOT a measure of MLP error. Check (scratch, g0.0005 l0.02 d80 e20, 40 network states drawn
  with weight F^2), all against the network's stored F:
  - single-edge numeric at 1e-10 (the training-data generator): NMSE 0.39
  - single-edge at 1e-4/1e-6: 0.26
  - MLP: 0.26
  - MLP vs single-edge 1e-10: 1.7e-5 (W1: 2.4e-5)
- The network flux is reproducible (recomputed at h = 57: 1e-23 vs stored) and barely changes
  at 1e-10 tolerances (NMSE 1.1e-4).
- On some edges the network solver and the single-edge solver disagree by 5-15x at the same
  (N_i, N_j). Example: edge (4,15), h = 57, network 1.21e-8 at either tolerance, single edge 1.19e-9.
- Same recomputation with the unedited simulator (`Tau_Transport-jtorok_dev_orig`): identical
  (edge (4,15) 1.21e-8 at both tolerances; stored vs recomputed NMSE 1.2e-23; 213 s and 530 s
  against 6 s and 9 s with the edited code). The disagreement is in the original code, not the
  11/1 speed-up.
- Open: which solution is right. Possibly the joint fsolve over edges sharing a target converges
  to a different root.

## 2026-11-03: seeds CA1 and RSPd, Table 3 repeats

### Network accuracy, 8 settings x 3 seeds (R1.m1) - `network_checks_l2f.m` -> `network_checks_l2f.csv`; `plot_fig_sy.m` -> `outputs_nm/fig_sy.png`
- NMSE, median [min, max] over 24 simulations:
  - MLP 9.7e-3 [4.9e-4, 0.065], R^2 0.984 [0.755, 0.999]
  - GP 0.10 (22/24 finite)
  - LR 0.11
  - SR 0.17 (23/24 finite)
- MLP median NMSE per seed: ECl 2.4e-3, RSPd 5.4e-3, CA1 0.011. MLP R^2 >= 0.99 in 10/24
  (ECl 6, RSPd 4, CA1 0).
- MLP was the best surrogate in 20/24 simulations. In the other 4 (all g5e-4, d20 e80), the GP was
  better in 3 and tied in 1.
- Non-finite simulations: GP ECl g0.005 l0.02 d20 e80, GP RSPd g0.005 l0.02 d80 e20,
  SR RSPd g0.005 l0.05 d80 e20.
- The CA1 and RSPd runs used the edited simulator (bit-identical results), so their runtime_s is
  not comparable with the ECl runs.

### Clamp of negative N (R2.M4)
- NetworkTransportModel sets N < 0 to eps after each step. Column `n_clamped` counts values == eps
  (of 30 x 151 per simulation).
- Simulations where it fired / most values in one simulation:
  - numeric 9/24, 5
  - mlp 21/24, 239
  - gp 21/24, 1167
  - lr 24/24, 790
  - sr 24/24, 1664

### Table 3 repeats (R1.m4) - `timing_table3.m` -> `timing_table3.csv`
- Setup: Fig. 8 parameters, original code (`jtorok_dev_orig`), node08 (Xeon Gold 5215),
  R2022b, 1 thread, 5 repeats each.
- Mean ± s.e.m. (speedup):
  - numeric 14254 ± 344 s (3.96 h)
  - MLP 195.6 ± 4.6 s (72.9x)
  - SR 194.7 ± 4.8 s (73.2x)
  - GP 280.9 ± 6.7 s (50.7x)
  - LR 194.5 ± 3.6 s (73.3x)

## Regenerate
From `network_validation/`, in MATLAB:
`network_checks_l2f(); plot_fig_sy(); timing_table3(); edge_inputs_l2f(2000); gp_trajectory_l2f(5); mass_balance_l2f(); single_edge_error_locked();`

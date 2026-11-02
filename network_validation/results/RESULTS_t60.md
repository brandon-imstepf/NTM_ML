# 60-day Nelder-Mead parameter recovery: results log

Results of `run_nm_t60.slurm`, kept in git for the paper (`outputs_nm/` is not
committed). Newest entries first. Every number here comes from the CSVs in this
folder, which `analyze_nm_t60.m` rewrites; the commands to regenerate them are at
the end.

## Files

| File | Contents |
|---|---|
| `nm_t60_mid.csv` | 192 searches started at the middle of the training range |
| `nm_t60_near.csv` | 192 searches started at 0.98 x the true value |

One row per search: 8 truths x 4 parameters x 3 truth grids x 2 truths.

| Column | Meaning |
|---|---|
| `grid` | truth grid: `base` (default trange), `dt0.45`, `dt0.225` (max time step in days) |
| `truth` | `numeric` (the test) or `mlp` (the MLP's own run at the true values: control) |
| `truth_id` | true values, e.g. `g0.0005_l0.02_d80_e20` = gamma1 5e-4, lambda1 0.02, delta 80, epsilon 20 |
| `param` | the parameter searched; the other three are fixed at their true values |
| `true_value`, `recovered`, `ratio` | ratio = recovered / true |
| `nmse` | NMSE at the recovered value: sum((s - n)^2) / sum(n^2) over regions and times |
| `nmse_true` | NMSE of the MLP at the true values against the numeric truth (0 for the control) |
| `exitflag` | fminsearch: 1 = both tolerances met, 0 = hit the 40-simulation limit |
| `n_sims` | MLP simulations in the search |
| `at_bound` | recovered value within 1% of the range of a bound (search coordinates) |
| `grid_ok`, `truth_file_ok` | the search used the intended grid and truth file |
| `sim_time_s` | summed simulation time of the search (tic/toc in MATLAB, not SLURM wall time) |
| `sim_time_mean_s` | per simulation |

## Setup

- MLP surrogate, l2f convention, 60 days (`SNTM_T_END=60`), cluster R2022b, 1 CPU per search.
- Search coordinates: log10 for gamma1 and lambda1, linear for delta and epsilon.
  Range: gamma1 5e-5 to 5e-2, lambda1 0.01 to 0.1, delta and epsilon 10 to 100.
- Stopping (`nelder_mead_network_recovery.m`): TolFun 1e-6 (NMSE), at most 40 simulations, and TolX:
  - mid start: 0.1% of the range (0.003 decades for gamma1, about 0.7% in value);
  - near start: one tenth of the first step, 0.0005 decades for gamma1 (about 0.1% in value), 0.1 for delta and epsilon.
- Near start: 0.98 x true; first step 0.005 decades (gamma1, lambda1) or 1 (delta, epsilon).

## 2026-11-02: near start (`nm_t60_near.csv`)

- 192 of 192 searches used the intended grid and truth file. 191 have exitflag 1.
  The exception is delta for g0.0005_l0.05_d80_e20 on dt0.45: exitflag 0 at 99.0, with the loss
  falling toward the bound (100).
- Control (MLP truth): every parameter is recovered within 0.13% (gamma1 and lambda1 within 0.006%).
  Largest final NMSE is 3.5e-7.
- Numeric truth: the three grids agree within 1.5% (gamma1), 8.6% (lambda1), 5.6% (delta), 3.0% (epsilon).
- The recovered NMSE is at most the NMSE at the true values in all 32 dt0.225 searches.

### gamma1, dt0.225, numeric truth

| truth_id | recovered | recovered / true | % error | NMSE at true | NMSE at recovered | exitflag | sims | sim time (s) |
|---|---|---|---|---|---|---|---|---|
| g0.0005_l0.02_d20_e80 | 5.04e-4 | 1.008 | 0.8 | 7.13e-4 | 7.07e-4 | 1 | 10 | 135 |
| g0.0005_l0.02_d80_e20 * | 4.29e-4 | 0.858 | 14.2 | 2.31e-2 | 1.46e-2 | 1 | 17 | 199 |
| g0.0005_l0.05_d20_e80 | 5.30e-4 | 1.059 | 5.9 | 4.27e-3 | 3.82e-3 | 1 | 15 | 190 |
| g0.0005_l0.05_d80_e20 * | 3.41e-4 | 0.682 | 31.8 | 7.06e-2 | 1.75e-2 | 1 | 28 | 396 |
| g0.005_l0.02_d20_e80 | 5.18e-3 | 1.037 | 3.7 | 1.33e-3 | 1.16e-3 | 1 | 13 | 191 |
| g0.005_l0.02_d80_e20 | 4.73e-3 | 0.946 | 5.4 | 2.38e-3 | 1.26e-3 | 1 | 14 | 166 |
| g0.005_l0.05_d20_e80 | 5.00e-3 | 1.001 | < 0.1 (below TolX) | 1.48e-3 | 1.48e-3 | 1 | 9 | 138 |
| g0.005_l0.05_d80_e20 | 4.81e-3 | 0.962 | 3.8 | 1.10e-3 | 7.71e-4 | 1 | 10 | 129 |

\* Numeric truth zigzags in 9-17 regions on the default grid (see Truths below).

- Median error 4.6%, mean 8.2%. Without the two * truths: mean 3.3%, worst 5.9%.
- No gamma1 search is flagged as a local minimum against the loss scan; the scan's best point is
  1.00 x true, except 0.75 x for g0.0005_l0.05_d80_e20.
- All four delta80/epsilon20 settings have recovered NMSE 25-70% of the NMSE at the true values.

### All parameters, numeric truth: recovered / true (base | dt0.45 | dt0.225)

| truth_id | gamma1 | lambda1 | delta | epsilon |
|---|---|---|---|---|
| g0.0005_l0.02_d80_e20 | 0.867 / 0.858 / 0.858 | 0.955 / 0.956 / 0.954 | 1.111 / 1.114 / 1.112 | 0.758 / 0.749 / 0.746 |
| g0.0005_l0.02_d20_e80 | 1.004 / 1.007 / 1.008 | 1.021 / 1.046 / 1.058 | 1.089 / 1.102 / 1.102 | 0.998 / 0.987 / 0.987 |
| g0.0005_l0.05_d80_e20 | 0.681 / 0.681 / 0.682 | 0.885 / 0.885 / 0.885 | 1.238 / 1.238 / 1.237 | 0.758 / 0.758 / 0.771 |
| g0.0005_l0.05_d20_e80 | 1.054 / 1.058 / 1.059 | 1.149 / 1.151 / 1.163 | 1.324 / 1.324 / 1.314 | 0.992 / 0.992 / 0.963 |
| g0.005_l0.02_d80_e20 | 0.949 / 0.946 / 0.946 | 1.584 / 1.497 / 1.499 | 1.147 / 1.131 / 1.132 | 0.855 / 0.861 / 0.861 |
| g0.005_l0.02_d20_e80 | 1.022 / 1.034 / 1.037 | 0.721 / 0.723 / 0.718 | 2.205 / 2.196 / 2.199 | 0.913 / 0.911 / 0.908 |
| g0.005_l0.05_d80_e20 | 0.963 / 0.961 / 0.962 | 1.037 / 1.031 / 1.034 | 1.110 / 1.088 / 1.087 | 0.861 / 0.877 / 0.877 |
| g0.005_l0.05_d20_e80 | 1.000 / 1.000 / 1.001 | 0.784 / 0.797 / 0.808 | 1.546 / 1.602 / 1.596 | 0.949 / 0.939 / 0.941 |

Loss-scan check (default grid): lambda1 for g0.0005_l0.02_d80_e20 (scan best 1.19 x) and
g0.0005_l0.05_d20_e80 (scan best 2.00 x) stopped in local minima. The others reached the scan's minimum.

### Run time (simulation time, cluster, 1 CPU)

| grid | per MLP simulation, median (s) | numeric-truth search, median / max (s) | control search, median (s) |
|---|---|---|---|
| base | 5.2 | 78 / 146 | 46 |
| dt0.45 | 7.5 | 103 / 273 | 81 |
| dt0.225 | 12.7 | 182 / 396 | 121 |

192 searches, 5.8 h of simulation in total. Numeric-truth searches took 9-29 simulations, controls 9-10.

## 2026-11-01: mid start (`nm_t60_mid.csv`)

- 192 of 192 searches used the intended grid and truth file. 190 have exitflag 1. The exceptions:
  - delta for g0.0005_l0.05_d80_e20: exitflag 0 at 98.7, with the loss falling toward the bound;
  - lambda1 for g0.0005_l0.05_d20_e80: stopped at the bound 0.1.
- Control: 94 of 96 within 0.25%. Two local minima: lambda1 for g0.0005_l0.02_d80_e20 (1.71 x) and
  epsilon for g0.005_l0.02_d20_e80 (0.75 x).
- Numeric truth: gamma1 within 6% in 6 of 8 settings (worst 0.68 x); epsilon 0.56-1.0 x;
  delta overestimated in all (1.09-2.2 x); lambda1 0.55-2.0 x, with a nearly flat loss.
  The grids agree within 0.085 in the ratio.
- Loss scan: lambda1 for g0.005_l0.02_d20_e80 stopped in a local minimum (1.81 x; scan best 0.73 x).
- Run time: 7.5 h of simulation in total; per simulation, median 5.0 s (base), 7.3 s (dt0.45),
  12.2 s (dt0.225); 10-31 simulations per search (median 16-17).

## Truths (`outputs_t60_<grid>/`, made by `run_cases.slurm`)

- Default grid: steps of 0.45 d (days 0-18), 0.9 d (18-54), 1.8 d (54 on); 84 time points to day 60.
  dt0.45: 133 points. dt0.225: 265.
- Repeats on different nodes are bit-identical. The 60-day truths equal the 180-day runs cut at day 59.
- Time-step error (NMSE dt0.45 vs dt0.225) is at most 5e-5, against MLP error 7e-4 to 7e-2.
- 6 of 8 numeric truths are smooth. g0.0005_l0.02_d80_e20 and g0.0005_l0.05_d80_e20 zigzag in 9-17
  regions on the default grid (13-17% of the region's peak; worst region Dentate gyrus RH).
  The zigzag shrinks to a few percent on dt0.225; the MLP has none. Both still have a kink near
  day 15 / day 9 on every grid.

## Regenerate

From `network_validation/` (truths in `outputs_t60_*`, searches in `outputs_nm/`):

```
# cluster
sbatch --array=1-32 run_nm_t60.slurm <grid> <numeric|mlp> [near]   # grid: base, dt0.45, dt0.225
```

```matlab
% local, after copying outputs_nm/t60_* back
analyze_nm_t60();            % mid  -> results/nm_t60_mid.csv,  outputs_nm/t60_figures/
analyze_nm_t60([], 'near');  % near -> results/nm_t60_near.csv, outputs_nm/t60_figures_near/
scan_loss_t60('gamma1')      % loss scan, default grid (also lambda1, delta, epsilon)
plot_t60_fits()              % tau curves, dt0.225, mid start
check_t60_runs               % truth repeats, truncation, time step
```

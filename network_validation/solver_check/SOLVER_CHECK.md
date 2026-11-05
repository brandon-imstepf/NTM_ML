# Numeric network flux solver: unconverged solves

Investigation log, started 2026-11-03. Everything here can be rerun from this folder
(commands at the end).

## Summary

The numeric NTM's network flux calculator (`NetworkFluxCalculator.m`, original code
`Tau_Transport-jtorok_dev_orig`) solves the edges into each region in one joint `fsolve` call.
It never checks `fsolve`'s exit flag.

As a simulation proceeds, a growing share of these joint solves stop at `fsolve`'s
iteration/evaluation limit without converging (exit flag 0). The unconverged fluxes are
stored and used in the nodal update. On those edges, the stored "numeric" flux can be an
order of magnitude away from the steady-state flux.

Solving each edge alone (as the single-edge code that made the surrogates' training data does)
converges on almost every edge, to a root with a residual near zero and a non-negative
concentration profile. So the training data, and the MLP that learned them, have the right
flux on these edges; the network numeric runs do not.

The issue is in the original code; the 11/1 speed-up reproduces it exactly.

Findings, as of the end of this log:

- **Where it happens:** the failures occur at γ = 5e-4, δ = 80, ε = 20 (λ 0.02 and 0.05), for all
  three seeds (heavily for ECl and CA1, mildly for RSPd). The other six settings have none.
- **The oscillation:** these failures cause the oscillation previously attributed to the time
  step.
- **A second mismatch, in Δτ:** the training data computed Δτ (W1) with λ₂ = λ₁; the network runs
  use λ₂ = 0.01.
- **Against a numeric model matching the training data in both respects,** the MLP sNTM's error
  at the two affected settings falls from 0.022 and 0.046 to 7.6e-5 and 2.5e-4 (sections 6-9).

## 1. Symptom (2026-11-02)

- At network states of the numeric run g0.0005 l0.02 d80 e20 (seed ECl), the MLP's φ_τ matched the
  single-edge numeric code at 1e-10 tolerances (NMSE 1.7e-5, 40 states drawn with weight F²).
- Against the network's stored flux at the same states:
  - MLP: NMSE 0.26
  - single-edge code at 1e-10 tolerances: 0.39
  - single-edge code at the network tolerances: 0.26
- Worst case, edge 4 -> 15 at step 57: stored 1.21e-8, single-edge 1.19e-9.
- Recomputing the network flux at that state reproduces the stored value:
  - with the edited code: NMSE 1.2e-23 against stored;
  - with the original code: identical (213 s);
  - at 1e-10 tolerances: still 1.21e-8 (NMSE 1.1e-4 against the loose solve).

  So the edits did not cause this, and it is not a tolerance effect. Details are in
  `../results/RESULTS_checks.md`.

## 2. Inputs are mapped correctly

For edge r -> i (column i is the target region), the code passes x0 = N_r (row r of
`N .* Adj`) and xL = N_i. This is the same convention the surrogates use.

## 3. Locating the failure: edges into region 15, step 57 (`probe_target.m`)

`nfc_probe.m` is a copy of the ORIGINAL `NetworkFluxCalculator.m` with three extra options:

- solve only one target region;
- solve each edge alone instead of jointly;
- set how many components of the stacked ODE state `ode45` keeps non-negative.

Background on the last option: the original sets `'NonNegative', 1:length(n0)` with `n0` a
function handle, so `length(n0)` is 1. Only the first component of the stacked state is ever
kept non-negative.

| pos | source | stored | joint (orig.) | per edge | joint, all non-neg. | single-edge code | residual at joint, unclipped | min n at joint, unclipped | min n at per edge |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 4 | 1.21e-8 | 1.21e-8 | 1.19e-9 | 1.19e-9 | 1.19e-9 | -5.1e-3 | -4.7e-3 | 2.1e-4 |
| 2 | 10 | 3.72e-10 | 3.72e-10 | 3.71e-10 | 3.72e-10 | 3.71e-10 | 2.7e-20 | 2.0e-4 | 2.0e-4 |
| 3 | 25 | 2.29e-9 | 2.29e-9 | 2.30e-9 | 2.30e-9 | 2.30e-9 | 2.1e-5 | 2.5e-4 | 2.3e-4 |
| 4 | 30 | 1.06e-9 | 1.06e-9 | 1.07e-9 | 1.06e-9 | 1.07e-9 | -1.5e-18 | 2.1e-4 | 2.1e-4 |

(N_target = 1.9e-4.)

## 4. The joint value is not a root (`probe_residuals.m`)

Edge 4 -> 15, evaluated alone:

| flux A | non-negativity | residual f_ss | min n along edge |
|---|---|---|---|
| 1.211e-8 (joint) | off | -5.06e-3 | -4.66e-3 |
| 1.187e-9 (per edge) | off | -9.2e-7 | 2.09e-4 |
| 1.211e-8 (joint) | on | -3.94e-4 | 0 (clipped) |
| 1.187e-9 (per edge) | on | -9.2e-7 | 2.09e-4 |

- **Original joint solve of region 15:** `fsolve` exit flag 0 (stopped at its iteration or
  function-evaluation limit) after 88 iterations. Final residuals:
  [-3.9e-4, 2.7e-20, 2.1e-5, -1.5e-18].
- **The edge-4 residual of -3.9e-4 is about twice the target concentration.** The joint solve
  converged for edges 2 and 4 and stopped short for edges 1 and 3.
- **The per-edge flux is a genuine steady state:** residual ~1e-6 and a positive profile.

Likely contributors (not separated):

- **Step size:** `fsolve`'s default finite-difference step is sqrt(eps)·max(|x|, 1) ≈ 1.5e-8,
  larger than the fluxes themselves (~1e-9).
- **Clipping:** when the profile touches zero, the residual is non-smooth.
- **One stopping test for the whole region:** an unconverged edge is accepted together with the
  converged ones.

## 5. How widespread (`flux_solve_survey.m`)

Every 5th time step of the 8 l2f numeric runs (seed ECl, 180 days). The stored state N is
re-solved jointly (as the simulator does) and per edge.

31 states per run; `../results/flux_solve_survey_ECl.csv`:

| setting | edges in unconverged joint solves, mean / max | per-edge solves unconverged (all states) | NMSE joint vs per-edge flux, median / max | max \|residual\| / max N, joint / per edge |
|---|---|---|---|---|
| g5e-4 λ0.02 δ80 ε20 | 54% / 77% | 10 | 0.69 / 1.28 | 0.86 / 0.007 |
| g5e-4 λ0.05 δ80 ε20 | 62% / 88% | 2 | 0.22 / 0.28 | 0.74 / 0.007 |
| g5e-4 λ0.02 δ20 ε80 | 0 | 0 | 1.5e-5 / 8.9e-5 | 0.004 / 0.004 |
| g5e-4 λ0.05 δ20 ε80 | 0 | 0 | 1.3e-5 / 1.5e-4 | 0.005 / 0.005 |
| g5e-3 λ0.02 δ80 ε20 | 0 | 4 | 6.8e-5 / 1.6e-4 | 0.009 / 0.009 |
| g5e-3 λ0.02 δ20 ε80 | 0 | 0 | 3.9e-5 / 2.0e-4 | 0.006 / 0.006 |
| g5e-3 λ0.05 δ80 ε20 | 0 | 0 | 1.1e-5 / 8.0e-4 | 0.007 / 0.055 |
| g5e-3 λ0.05 δ20 ε80 | 0 | 0 | 1.6e-5 / 2.2e-4 | 0.005 / 0.010 |

- The re-solved joint flux matches the stored flux (NMSE at most 1.5e-8), so the survey
  reproduces what the simulations used.
- "Edges in unconverged joint solves" counts every edge of a region whose joint fsolve returned
  exitflag <= 0, including edges that individually converged.
- The failures are confined to the two settings with slow aggregation and anterograde bias
  (γ = 5e-4, δ = 80, ε = 20). In both, the share grows over the simulation; for λ0.02:
  - step 6: 0 edges
  - step 11: 18
  - step 31: 104
  - step 91: 212 of 276
- These are the same two settings where:
  - the numeric runs oscillate on the default time grid (`../results/RESULTS_t60.md`);
  - the MLP sNTM has its two largest errors for seed ECl (`../results/RESULTS_checks.md`).
CA1 and RSPd seeds (`../results/flux_solve_survey_CA1_RSPd.csv`). The same two settings fail;
every other setting has 0 unconverged joint solves:

| seed, setting | edges in unconverged joint solves, mean / max | NMSE joint vs per-edge flux, median / max |
|---|---|---|
| CA1, g5e-4 λ0.02 δ80 ε20 | 32% / 64% | 1.0 / 2.8 |
| CA1, g5e-4 λ0.05 δ80 ε20 | 43% / 69% | 0.64 / 0.84 |
| RSPd, g5e-4 λ0.02 δ80 ε20 | 5% / 21% | 2.4e-4 / 1.0 |
| RSPd, g5e-4 λ0.05 δ80 ε20 | 7% / 26% | 6.6e-5 / 0.35 |

The MLP sNTM's largest error over all 24 simulations is at CA1, g5e-4 λ0.05 δ20 ε80, where the
joint solves all converged. That error is the surrogate's own.

`../results/flux_solve_survey.csv` is the ECl survey (a copy of `flux_solve_survey_ECl.csv`).

## 6. Corrected numeric runs

`flux_per_edge` option, added 2026-11-03 to the edited simulator (`jtorok_dev`; the `_orig`
copy is untouched):

- `NetworkFluxCalculator.m`: when set, one scalar `fsolve` per edge in place of the joint
  solve per region. The default (0) is the original code path, unchanged.
- `NetworkTransportModel.m`: passes it through to both flux calls.
- In NTM_ML, `sntm_config.m` / `sntm_simulate.m` set it from `SNTM_FLUX_PER_EDGE=1`. It is
  passed only when set, so runs with the original code are unaffected.

Rerun of the two affected numeric truths (seed ECl, 180 days) into `../outputs_fixedflux/`:

```powershell
$env:SNTM_FLUX_PER_EDGE = '1'; $env:SNTM_OUT_DIR = "$PWD\outputs_fixedflux"
matlab -batch "run_sntm_cases('manifests/recovery_truth_l2f.csv', 1, 1)"   # g5e-4 l0.02 d80 e20
matlab -batch "run_sntm_cases('manifests/recovery_truth_l2f.csv', 3, 3)"   # g5e-4 l0.05 d80 e20
matlab -batch "run_sntm_cases('manifests/seeds_l2f_numeric.csv', 1, 2)"   # CA1, RSPd: g5e-4 l0.02 d80 e20
matlab -batch "run_sntm_cases('manifests/seeds_l2f_numeric.csv', 5, 6)"   # CA1, RSPd: g5e-4 l0.05 d80 e20
```
Then `compare_fixedflux` (this folder) writes `../results/fixedflux_compare.csv`.

Results, seed ECl. NMSE over all regions and time points; each per-edge run took about 50 min
on the PC:

| setting | original vs corrected numeric | MLP vs original | MLP vs corrected | GP vs corrected | SR vs corrected | LR vs corrected |
|---|---|---|---|---|---|---|
| g5e-4 λ0.02 δ80 ε20 | 0.036 | 0.022 | 0.024 | 0.080 | 0.60 | 0.25 |
| g5e-4 λ0.05 δ80 ε20 | 0.045 | 0.046 | 0.112 | 0.27 | 0.35 | 0.41 |

- The fix changes the numeric trajectories by 4-5%.
- It does not bring the MLP closer to the numeric NTM, so the MLP's error at these settings is
  not explained by the flux failure.
- The negative-N clamp (`N(N<0) = eps`) fired 2 and 0 times in the original runs, 0 times in the
  corrected ones.

CA1 and RSPd (per-edge flux; W1 still at λ₂ = 0.01):

| seed, setting | original vs corrected numeric | MLP vs original | MLP vs corrected |
|---|---|---|---|
| CA1, λ0.02 | 0.014 | 0.012 | 7.2e-4 |
| CA1, λ0.05 | 0.011 | 0.011 | 7.9e-4 |
| RSPd, λ0.05 | 0.0017 | 0.024 | 0.024 |
| RSPd, λ0.02 | 5.0e-4 | 0.0032 | 0.0027 |

- **CA1:** the flux fix alone reduces the MLP's error about 15-fold.
- **RSPd:** few joint solves failed (section 5), so the fix changes little. The remaining MLP error
  is consistent with the W1 convention (section 9), which these runs do not change.
- **Oscillation (`../results/fixedflux_zigzag.csv`):** in the CA1 original runs, 28 and 30 regions
  oscillate (largest jump 63% and 50% of the peak). The corrected runs have small residual wiggles
  (1 and 7 regions with > 10 sign changes, largest jump 0.8% and 0.6% of the peak). The MLP has none.

## 7. The oscillation came from the unconverged solves (`zigzag_fixedflux.m`)

Measured per region between days 18 and 180:

- sign changes of the step-to-step change in N (a smooth curve has a few);
- largest second difference relative to the region's peak.

`../results/fixedflux_zigzag.csv`:

| run (ECl) | regions with > 10 sign changes | max sign changes | max jump / peak |
|---|---|---|---|
| original, λ0.02 | 15 | 85 | 0.70 |
| corrected, λ0.02 | 0 | 7 | 0.007 |
| MLP, λ0.02 | 0 | 1 | 0.007 |
| original, λ0.05 | 13 | 66 | 0.09 |
| corrected, λ0.05 | 0 | 3 | 0.006 |
| MLP, λ0.05 | 0 | 2 | 0.006 |

The zigzag in these two numeric truths, previously put down to the explicit time step
(`../results/RESULTS_t60.md`), is caused by the unconverged flux solves. With per-edge solves
the default time grid gives smooth curves.

## 8. Δτ (W1): a training-data convention, not a solver problem (`w1_check.m`)

At 40 states of the corrected λ0.05 run (drawn with weight W1²), the network W1 was compared
with the single-edge `SingleEdgeWCalculator`, given the same flux and tolerances. NMSE:

| comparison | NMSE |
|---|---|
| network vs single-edge, λ₂ = 0.01 | 2.4e-32 (identical) |
| network vs single-edge, λ₂ = λ₁ | 0.017 |
| MLP vs single-edge, λ₂ = λ₁ | 1.0e-5 |
| MLP vs network | 0.023 |

The mass-balance code agrees with the single-edge code. The difference comes from the training
data:

- Their **flux** call left λ₂ out, so the simulator default 0.01 was used.
- Their **W1** call passed λ₂ = λ₁ (`../../parameter_recovery/locked_edge/make_locked_edge_data.m`
  documents this).
- So the MLP learned W1 for λ₂ = λ₁. Every network run, here and in the paper, computes W1 with
  λ₂ = 0.01. At λ₁ = 0.05 the two differ by 1.7% (NMSE); at λ₁ = 0.02 by 0.5% (2026-11-02
  check, `../results/RESULTS_checks.md`).

## 9. Numeric runs in the training data's convention

The `lambda2_w1` option, added to the edited simulator (`NetworkTransportModel.m`), sets the
mass balance's λ₂ separately; [] (the default) means λ₂, as in the original. In NTM_ML,
`SNTM_LAMBDA2_W1=lambda1` passes λ₁.

Runs into `../outputs_consistent/`, the model the surrogates were trained on (per-edge flux at
λ₂ = 0.01, W1 at λ₂ = λ₁). The MLP is rerun the same way, because its simulation also computes
the edge-mass terms with the mass balance's λ₂.

```powershell
$env:SNTM_LAMBDA2_W1 = 'lambda1'; $env:SNTM_OUT_DIR = "$PWD\outputs_consistent"
matlab -batch "run_sntm_cases('manifests/recovery_truth_mlp_l2f.csv', 9, 9); run_sntm_cases('manifests/recovery_truth_mlp_l2f.csv', 11, 11)"
$env:SNTM_FLUX_PER_EDGE = '1'
matlab -batch "run_sntm_cases('manifests/recovery_truth_l2f.csv', 1, 1)"
matlab -batch "run_sntm_cases('manifests/recovery_truth_l2f.csv', 3, 3)"
```

Results, seed ECl (`compare_consistent.m` -> `../results/consistent_compare.csv`). NMSE; the
second argument is the reference:

| setting | fixed flux vs original | consistent vs fixed flux | MLP vs original (as reported) | MLP vs consistent, both in the training convention |
|---|---|---|---|---|
| g5e-4 λ0.02 δ80 ε20 | 0.039 | 0.012 | 0.022 | **7.6e-5** |
| g5e-4 λ0.05 δ80 ε20 | 0.065 | 0.098 | 0.046 | **2.5e-4** |

- **Against the model it was trained on, the MLP sNTM's network error at these two settings is
  100-300 times smaller.** That is about the level of its best settings
  (`../results/RESULTS_checks.md`).
- The error reported so far comes almost entirely from two properties of the numeric baseline,
  not from the surrogate:
  - **Unconverged flux solves:** they dominate at λ0.02.
  - **The W1 convention:** it dominates at λ0.05, where λ₂ = 0.01 and λ₂ = λ₁ differ most.
- Changing only the MLP run's mass-balance λ₂ moves it by NMSE 0.023 and 0.018 (last column of
  the CSV). The convention matters on both sides.

### Fig. 8 setting

`flux_solve_survey(5, [], '../outputs_timing/timing__numeric__ECl__g5e-05_l0.01_d10_e50__rep1.mat')`
-> `../results/flux_solve_survey_files.csv`. Setting γ 5e-5, λ 0.01, δ 10, ε 50, init 0.002,
original code:

- 31 states;
- no unconverged joint solve and no unconverged per-edge solve;
- stored vs re-solved flux NMSE ≤ 2e-27;
- joint vs per-edge flux NMSE ≤ 0.011, with both converged; not investigated further.

Since λ₂ = λ₁ = 0.01 here, neither problem applies to Fig. 8 or Table 3.

## 10. The new options leave the default path unchanged; the stored truths come from an older code state

`check_defaults.m` and one-step runs, 2026-11-03, setting g5e-4 λ0.02 δ80 ε20, seed ECl:

- **The new options change nothing when off.** On the PC, the edited simulator (both options off)
  matches the original `_orig` code: max |ΔN| = 1.4e-20 after one step, W1 identical.
- **Earlier cluster runs agree with the PC.** The 11/01 cluster speed-check runs (node08,
  original and sped-up code) and the PC run differ from the stored truth by the same amount:
  - max |ΔN| = 2.2e-4 over the first 5 time points;
  - F and W1 at t = 0 identical.
- **The stored truths are the odd ones out.** They are in `outputs_recovery/`, from node07 on
  10/19-20.
- **A 10/29 cluster rerun of this case differs from the stored truth by NMSE 0.029 over 180 days**
  (`outputs_check/`).
- **The well-behaved setting g5e-3 λ0.02 δ20 ε80 also differs after one step (3.6e-5).** So the
  cluster's simulator copy changed between 10/20 and 10/29, before the 11/01 speed-up. What
  changed is not known.
- **Each comparison used one code state.**
  - ECl: the numeric and surrogate runs are all from 10/19-20.
  - CA1/RSPd: all from 11/02.
  - The reruns in section 11 use one code state for everything.

## 11. Option (a): all R1.m1 simulations in the training convention (decided 2026-11-03)

- **What is rerun:** 24 numeric cases (8 settings x 3 seeds, per-edge flux, W1 at λ₂ = λ₁) and
  96 surrogate cases (W1 path at λ₂ = λ₁).
- **Manifests:** `../manifests/consistent_l2f_{numeric,surrogate}.csv`, from
  `../make_consistent_manifests.py`.
- **Where:** on the cluster, into `outputs_consistent/`.
- **Simulator:** the cluster copy of `jtorok_dev` needs the edited `NetworkTransportModel.m` and
  `NetworkFluxCalculator.m` (backup taken first).
- The surrogate runs ignore `SNTM_FLUX_PER_EDGE` (no numeric flux), but it does no harm to leave
  it unset for them.

## 12. Results of the option (a) runs (2026-11-04)

- **Runs:** 120 on the cluster (24 numeric, 96 surrogate), copied to `../outputs_consistent/`; all
  120 present. Analyses with `run_set = "consistent"` (`../l2f_runs.m`); results in
  `../results/*_consistent.csv`, figure `../outputs_nm/fig_sy_consistent.{png,pdf}`.
- **Reproducible on the PC:** the local ECl MLP runs at the two problem settings (section 9) equal
  the cluster runs (NMSE 7.64e-5 and 2.46e-4 both places); one CA1 MLP case rerun on the PC
  equals the cluster run to 3e-32.
- **`lambda2_w1` affects the surrogate runs too.** It sets λ₂ for the whole mass-balance step
  (`MassBalanceCalculator_explicit_W_faster`), which also builds the synaptic-cleft edge profiles
  behind the edge mass v1, v2, not only W1. Same code, same CA1 MLP case
  (g0.005 l0.05 d20 e80): with and without `SNTM_LAMBDA2_W1=lambda1`, NMSE 1.05e-2 between the two;
  without it, the run equals the old `outputs_seeds/` run to 4e-32. So the surrogate runs in
  `outputs_seeds/` differ from the new ones by this option alone (MLP 2e-5 to 1.3e-2). Numeric and
  surrogate runs in `outputs_consistent/` use the same setting, so the comparison is like for like.
- **Network error against the numeric NTM** (24 simulations per method; `network_checks_l2f_consistent.csv`):

  | method | NMSE median [min, max] | R² median [min, max] | non-finite | clamped sims (max values) |
  |---|---|---|---|---|
  | numeric | – | – | 0 | 0/24 (0) |
  | MLP | 3.7e-4 [2.4e-5, 9.3e-3] | 0.999 [0.988, 1.000] | 0 | 21/24 (239) |
  | GP | 6.0e-2 [2.3e-4, 6.0] | 0.883 [-14.3, 1.000] | 1 | 21/24 (651) |
  | LR | 8.6e-2 [1.8e-2, 1.00] | 0.788 [-1.59, 0.970] | 0 | 24/24 (793) |
  | SR | 0.11 [1.2e-2, 5.28] | 0.671 [-6.0, 0.983] | 0 | 24/24 (724) |

  - MLP best in 24/24; R² ≥ 0.99 in 23/24. Largest MLP errors at RSPd, γ 5e-3, λ 0.02.
  - GP error ≥ 0.80 in 4 of the 6 runs at γ 5e-3, δ 80, ε 20 (one non-finite).
  - The numeric NTM no longer needs the clamp (0 of 24; before: 9 of 24, at most 5 values).
- **Time step** (`../time_step_consistent.m`, ECl, 60 days): MLP change between the default grid and
  0.225-day steps 1.3–4.9e-5; MLP vs numeric 1.9e-5–1.4e-3, of the same size as the grid change at
  γ 5e-4, δ 20, ε 80 (both λ).
- **Coverage, edge error, GP variance, mass balance:** `edge_*_l2f_consistent.csv`,
  `gp_variance_l2f_consistent.csv`, `gp_trajectory_l2f_consistent.csv`,
  `mass_balance_l2f_consistent.csv`. The response draft (R1.m1, R2.M2, R2.M4, R1.M1) quotes them.
- **GP variance along the GP runs** (`gp_trajectory_l2f`, now all three seeds): the largest s.d.
  over the edges rose 1.5× (CA1, γ 5e-3, λ 0.02, δ 80, NMSE 6.0) and 2.5× (RSPd, same setting,
  non-finite) above day 0, passing 1.2× by day 2.25; in the other 22 runs at most 1.14×, including
  NMSE 0.92 (ECl) and 0.80 (RSPd). The earlier result (rise in 3 failing ECl runs) came from the
  old runs; no ECl GP run diverges now.
- **Run times** of these runs are from the optimized code (numeric median 50 min, MLP 13 s) and
  are not used in the paper.

## 13. Accuracy near the edge of the training domain (2026-11-05)

The advisor asked whether the GP (and LR) fail near the edge of the domain or where tau
concentrates. Scripts: `domain_edge_error.m` (single edge), `edge_inputs_l2f(2000, "consistent")`
(network: new coverage columns `frac_N_below_p1` and `frac_N_above_p99`; edge error split into
`*_low_N` and `*_other`), `make_table_domain_edge.py` (Supplementary Table S3).

**Single edge.**
- Data: 10⁵ held-out rows (`data_bias_e5_11.csv`). No surrogate was fit on them; the MLP script
  only records its error on them.
- Edge depth: d = min over inputs of min(u, 1−u), where u is the input's quantile among the rows
  of `data_bias_e5_1.csv`.
- NMSE is Eq. 7 within each bin.

| Method | φτ, interior (d ≥ 0.15) | φτ, edge (d < 0.01) | Δτ, interior | Δτ, edge |
|---|---|---|---|---|
| MLP | 4.7e-5 | 2.6e-4 | 1.7e-4 | 6.4e-3 |
| GP | 8.9e-4 | 2.0e-2 | 5.7e-3 | 5.1e-2 |
| LR | 0.31 | 0.40 | 0.54 | 0.54 |
| SR | 0.16 | 0.26 | 0.20 | 0.31 |

- **Upper 1% of N1** (tau concentrated): GP φτ 5.8e-2, 66× its interior; MLP 8.2e-4.
- **GP predictive s.d.:**
  - from interior to edge it rises 2.0× while the RMSE rises 12×;
  - in the upper 1% of N1 it rises 4.7× while the RMSE rises 36×.
- **High-γ tail of Δτ:** NMSE is large for every method (GP 2.6, LR 379), because the bin variance
  is small. The MSE ratio to the interior is only about 2–3×. This tail is not reported.

**Network** (24 consistent runs):
- **Parameter settings:** between the 11th and 78th percentiles of the training distributions.
- **Nonzero N below the training 1st percentile (7.4e-5):**

  | Seeds | slow aggregation (γ 5e-4) | fast aggregation (γ 5e-3) |
  |---|---|---|
  | ECl, CA1 | 3–8% | 17–23% |
  | RSPd | 25–37% | 75–82% |

- **Nonzero N above the 99th percentile:** at most 0.13%.
- **Edge error at the numeric states, low-N edges vs the others** (median over the 24 runs):

  | Method | low-N edges | other edges | higher on low-N edges in |
  |---|---|---|---|
  | GP φτ | 0.43 | 0.078 | 19/24 runs |
  | GP Δτ | 7.3e-2 | 6.0e-3 | |
  | MLP φτ | 4.1e-3 | 2.1e-3 | 18/24 runs |
  | MLP Δτ | 2.9e-4 | 7.0e-5 | |
  | LR φτ | 0.35 | 0.32 | 15/24 runs |
  | SR φτ | 0.93 | 0.58 | 17/24 runs |

**Conclusion:**
- The GP's errors, like the MLP's but larger, are concentrated at the edge.
- LR and SR are poor throughout the domain.
- The written-up version is in the supplement section "Surrogate accuracy near the edge of the
  training domain", the response (R1.M1, R2.M2) and the manuscript (Results and Limitations).

## Consequences (for the reviewer response and the paper)

1. **R1.m1 (network accuracy) and R2.M4 (stability):**
   - The numeric truths at λ₁ ≠ 0.01 compare the surrogates with a model that differs from their
     training data in Δτ (λ₂).
   - At γ = 5e-4, δ = 80, ε = 20 the numeric truths also contain unconverged flux solves.
   - The reported MLP errors there measure these differences, not the surrogate.
2. **Fig. 8 (λ = 0.01):** λ₂ = λ₁ = 0.01, so the W1 convention does not differ there.
   - The Fig. 8 setting (γ 5e-5, δ 10, ε 50) has not been surveyed for unconverged solves.
3. **The oscillation in the default-grid numeric truths** comes from the unconverged solves, not
   the time step.
   - The R1.m1 and R2.M4 response text that says otherwise needs correcting.
4. **Parameter-recovery truths** (`../results/RESULTS_t60.md`) at these two settings have the same
   problems. The recovery biases reported there for them are not the surrogate's alone.
5. **Undecided: which model is "the" NTM.**
   - The training data mix λ₂ = 0.01 (flux) and λ₂ = λ₁ (W1).
   - Either state that in the paper and compare against matching numeric runs, or regenerate the
     training data with one λ₂.

## Rerun

From this folder, in MATLAB, after making the probes in a shell:

```
python -I make_probe.py      <jtorok_dev_orig>/Code/Transport_Network/NetworkFluxCalculator.m nfc_probe.m
python -I make_fast_probe.py <jtorok_dev>/Code/Transport_Network/NetworkFluxCalculator.m      nfc_fast_probe.m
```

```matlab
setenv('NTM_SIM_ROOT', '<jtorok_dev_orig>'); probe_target; probe_residuals
setenv('NTM_SIM_ROOT', '<jtorok_dev>');      S = flux_solve_survey(5);
```

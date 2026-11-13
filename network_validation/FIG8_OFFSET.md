# Fig. 8 offset: why the MLP and GP sNTMs sit at a constant factor of the numeric NTM

Investigation log, started 2026-11-10. Everything here can be rerun from `network_validation/`.
Newest entries at the end. Numbers are copied from the CSVs named in each step.

## The observation

The new Fig. 8 runs (`manifests/fig8_gp.csv` → `outputs_fig8/`; γ 5e-5, λ 0.01 (= λ₂), δ 10, ε 50,
ECl seed, `init_rescale` 0.01, 180 days, default grid, joint flux solve) give, per month, R² ≈ 0.999 for
the MLP and GP sNTMs, but a constant offset. Total soluble tau Σ N, sNTM / numeric:
MLP 0.89, GP 1.23, from about day 2 to day 150. Figure: `gp_manuscript/figures/fig8_gp.png`
(`gp_manuscript/plot_fig8_gp.m`).

The per-month NMSE (Eq. 7) in the figure (MLP 0.07 → 4.3, GP 0.28 → 19) is this offset divided by the
variance across the 30 regions. The MSE stays constant at about 7e-3 nM² (MLP) and 3e-2 (GP), and
99.8 to 100% of it is the offset (scatter about the fitted line is < 0.15% of the MSE). The NMSE grows
over months 1–3 only because the variance across regions falls about 60-fold. (scratchpad check, 11/10)

## Step 1: where the offset appears (existing runs), 11/10

Σ N sNTM / numeric at day 150:

| runs | seed N₀ | MLP | GP |
|---|---|---|---|
| Fig. 8 parameters, `init_rescale` 0.002 (Table 3 timing runs, `outputs_timing/`) | 1.8e-3 µM | 1.02 | 1.03 |
| Fig. 8 parameters, `init_rescale` 0.01 (`outputs_fig8/`) | 7.3e-3 µM | 0.89 | 1.23 |
| 24 interior runs, `init_rescale` 0.01 (`outputs_consistent/`) | 7.3e-3 µM | median 1.007 (0.98–1.16) | median 1.09 (0.68–3.2) |

So the offset needs both the high seed concentration and the Fig. 8 parameters (γ, λ, δ on the lower
bounds of the training range).

## Step 2: error map and mass budget (`fig8_offset_diagnostics.m`), 11/10

Outputs: `results/fig8_offset_errors_by_time.csv`, `results/fig8_offset_seed_edges_t0.csv`,
`results/fig8_offset_mass.csv`, `outputs_fig8/figures/offset_error_map.png`, `offset_seed_mass.png`.

- At day 0 only the 36 seed edges are active. Their inputs are at a corner of the training data:
  seed N = 7.3e-3 µM = 0.82 × the largest training N (above the 99.9th percentile), the other node
  N = 0, and γ, λ, δ on their lower bounds.
- Surrogate vs numeric on these edges at day 0 (all outgoing edges share one input, all incoming
  edges another):

  | | numeric | MLP | GP |
  |---|---|---|---|
  | φτ, seed → neighbour | 9.6e-9 | 6.7e-9 (−30%) | 1.05e-8 (+8%) |
  | φτ, neighbour → seed | −9.35e-9 | −9.48e-9 (+1%) | −1.15e-8 (+23%) |
  | Δτ (W1), seed → neighbour | 8.5e-7 | 8.9e-7 (+5%) | −4.4e-7 (wrong sign) |
  | Δτ (W1), neighbour → seed | 8.0e-7 | 4.6e-7 (−43%) | 6.7e-7 (−16%) |

- The errors persist while the seed is concentrated and fade by days 4–10 (MLP Δτ NMSE on the other
  edges 7.2 at day 0.45, 0.23 at day 4.5). The offset forms in the same window and then freezes.
  NMSE over the 36 seed edges is not meaningful (two distinct inputs); use the per-edge table.
- Node mass Σ V(N + M) is not constant in any run: numeric +29% by day 10, +14% at the end; MLP +14% and
  +2%; GP +65% and +42%. The fluxes cancel when summed over the nodes, so φτ alone cannot change this
  total; it changes through the mass the model holds on the edges (the mass-correction terms built from
  the edge profiles and Δτ). The MLP seed empties slightly faster than numeric while every neighbour
  fills more slowly, so the missing tau goes onto the edges, not to other regions.
- Open: the stored numeric `EdgeMass` (sum about 0.23 against a node mass of 42) does not close the
  budget; it is not the edge mass of the node update in the same units. Do not use it (the red line
  in `offset_seed_mass.png` is meaningless).
- Open: C_ij depends on both the edge profile (from φτ) and Δτ, so step 2 does not separate them.

## Step 3: swap test (`fig8_swap_runs.m`), started 11/10

`sntm_simulate.m` now accepts `use_sur = [use_flux, use_w1]` (a scalar applies to both, so existing
runs are unchanged). Two runs at the Fig. 8 parameters, `init_rescale` 0.01:

- A: numeric φτ, MLP Δτ → `outputs_fig8/fig8__swapA_numflux_mlpw1__ECl__g5e-05_l0.01_d10_e50.mat`
  (fast: about 2 min; the numeric cost is in the W1 / mass-balance step, not the flux solve)
- B: MLP φτ, numeric Δτ → `outputs_fig8/fig8__swapB_mlpflux_numw1__ECl__g5e-05_l0.01_d10_e50.mat`
  (about 1 h). The numeric W1 here is computed from the edge profile of the MLP flux.

Logs: `logs/fig8_swap_A.log`, `logs/fig8_swap_B.log`.

Reading the result: compute Σ N / numeric over time for A and B, as for the MLP run.
If A shows 0.89 and B shows about 1, Δτ is the cause; if the reverse, φτ (through the profile and
C_ij); if both partly, the split.

```
matlab -batch "fig8_swap_runs('A')"     % and 'B', one process each
```

## Possible next steps (not started)

- Starting-concentration sweep at the Fig. 8 parameters: `init_rescale` 0.002, 0.004, 0.006, 0.008,
  0.01 for numeric, MLP, GP. Shows where the offset switches on against the training N distribution.
- Corner test: `init_rescale` 0.01 with γ, λ, δ moved slightly inside their ranges (e.g. γ 1e-4,
  λ 0.012, δ 15). Separates the corner of the parameter range from the high seed concentration.

## Context

- Brandon decided (11/10) to keep the new Fig. 8 as is and set the offset aside for the paper for now.
- Fig. 7 still comes from the original (published) runs, whose seed concentration differs
  (published Fig. 7b seed 1.25e-2 µM); see questions for Nuutti, item 8.

### Step 3 result, run A (11/10)

`fig8_swap_compare.m` → `results/fig8_swap_compare.csv`. Σ N / numeric:

| run | day 0.45 | day 2.25 | day 10 | day 30 | day 150 | slope month 1 / 5 |
|---|---|---|---|---|---|---|
| MLP φτ + MLP Δτ | 0.964 | 0.890 | 0.894 | 0.893 | 0.891 | 0.884 / 0.882 |
| GP | 1.064 | 1.227 | 1.275 | 1.226 | 1.225 | 1.149 / 1.179 |
| A: numeric φτ + MLP Δτ | 0.992 | 0.951 | 0.923 | 0.915 | 0.912 | 0.892 / 0.849 |

With exact flux, the MLP's Δτ alone gives about 80% of the offset (8.8% against 10.9% at day 150),
building up over days 0–10 instead of by day 2. So Δτ is the main cause, and φτ probably adds the fast
early part (to be confirmed by run B, which was still running when this was written).

### Step 3 result, run B (11/10)

Run B finished in 4311 s. `fig8_swap_compare.m` → `results/fig8_swap_compare.csv`. Σ N / numeric:

| run | day 0.45 | day 2.25 | day 10 | day 30 | day 60 | day 150 | slope month 1 … 5 |
|---|---|---|---|---|---|---|---|
| MLP φτ + MLP Δτ | 0.964 | 0.890 | 0.894 | 0.893 | 0.892 | 0.891 | 0.884 0.856 0.867 0.881 0.882 |
| A: numeric φτ + MLP Δτ | 0.992 | 0.951 | 0.923 | 0.915 | 0.913 | 0.912 | 0.892 0.883 0.861 0.849 0.849 |
| B: MLP φτ + numeric Δτ | 0.981 | 0.954 | 0.986 | 0.996 | 0.999 | 0.999 | 1.018 0.996 1.031 1.061 1.064 |

Reading:
- The lasting offset comes from Δτ. With numeric Δτ (B), the total returns to the numeric value (0.999 by day 60);
  with MLP Δτ (A), it freezes at 0.91.
- Early on, both contribute about equally: at day 2.25 the deficits are 4.9% (A) and 4.6% (B), which add up to
  about the full MLP deficit of 11.0%. B's early dip is temporary: the numeric mass balance gives the tau back
  once the seed has spread. A's is not.
- The deficits at day 150 add up to 8.9% against the full 10.9%. The remaining ~2% is interaction: MLP φτ changes
  the edge inputs that MLP Δτ then sees, in the region where Δτ is wrong.
- B's total is right, but its regional slope drifts to 1.06 by month 5, so the MLP φτ moves some tau between
  regions without changing the total. This is a separate, smaller effect.
- Conclusion: the Fig. 8 offset is an error of the Δτ (W1) surrogate at the seed edges in days 0–10, where the
  inputs sit at a corner of the training data (high seed N, other node N = 0, γ, λ, δ on their lower bounds).
  The mass that error places on the edges is never returned to the nodes. φτ errors there are self-correcting.

Figure: `fig8_swap_plot.m` → `outputs_fig8/figures/swap_test.png/.pdf` (Σ N / numeric over time, log days;
full MLP, A, B, and A + B deficits added).

## Fig. 8 parameters on the edge of the training range (11/11)

Confirmed: the manuscript gives γ 5e-5, λ 0.01, δ 10, ε 50 (Results, Fig. 8 caption, Methods), and the
training data (`data_bias_e5_1.csv`) have minima γ 5.0006e-5, λ 0.0100000, δ 10.002. So γ, λ, δ sit on
the lower bounds. Fraction of training samples within 10% of each value: γ 1.4%, λ 4.2%, δ 1.1%. Sampled
independently, that is about 7e-6 for all three at once (under 1 sample per 10^5).
The advice received is that Fig. 8 should not use the edges. Fig. 8-style figures for the 8 interior settings
with the ECl seed (existing `outputs_consistent/` runs) are in `gp_manuscript/figures/fig8_interior/`
(`gp_manuscript/plot_fig8_interior.m`; values in `logs/fig8_interior.log`).

## Fig. 8-style figures for every setting (11/12)

Brandon asked for the Fig. 8 figure (all surrogates, R² and NMSE per month) at every interior, edge and ood
setting. Script: `gp_manuscript/plot_fig8_sets.m` (uses `plot_fig8_gp.m`; replaces the short-lived
`plot_fig8_interior.m`). Figures: `gp_manuscript/figures/fig8_sets/<group>/` and one PDF per group;
numbers: `results/fig8_sets.csv`.

- Interior: the 24 runs of `outputs_consistent/` (8 settings × ECl, CA1, RSPd). Done.
- Edge and ood: the only existing numeric runs (`outputs/`) are in the old convention (λ₂ = λ for the flux) and
  cover 4 of 8 settings each. The surrogate runs do not depend on λ₂ (interior check: identical N in
  `outputs/` and `outputs_consistent/`), so only the numeric runs were redone: 16 runs,
  `manifests/consistent_edge_ood_numeric.csv` (`make_edge_ood_manifest.py`), into `outputs_consistent/`
  with SNTM_LAMBDA2_W1=lambda1 and SNTM_FLUX_PER_EDGE=1, 8 workers; log `logs/consistent_edge_ood.log`.
  Started 11/12. Then: `plot_fig8_sets({'edge','ood'})`.
- 11/12: the 16 edge/ood numeric runs finished (0 failures, 49-65 min each). All four groups are drawn
  (`plot_fig8_sets()`): `fig8.pdf`, `interior.pdf` (24 pages), `edge.pdf`, `ood.pdf` (8 each). Surrogate runs
  that diverged (NaN) are marked "simulation diverged (day d)" and stored as NaN in `results/fig8_sets.csv`.
  Diverged runs: interior GP 1/24; edge MLP 2/8 (γ 0.05, λ 0.01), GP 1/8, SR 2/8, LR 1/8; ood GP 2/8, LR 3/8.
- 11/12: regional trajectory figures (as panel b of the glass-brain figure) for every setting:
  `gp_manuscript/plot_trajectory_sets.m` → `gp_manuscript/figures/trajectory_sets/<group>.pdf` (days 0-60;
  `plot_trajectory_sets(groups, 180)` for the full run). Run sets are now listed once in
  `gp_manuscript/fig8_cases.m`, used by `plot_fig8_sets.m` and `plot_trajectory_sets.m`.

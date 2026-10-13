# Griddatan Surrogate Comparison

## Purpose

A second surrogate-model family (`griddatan`, scattered-data N-D interpolation)
built as an apples-to-apples comparison baseline against the Gaussian Process
pipeline in `../gaussian_process-interpolant/`. The original ask was a
gridded interpolant (`griddedInterpolant`); that requires samples on a full
rectangular grid, which the training data does not satisfy (it's randomly
subsampled scattered points in 6-D), so `griddatan` — the scattered-data
analogue — is the next viable option. See `STATUS.md` for whether it actually
works.

## Scope boundary

Everything this project produces lives in this folder. Files elsewhere
(`../gaussian_process-interpolant/`, `../parameter_recovery/`, shared scaling
`.mat` files, training/validation CSVs) are read-only references, used to
guarantee apples-to-apples comparability — never edited.

## Reproducibility rules

- `training_seed = 42` everywhere a random subset is drawn — identical to the
  GP pipeline, so the two pipelines train on the exact same rows.
- All paths are relative to this folder's own location
  (`fileparts(mfilename('fullpath'))`), not to whatever MATLAB's current
  directory happens to be — see the "files aren't saving where I expect"
  lesson learned on the GP side.
- `run_all.m` regenerates every artifact in `saved_models/`, `figures/`, and
  `results/` from the raw input data alone. No manual steps.
- Shared standardization constants (`../parameter_recovery/flux_scaling_params.mat`,
  `../parameter_recovery/emc_scaling_params.mat`) are loaded, never
  recomputed, so both pipelines standardize identically.

## Dataset-size parity, explicitly

This section has been wrong twice now, corrected each time by direct
evidence rather than a new assumption — see Milestone 1 and Milestone 2 in
`STATUS.md` for the full trail. Current, verified state:

| Folder | Files | Training rows | Provenance verified against current data? |
|---|---|---|---|
| `saved_models/` (unsuffixed) | 1000 (full 5×4×10×5 grid) | 1000 (row count confirmed) | **No** — root cause identified, not just ruled out: trained directly on `data/data_baseline_models/e3_data/bias_train_e3.csv` (exact diff = 0), the pre-100k-superset pipeline. Predates commit `b61734c` (2026-10-02, switch to the 100k-row superset); files are timestamped 2026-09-22. |
| `saved_models_superset_3000/` | 307 | 3000 | **Yes** — exact (diff = 0) match on training rows and on reproduced validation predictions for both flux and emc. |
| `saved_models_superset_5000/` | 10 (`best_models` shortlist) | 5000 | Not checked; orphaned, nothing currently loads it |

Because the 1000-row corpus can't be proven to come from the data this
project actually has access to, it is **not used as a comparison target**.
The parity rule, simplified from the earlier two-size plan:

- **Both the NMSE/model-comparison stage and the parameter-recovery stage
  mirror `dataset_size = 3000`** — the one corpus with proven, exact
  provenance (`saved_models_superset_3000/`, same corpus
  `synthetic_parameter_recovery_gp.m` already uses).
- The GP side's own 1000-combination grid-search NMSE comparison
  (`visualize_gp.m`) is therefore a known limitation for any head-to-head
  claim at Milestone 12: it rests on training data that can't currently be
  verified to match what's on disk. This is reported as-is, not silently
  worked around — it is the user's call whether to regenerate that corpus
  with the current pipeline if a true 1000-combination comparison is wanted
  later.

## Known risks going in

1. **No persistent model object.** Unlike `fitrgp`, `griddatan` re-triangulates
   from scratch on every call. Naively called inside the parameter-recovery
   optimization loop (hundreds of repeated queries per row), this would be
   intractable in 6-D. **Resolved at Milestone 5**: the reusable predictor is
   a k-d tree (`createns`/`knnsearch`), proven exact (diff=0) against raw
   `griddatan` on both toy and real data. Measured impact: a realistic
   repeated-single-query workload (the shape the recovery loop actually
   uses) runs in 0.086 s where calling `griddatan` fresh each time would
   have cost roughly 4.6 hours (± ~25%; single-call cost varies 17.8–29.8 s
   across repeated measurements).
   **Correction to the original resolution note:** `griddatan` does use
   Delaunay/Qhull internally for *both* `'nearest'` and `'linear'` —
   `'nearest'` locates the nearest vertex of the enclosing simplex via
   `dsearchn`, not a global nearest-neighbor search, and pre-averages
   exact-duplicate training rows before triangulating. (Milestone 4's
   timing test showed no interior/exterior asymmetry for `'nearest'`
   because `dsearchn`'s simplex walk is cheap either way — not because no
   triangulation happens.) Consequence, confirmed by a forced-tie stress
   test: when a query is exactly equidistant from 2+ training points, the
   k-d tree and `griddatan` disagree roughly half the time (101/200 in a
   randomized 6-D test). This is **not observed** on any toy or real query
   tested (diff=0 throughout) — the real 3000-row training corpus has zero
   exact-duplicate rows and a minimum pairwise distance of 0.00994 — but it
   is a genuine, confirmed limitation of the k-d tree as a general
   substitute, not a hypothetical one. `fit_nearest_predictor` now warns if
   a future training set contains duplicate rows. Full detail in
   `STATUS.md`'s Milestone 5 entry.
2. **No native predictive variance.** The GP's delta-method confidence
   interval (used on the gamma-recovery plot) relied on `predict()`'s
   uncertainty output; `griddatan` has no equivalent. **Resolved at
   Milestone 9**: no substitute metric was fabricated; the illustrative
   recovery plots carry an explicit on-plot note that no
   predictive-uncertainty band is available, rather than a silently
   missing or fake one.
   **Bigger finding at the same milestone**: gamma recovery itself fails
   outright on this data (mean `log_error_gamma` = 1.53 across the full
   1000-row holdout, i.e. off by ~34x on average) — root-caused to raw,
   unstandardized feature distance: delta/epsilon (scale ~10-100) dominate
   the nearest-neighbor search over gamma1 (scale ~1e-5 to 5e-2) by
   multiple orders of magnitude (confirmed on 4 rows at Milestone 9's
   review: 2.5 to 10 orders of magnitude depending on the row), making the
   search functionally blind to gamma1 regardless of which row is checked.
   The same mechanism explains Milestones 6-8's NMSE>1 results. A
   diagnostic-only re-run with standardized features improves but does not
   fix this (mean drops to 1.17, still a failure) — full detail in
   `STATUS.md`'s Milestone 9 entry.
3. **Convex-hull boundary failure.** The existing draft (`interpolant.m`)
   already found `'linear'` produces `NaN` for query points outside the
   training convex hull, and dropped it in favor of `'nearest'`. This project
   re-verifies that independently rather than assuming it still holds.
4. **Prior-art red flag.** Two existing one-off `'nearest'`-method runs
   already on disk (`../gaussian_process-interpolant/griddatan_saved_models/`,
   `FM_nearest.mat` and the mangled `FM_nearest_✐.mat` — see the filename
   bug below) report flux/emc NMSE of **1.65–1.90** — worse than predicting
   the mean (NMSE = 1). Their source script (`interpolant.m`) sets
   `training_dataset_size = 10000` in code but its own comment says "5000
   points" when explaining the convex-hull NaN issue — an internal
   inconsistency, so the training size behind these specific NMSE numbers
   isn't fully trustworthy either. **Resolved at Milestone 6**: re-run
   directly (same seed, same `training_dataset_size=10000`, same
   `qhull_opts`) — CONFIRMED, not just carried over: flux NMSE 1.6516, emc
   NMSE 1.8578, closely matching the prior-art claim. The parity-anchored
   `dataset_size=3000` corpus (Milestones 6-7) is even worse: flux NMSE
   1.8309, ranking dead last against all 307 verified GP configs (249x
   worse than the best, ~5.9x worse than the GP's own worst config).
   `qhull_opts` itself confirmed to make zero difference to `'nearest'`
   predictions (max abs diff = 0).
5. **Bootstrap parity.** `gaussian_process.m` uses `n_boot = 1000`,
   `boot_size = 1000` for its bootstrap NMSE distributions. `griddatan`'s
   training/validation stage (Milestone 6) commits to the same two numbers —
   not left implicit.
6. **`readtable` delimiter auto-detection is unreliable on this project's
   GP-family result CSVs.** Found at Milestone 7: a GP results CSV's `tags`
   column (e.g. `FM_exact__BF_linear__KF_ARDMatern32__OPT_fminunc`) has far
   more underscores than the file has commas, which fools MATLAB's
   delimiter auto-detection into silently shredding columns AND dropping a
   row (confirmed: `gp_comparison_results_superset_3000.csv` reads as
   306 rows / 11 bogus columns without a delimiter, vs. the correct 307
   rows / 3 columns with one). Confirmed to also affect
   `gp_comparison_results_superset_5000.csv`,
   `gp_comparison_results_superset.csv`, and
   `gp_comparison_results.csv` (the 1000-combination corpus in risk 4's
   first bullet). **Standing rule: every `readtable` call on a GP-family
   results CSV in this project must pass `'Delimiter', ','` explicitly** —
   do not assume a file is safe to auto-detect just because another one
   parsed correctly (`comparison_table.csv` parses fine without it, purely
   by luck of having more columns).

## Planned file layout

Numbered to run in order; see `STATUS.md` for which of these actually exist
and have passed adversarial review versus which are still planned.

Prefixed `s##_` (not a bare leading digit) — MATLAB can't execute a script
whose filename starts with a digit, whether invoked via `run('00_x.m')` or
by typing the bare name directly; a leading digit isn't a valid identifier,
and MATLAB parses a bare invocation as one (confirmed at Milestone 2: a
trivial `00_test.m`, content `disp(1+1)`, throws a misleading "Invalid text
character" error via both call paths; `s00_test.m` with identical content
runs fine either way). The original `00_`/`01_`/... scheme in an earlier
draft of this table was wrong for that reason.

| File | Milestone | Purpose |
|---|---|---|
| `s00_data_and_scaling_parity.m` | 2, 3 | Load identical train/val subsets + scaling constants; verify parity with GP side |
| `s01_toy_validation.m` | 4 | Sanity-check `griddatan` on a synthetic function with a known answer |
| `s02_fast_predictor.m` | 5 | Reusable `delaunayn`+`tsearchn` predictor, parity-tested against raw `griddatan` |
| `s03_train_validate_models.m` | 6 | Full train+validate sweep, NMSE + bootstrap, saved per-run |
| `s04_build_comparison_table.m` | 7 | Comparison CSV, joinable with the GP side's |
| `s05_compare_nmse_boxplot.m` | 8 | Combined GP + griddatan NMSE visualization |
| `s06_parameter_recovery_gamma.m` | 9 | Gamma recovery using the griddatan surrogate |
| `s07_parameter_recovery_joint.m` | 10 | Joint (gamma, lambda) recovery |
| `run_all.m` | 11 | Single entry point, clean-run reproducibility check |
| `s08_final_report.m` | 12 | Head-to-head GP vs. griddatan verdict |

## Status

See `STATUS.md` for the live milestone log, including any failures.

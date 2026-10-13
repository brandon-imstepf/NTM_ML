# Milestone Status Log

Live record of what's been attempted, what happened, and the adversarial
review verdict for each milestone. A failure is reported here, not hidden —
downstream milestones that depend on a failed one are marked **blocked**,
not silently skipped.

Verdict legend: **PASS** / **PASS WITH CAVEATS** / **BLOCKED-FAILED** / **IN PROGRESS**

---

## Milestone 12 — Final head-to-head report

**Status:** PASS WITH CAVEATS (final). The caveat is the same one carried
through the whole project: griddatan loses decisively, reported as the
honest result, not a shortcoming of the report.

**What was done:** `s08_final_report.m` synthesizes Milestones 6-10's
already-saved results into one head-to-head comparison, with one new
computation: the GP side's **full-dataset** gamma-recovery accuracy was
never persisted to disk by `synthetic_parameter_recovery_gp.m` (read-only
reference) — that script computes `log_error_gamma` in memory but only
ever uses it to pick its best-10 rows to plot, never saves the aggregate.
This milestone reproduces that same existing computation (identical
formula, identical full 1000-row holdout, same GP model
`FM_exact__BF_linear__KF_ARDMatern32__OPT_fminunc`, read-only against its
saved `.mat`) so the comparison is genuinely apples-to-apples on full
datasets on both sides, not griddatan-full-dataset-vs-GP-best-10. Result
saved to `results/gp_gamma_recovery_full_reference.csv`. A transcript of
the full report is saved to `results/final_report.txt` via `diary()`.

**Headline numbers, all pulled from already-adversarially-reviewed
milestones except the newly-added GP gamma-recovery figure (not yet
independently reviewed — see below):**

| Axis | GP | griddatan |
|---|---|---|
| Flux NMSE (dataset_size=3000) | best 0.00736, worst 0.31015, median 0.29619 (307 configs) | 1.8309 (1 config) |
| Rank among all 308 configs | 1-307 | **308 (dead last)** |
| Gamma recovery, full 1000-row holdout | mean log_error_gamma = **0.2942** (median 0.0687) | mean = **1.5287** (median 1.5091) |
| Joint recovery, full 1000-row holdout | not computed (see caveat) | mean log_error_joint = 1.5721; 997/1000 rows collapse to one fixed point |
| "Training" wall time (one model) | 295.64 s (`fitrgp`, this model) | ~0.02-0.04 s (k-d tree fit) |

Even GP's single *worst* config (out of 307) beats griddatan's NMSE by
5.9x; GP's best beats it by 249x. On gamma recovery, GP is off by a mean
of 0.29 decades (~2.0x) vs. griddatan's 1.53 decades (~34x) on the
identical holdout set.

**Caveats stated plainly, not hidden:**
- The GP side's OWN 1000-combination grid-search NMSE comparison
  (`visualize_gp.m`, the unsuffixed `saved_models/`) is NOT used anywhere
  in this report — its training-data provenance could not be verified
  (Milestone 1). Every GP number here comes only from the verified
  `saved_models_superset_3000/` corpus (307 configs).
- GP-side joint (gamma,lambda) recovery was **not computed** for this
  report (scope/time: a 1000-row x 2500-point grid-search+polish via
  `predict()` is substantially larger than the gamma-only pass). The
  report explicitly says GP joint recovery is *expected* to also
  outperform griddatan given the other results, but states this is not an
  independently verified number and should not be quoted as one.
- griddatan's raw (non-fast-predictor) per-query cost was never
  benchmarked head-to-head against GP's `predict()` at the single-query
  level; the report notes this rather than implying a comparison that
  wasn't actually measured.

**Root cause recap, carried into the final verdict rather than just
asserted:** griddatan's raw, unstandardized feature distance is dominated
by delta/epsilon (scale ~10-100) over gamma1/lambda1 (scale ~1e-5 to
1e-1) — independently confirmed at Milestones 9 and 10, not a tuning
problem but a structural one inherent to how griddatan (as this project's
own `interpolant.m` originally built it) computes nearest-neighbor
distance.

**Adversarial review:** complete (independent agent). Verdict: **PASS** —
"the strongest-verified milestone of the series" in the reviewer's words.
The one genuinely new computation (GP full-dataset gamma recovery) was
reproduced via a fully independent from-scratch script: a full 1000-row
rerun matched mean=0.2942/median=0.0687 to 4 decimals, and a row-level
cross-check against the saved CSV found max abs diff = 5.03e-17 (pure
floating-point noise). The GP model file and its `wall_time=295.6444s`
field were confirmed genuine (not stale). Every table number (NMSE
best/worst/median, rank 308/308, griddatan's gamma/joint figures) was
independently recomputed and matched exactly, including a second,
independent method for the rank calculation. All four final-verdict
multipliers (5.9x, 249x, ~34x, ~2.0x) were recomputed by hand and
confirmed correct. Grepped the script directly to confirm it never reads
the unverified `gp_comparison_results.csv` corpus despite referencing it
only in a caveat — the caveat is enforced in code, not just claimed in
a comment. `results/final_report.txt` confirmed to contain the complete
transcript. No scope violations, no overclaiming, no softened or spun
language found. One unrelated observation (not a violation): `git status`
showed concurrent, unrelated changes to `network_validation/` files
(an ODE-tolerance env-var override) made the same day — confirmed no
griddatan script references that folder; flagged for the user's own
awareness, outside this project's scope to act on.

---

## Milestone 11 — End-to-end reproducibility (`run_all.m`)

**Status:** PASS WITH CAVEATS (final)

**What was done:** `run_all.m` runs all 8 scripts built so far
(`s00` through `s07`) in dependency order from a genuinely clean state.
Each script is spawned as its own fresh `matlab -batch` subprocess rather
than chained as plain script calls in one shared workspace — every
`s##_*.m` script starts with `clc; clear;`, and a script-calling-a-script
shares one workspace in MATLAB (unlike a function call), so chaining them
directly would have each script's own `clear` wipe out `run_all.m`'s own
loop state. Spawning a subprocess per script sidesteps that and is a
stronger clean-session test besides: it proves each script works from
nothing, not just that it survives running right after another script
already populated the workspace.

**Result: full pipeline reproduces cleanly, all 8/8 scripts passed:**

| Script | Exit code | Wall time |
|---|---|---|
| `s00_data_and_scaling_parity` | 0 | 10.9-13.5 s |
| `s01_toy_validation` | 0 | 205.4-252.1 s |
| `s02_fast_predictor` | 0 | 113.1-114.0 s |
| `s03_train_validate_models` | 0 | 714.7-742.6 s |
| `s04_build_comparison_table` | 0 | 9.8-12.3 s |
| `s05_compare_nmse_boxplot` | 0 | 39.0-41.9 s |
| `s06_parameter_recovery_gamma` | 0 | 73.9-74.5 s |
| `s07_parameter_recovery_joint` | 0 | 72.3-85.7 s |
| **Total** | — | **1258-1318 s (~21-22 min), across 3 independent full runs** |

`s03`'s 742.6 s is expected, not a regression: its Part B calls raw
`griddatan` (not the fast predictor) four times on the 10,000-row corpus,
each individually slow (30-200s, consistent with every other raw
`griddatan` timing measured since Milestone 4).

**Every number re-verified identical to its original milestone run** —
spot-checked against this log: Milestone 2/3's parity diffs (all 0),
Milestone 5's exact-parity diffs (all 0), Milestone 6's NMSE (1.8309 /
1.8158 / 1.6516 / 1.8578), Milestone 7's ranking (309/309, 249x, 5.9x),
Milestone 9's `log_error_gamma` (1.5287 / 1.1717 diagnostic), Milestone
10's corner-collapse (997/1000) — all bit-identical to the values already
committed to this document. Nothing drifted on a from-scratch re-run.

Saved to `results/run_all_summary.csv` and per-script logs in
`results/run_all_logs/`.

**Adversarial review:** complete (independent agent). Verdict: **PASS
WITH CAVEATS**. The 8/8 pass, bit-identical numbers, clean subprocess
isolation (watched via `tasklist` across the whole run), and path
independence (spot-checked `s00`/`s03`/`s04`/`s06`/`s07`) were all
independently reconfirmed via the reviewer's own full `matlab -batch
"run_all"` re-run (1178.9s, same shape, same bottleneck script). The
"`clear` wipes the caller's workspace" design rationale was empirically
demonstrated, not just taken on faith, with a minimal parent/child test.

**Real bug found and fixed:** the original log-capture mechanism
(`system('matlab -batch "name" > log.txt 2>&1')`, shell-level redirect
around a nested `matlab -batch` call) silently never wrote anything —
confirmed by the reviewer three independent ways (a standalone
good/bad-child synthetic test, watching `results/run_all_logs/` stay
empty through an entire real run, and isolating the failure to nested
`matlab -batch` specifically, not `system()` redirection in general).
Exit-code capture was unaffected; only the log artifact was silently
missing, while STATUS.md had claimed it as verification evidence
("spot-checked against this log") — a real, if narrow, overclaim.
**Fixed** by switching to the two-output form of `system()`
(`[status, cmdout] = system(cmd)`), which captures the subprocess's
combined stdout/stderr directly into a MATLAB variable regardless of
shell-level redirection, then writing `cmdout` to the log file explicitly
via `fopen`/`fprintf`/`fclose`. Verified working: a targeted sanity check
confirmed a populated log file before committing to a full re-run, then a
full `run_all.m` re-run (results above) produced all 8 log files in
`results/run_all_logs/` with correct, readable content (spot-checked
`s06`'s log directly). Also addressed two minor code-quality items from
the same review: wrapped the `system()`/log-write block in `try`/`catch`
so `cd(old_dir)` always runs even if something throws (re-thrown, not
swallowed), and added an explicit comment acknowledging `run_all.m` has
no subprocess-timeout/hang protection (a real but accepted limitation
given every script currently completes in well under 15 minutes). No
scope violations by the reviewer; unrelated pre-existing changes in
`network_validation/` were present before the review started and are not
attributable to this milestone.

---

## Milestone 10 — Joint (gamma, lambda) recovery port

**Status:** PASS WITH CAVEATS (final). Same framing as Milestone 9: the
engineering is sound, the recovery RESULT is a clear failure, reported
plainly.

**What was done:** `s07_parameter_recovery_joint.m` ports
`synthetic_parameter_recovery_gp.m`'s `jointRecoveryFit`/`jointRecoveryPlot`
(read-only reference): a 50x50 log-spaced (gamma1, lambda1) grid, coarse
grid-search for a starting point, `fminsearch` polish in log-space. Run
over the **full** 1000-row holdout set — a deliberate deviation from the
GP script, which only covers its first 100 rows (`size(dataset,1)/10`) as
a performance shortcut; griddatan's fast predictor (Milestone 5) makes
that shortcut unnecessary here.

**Predicted outcome, stated before running, then checked rather than
assumed:** Milestone 9 found griddatan's nearest-neighbor search
functionally blind to gamma1 because delta/epsilon's raw scale dominates
distance. Since lambda1's raw scale (~0.01-0.1) is just as dwarfed, the
same mechanism predicts the 2-D (gamma1, lambda1) loss surface will also
be a flat plateau. Confirmed true.

| Statistic | Value |
|---|---|
| mean `log_error_joint` (1000 rows) | **1.5721** |
| median `log_error_joint` | **1.5443** |
| std | 0.7919 |
| Recovery wall time | 6.2-7.3 s for all 1000 rows (~0.007 s/row) |
| Rows collapsed to the grid's starting corner | **997/1000 (99.7%)** |
| Rows landing on the physical range boundary (gamma or lambda) | 0/1000 (a different, looser check — see below) |

**Confirmed directly, not just inferred:** all 10 illustrative plots (the
best-10 rows by `log_error_joint`) show a **completely flat** 50x50 loss
heatmap — a single solid color across the entire grid, log10(squared
error) varying by less than the plot's own color resolution. `contour()`
had nothing to draw on any of them (see bug below).

**"Best 10" selection is not genuine recovery — independently confirmed,
not just deduced:** on a genuinely flat loss surface, `min()`
deterministically returns the *first* tied index in column-major order —
i.e. `(gamma_grid(1), lambda_grid(1))`, the grid's own
smallest-gamma/smallest-lambda corner — regardless of the actual data.
Confirmed directly (now also checked in-script, see the added
`collapsed_to_corner` diagnostic): **997 of all 1000 holdout rows recover
a `(gamma1, lambda1)` that is bit-identical to that one fixed corner** —
only 4 unique `recovered_gamma1` values exist across the entire 1000-row
output. The "best 10 by error" rows are therefore not genuine successful
recoveries — they are simply the holdout rows whose TRUE (gamma1,
lambda1) happen by coincidence to be closest (in log-space) to that one
fixed point (e.g. row 516: true=(5.999e-5, 1.018e-2) vs. corner
(5.764e-5, 1.013e-2) → small error; row 822: true further away → larger
error). Only 3 of 1000 rows (568, 657, 699, all with unusually small true
gamma1) landed anywhere else, 2-2.8 decades off the corner.

**Bug found and fixed:** the first run's `std_log_loss > 0` flatness check
was too strict — floating-point jitter left `std_log_loss` tiny but
nonzero, so the code attempted to draw contour lines anyway. MATLAB's
`contour()` has its own near-equality tolerance for "constant data" and
silently skipped drawing anyway, firing a warning on every one of the 30
contour calls (10 plots x 3 contour levels) AND leaving the legend with
two entries ("mean log-loss", "mean ± 1σ log-loss") pointing at lines that
were never drawn. Fixed with a relative tolerance
(`std_log_loss < 1e-6 * max(|mean_log_loss|,1)`) and, when flat, skipping
the contour calls entirely (not just letting them fail silently) and
showing an on-plot text note instead — confirmed clean re-render, zero
warnings, no dangling legend entries.

**Full per-row results** saved to `results/joint_recovery_full.csv` (1000
rows), 10 illustrative plots saved to `figures/joint_recovery/`.

**Adversarial review:** complete (independent agent). Verdict: **PASS
WITH CAVEATS**. Every claim independently reproduced via a fully separate
code path, including the single most important one (corner-collapse):
reconstructed `gamma_grid(1)`/`lambda_grid(1)` from scratch and checked
all 1000 CSV rows directly (not just the best-10), confirming 997/1000
bit-identical to the corner — the review found this story *stronger* than
originally written, not weaker, and recommended upgrading the language
from "deduction, pending verification" to a confirmed fact, which this
entry now does. Headline stats matched to 1e-14 (float noise only) via an
independent sample and a from-scratch CSV recomputation (cross-checked in
Python, not just MATLAB). The flat-surface finding was confirmed on 8
additional rows never in the original best-10 (std of log10-loss
2.8e-17-8.9e-16 — genuinely flat, not cherry-picked). The contour-warning
fix was confirmed clean (zero warnings on a full re-run). **One fix
applied from review:** the "0/1000 at range boundary" check tests against
`gamma1_range`/`lambda_range` (the physical parameter bounds), a
different and much looser reference than `gamma_grid(1)`/`gamma_grid(end)`
(the actual corner recovery collapses to) — the original in-code comment
did not make this distinction clearly enough and could have misled a
future reader into thinking that check captured the corner-collapse
story. Added a proper `collapsed_to_corner` diagnostic (now 997/1000 is
computed directly in-script, not only in this document) and reworded both
comments to distinguish the two checks explicitly. No scope violations.

---

## Milestone 9 — Gamma-recovery port

**Status:** PASS WITH CAVEATS (final). "Pass" refers to the engineering
(the port runs correctly and the feasibility metric is computed honestly);
the actual recovery RESULT is a clear failure, reported as such per the
user's explicit instruction to report griddatan failing at any stage, not
hide it.

**What was done:** `s06_parameter_recovery_gamma.m` ports the gamma-recovery
procedure from `../parameter_recovery/synthetic_parameter_recovery_gp.m`
(read-only reference): hold (lambda1, delta, epsilon, N1, N2) at their true
per-row values, `fminbnd` over gamma1 to match observed standardized flux,
over the **full** 1000-row holdout set (not a cherry-picked subset) —
`log_error_gamma` computed for every row, matching how the GP side already
uses this statistic. Ported to the k-d tree fast predictor instead of
`predict()`; same `gamma1_range=[5e-5,5e-2]`. Only one griddatan "model"
exists, unlike the GP script's 4 near-tied hyperparameter configs.

**Headline result — griddatan parameter recovery is infeasible on this
data:**

| Statistic | Value |
|---|---|
| mean `log_error_gamma` (1000 rows) | **1.5287** |
| median `log_error_gamma` | **1.5091** |
| std | 0.8764 |
| Recovery wall time | 1.84 s for all 1000 rows (0.0018 s/row) |

`log_error_gamma=1` means the recovered gamma1 is off by a full decade
from the true value. A mean of 1.53 means griddatan's recovered gamma1 is,
on average, off by **more than 1.5 decades (~34x)** from the truth — not a
noisy-but-usable estimate, a non-recovery.

**Root cause, independently confirmed, not just inferred from the bad
number:** the illustrative plots (below) showed a completely FLAT
loss-vs-candidate-gamma curve for every inspected row. Checked directly in
the script, not just eyeballed: for a representative row, the
nearest-neighbor index returned by the k-d tree is **identical** across
`gamma1 ∈ {5e-5, 0.001, 0.01, 0.05}` — a 1000x span of the entire allowed
range. Squared per-dimension distance contributions for that match:
gamma1 = 3.3e-6, lambda1 = 5.1e-9, N1 = 5.1e-8, N2 = 1.2e-8, **delta =
0.333, epsilon = 1.49**. Delta and epsilon (raw scale ~10-100) dominate
the raw Euclidean distance by 5-6 orders of magnitude over
gamma1/lambda1/N1/N2 (raw scale ~1e-5 to ~1e-1). The nearest training
point is selected almost entirely by matching delta/epsilon; gamma1 is
functionally invisible to the search. **This is the same mechanism behind
the NMSE>1 results already reported in Milestones 6-8** — not a new,
separate problem, the same one showing up a second way.

**Diagnostic only (explicitly NOT part of the committed apples-to-apples
pipeline — interpolant.m's original design never standardized features,
only targets, and this project has followed that choice for parity since
Milestone 4):** re-ran the identical recovery loop with features z-scored
before building the k-d tree. Result: mean `log_error_gamma` drops from
1.5287 to **1.1717** — standardizing substantially helps, confirming the
root-cause mechanism — but recovery is **still a failure** even then (mean
still > 1 decade off). Scale disparity is a real, significant contributing
factor, not the sole cause; something about griddatan's nearest-neighbor
method, even on equal footing feature-wise, still cannot support gamma
recovery to a useful precision on this data. Reported honestly rather than
either oversold ("fixed by standardizing") or glossed over.

**Risk 2 (no native predictive variance) resolved by omission, stated
plainly:** the illustrative recovery plots carry an explicit on-plot note
("No predictive-uncertainty band: griddatan has no native predictive
variance") rather than a fabricated confidence band or a silently missing
one.

**Full per-row results** saved to `results/gamma_recovery_full.csv` (all
1000 rows: true gamma1, recovered gamma1, min loss, log_error_gamma) for
Milestone 12's final report to use directly, and 10 illustrative plots
(best-10 rows by `log_error_gamma`) saved to `figures/gamma_recovery/`.

**Adversarial review:** complete (independent agent). Verdict: **PASS
WITH CAVEATS**. Every number independently reproduced via a fully separate
code path (its own fminbnd/knnsearch loop, not calling this script's
functions): headline mean/median/std matched to 4 decimals; the root-cause
distance breakdown for row 70 matched exactly, and held qualitatively
(nearest-neighbor index pinned regardless of gamma1) on 3 additional rows
tested, though the specific "5-6 orders of magnitude" figure is row-70-
specific — it ranged 2.5 to 10 orders of magnitude across the rows tested,
still always large enough to make gamma1 functionally invisible, and the
STATUS.md text was already correctly scoped to "that match" rather than
claimed as universal; the standardized-feature diagnostic (1.5287→1.1717)
matched exactly via independent re-fit; the saved CSV cross-checked
bulletproof against both the console output and a fresh independent re-run
(diffs at floating-point noise level only). **Two fixes applied from this
review:** (1) the `xline` "True Gamma" label was clipped at the plot's top
edge in most of the 10 illustrative plots (found by viewing 5 of them
directly) — fixed via `'LabelVerticalAlignment','bottom'`, confirmed
legible on re-render; (2) this entry now notes the row-to-row variability
in the "orders of magnitude" figure rather than implying it's a fixed
constant. No scope violations, no overclaiming or softening found in the
writeup (language was independently judged "appropriately blunt").

---

## Milestone 8 — Cross-method NMSE visualization

**Status:** PASS (final)

**What was done:** `s05_compare_nmse_boxplot.m` extends the GP side's own
faceted-boxplot approach (`../gaussian_process-interpolant/visualize_gp.m`,
read-only reference) — same four facets (FitMethod, BasisFunction,
KernelFunction, Optimizer), same log-scale axis — to show griddatan
alongside it. Reads Milestone 7's `results/comparison_table.csv` rather
than re-parsing the GP CSV a second time. Filters to `dataset_size==3000`
only (excludes griddatan's dataset_size=10000 prior-art row — a different
corpus, not meant to be plotted next to this one). Produces two figures,
`figures/flux_nmse_boxplots_vs_griddatan.png` and the `emc_` equivalent.

**Design choice, stated rather than silently decided:** griddatan's single
result is drawn as a horizontal reference line (`yline`), not a fabricated
"box." It is one model, not a grid of 307 — giving it artificial spread
(e.g. by plotting its own `n_boot=1000` bootstrap distribution as a box)
would conflate "uncertainty from resampling one fixed model's residuals"
with "variability across many different model configurations," which
isn't an honest apples-to-apples comparison.

**Two bugs found and fixed while getting this to run:**
- `regexp` on a MATLAB `string` scalar (vs. `char`) returns the `'names'`
  struct wrapped in a 1×1 cell, not the struct itself — `tok.fm` then
  fails with "Dot indexing is not supported for variables of type cell."
  Fixed by wrapping with `char(...)` before the regex call, matching
  `visualize_gp.m`'s original (char-based) behavior exactly.
- The `yline` label ("griddatan (nearest)") was invisible on first render
  — the axis auto-scaled to exactly the line's height, clipping the label
  at the plot's top edge. Fixed by extending `ylim` to give ~60% headroom
  above the line; confirmed visually on the re-rendered figure that the
  label is now readable on all 8 tiles (4 facets × 2 targets).

**Result:** visually unambiguous — the griddatan reference line sits above
every single GP box in both figures, across all four facets. Consistent
with Milestone 7's numeric finding (GP's worst config still beats
griddatan by ~5.9x) — nothing here contradicts that, it's the same result
shown a different way.

**Adversarial review:** complete (independent agent). Verdict: **PASS**, no
bugs found. Independently re-verified via a second, non-regexp parsing
method (307/307 field matches), a fresh re-run of the script, and direct
visual inspection of both output PNGs (confirmed: griddatan's line clears
every GP box with zero ties — max GP flux NMSE 0.310148 vs. griddatan's
1.830932; max GP emc NMSE 0.344771 vs. 1.815821 — and the label is legible
on all 8 tiles, not clipped). One trivial, non-blocking style nit noted
(`ylim(gca)` read twice on one line) — not worth a follow-up edit. No
scope violations.

---

## Milestone 7 — Comparison table

**Status:** PASS WITH CAVEATS (final)

**What was done:** `s04_build_comparison_table.m` reads the GP side's own
`gp_comparison_results_superset_3000.csv` (read-only) and this project's
Milestone 6 output, and writes a single joined CSV
(`results/comparison_table.csv`) with columns `tags, flux_nmse, emc_nmse,
family, dataset_size` — `family` (`"GP"` / `"griddatan"`) and
`dataset_size` added explicitly so a future reader never has to guess
which rows are on the same corpus.

**Bug found and fixed before this could run at all:** `readtable` on the
GP CSV silently mis-detected the delimiter as `_` instead of `,` — the
`tags` column's values (e.g.
`FM_exact__BF_linear__KF_ARDMatern32__OPT_fminunc`) contain far more
underscores than the file has commas, which fooled MATLAB's auto-detection
into shredding every tag into a dozen bogus columns (`FM`, `exact`, `Var3`,
`BF`, ... one of which silently absorbed the actual `flux_nmse` value into
a column named `fminunc_0_007362`). First run failed outright
(`Unrecognized table variable name 'tags'`) rather than silently
corrupting data, which is how this was caught. Fixed by passing
`'Delimiter', ','` explicitly to every `readtable` call reading one of
these files — confirmed via a standalone test before and after the fix.

**Adversarial review widened this finding, independently reproduced here
before being written down:**
- The no-delimiter read doesn't just shred columns — it also silently
  **drops a row** (306 of the file's true 307 data rows come back).
  Confirmed directly: `readtable(...)` with no delimiter on
  `gp_comparison_results_superset_3000.csv` gives `height=306`; with
  `'Delimiter',','` gives the correct `height=307`.
- **Three more files in `gaussian_process-interpolant/` have the identical
  vulnerability and are not yet guarded anywhere:**
  `gp_comparison_results_superset_5000.csv`,
  `gp_comparison_results_superset.csv`, and — most importantly —
  `gp_comparison_results.csv` (the 999-row, 1000-combination grid-search
  corpus README.md already flags as unverified provenance). All three
  confirmed to shred to `width=11`, no `tags` column, when read without an
  explicit delimiter. `griddatan_comparison_results.csv` (short tags) is
  unaffected. **Standing rule for every future milestone that reads a
  GP-family results CSV (Milestone 12 will likely read
  `gp_comparison_results.csv`): always pass `'Delimiter', ','`
  explicitly — do not assume a file is safe just because another one
  parsed correctly** (`comparison_table.csv` itself parses fine without
  it, incidentally, because its extra `family`/`dataset_size` columns
  change MATLAB's auto-detection — that's luck, not a guarantee).

**Result:** 309 rows (307 GP configs + 2 griddatan runs:
dataset_size=3000 and 10000). Griddatan's dataset_size=3000 row ranks
**309th of 309** by flux NMSE — dead last, not close (its own
dataset_size=10000 row is 308th, so griddatan holds both of the bottom two
slots). Best GP config: flux NMSE 0.00736. Griddatan (3000): 1.8309.
**249x worse than the best GP model.** GP's *worst*-performing config in
its own grid search (flux NMSE = 0.310148,
`FM_none__BF_constant__KF_ARDExponential__OPT_fmincon`) still beats
griddatan outright, by **~5.9x** — independently confirmed (previously
flagged here as "visibly true but not confirmed"; now resolved).

**Fix applied from review:** the script's summary line originally read the
best GP NMSE via positional indexing (`T_gp.flux_nmse(1)`), correct only
because the source CSV happens to already be sorted ascending on disk — a
silent-misreport risk if that ever changes. Replaced with
`min(T_gp.flux_nmse)` / `max(T_gp.flux_nmse)`; re-run confirmed identical
numeric output.

**Adversarial review:** complete (independent agent). Verdict: PASS WITH
CAVEATS. Confirmed via its own independent MATLAB execution: the join is
exact (no row misalignment, no column swap), the rank-309 claim, the
249x and ~5.9x figures, no scope violations (`git status` at `NTM_ML` root
showed only this project's own untracked files plus two unrelated
pre-existing changes elsewhere, neither touched by this work). All
consequential claims independently re-verified here as well before being
written into this entry.

---

## Milestone 6 — Full train+validate pipeline

**Status:** PASS WITH CAVEATS (final)

**What was done:** `s03_train_validate_models.m` runs two independent parts.

**Part A (primary, parity-anchored — the Milestone 7+ comparison target):**
`dataset_size = 3000`, the same corpus Milestone 2 proved matches the GP
pipeline exactly. Flux and emc share one k-d tree (same `x` for both), with
`n_boot = 1000`, `boot_size = 1000` matching `gaussian_process.m` exactly.

| Metric | Flux | Emc |
|---|---|---|
| Validation NMSE | **1.8309** | **1.8158** |
| Training (in-sample) NMSE | 0.000000 | 0.000000 |

- **Validation NMSE > 1 for both targets** — `griddatan('nearest')` at
  `dataset_size=3000` performs *worse than predicting the mean* on this
  project's actual parity-anchored corpus. Reported as-is; this is the
  real number Milestone 12's head-to-head report will use, not a
  cherry-picked or best-case run.
- **Training NMSE = exactly 0**, asserted in-script (not just observed):
  `'nearest'` is an exact interpolator — every training point is its own
  nearest neighbor, since Milestone 5 close-out already confirmed this
  corpus has zero duplicate/near-duplicate rows. The GP, by contrast,
  smooths and generally has nonzero training residual. A real, structural
  difference between the two surrogate families, not a bug.
- A standardize/un-standardize round-trip no-op check was asserted
  (tolerance `1e-9`, not bit-exact `0` — floating-point division then
  multiplication by the same constant doesn't round-trip bit-exact; this
  tripped the script's first run, fixed immediately before re-running).

**Part B (prior-art re-check — resolves README.md's Risk 4):**
Reproduces `interpolant.m`'s own original setup directly:
`training_dataset_size = 10000`, same seed, same `qhull_opts =
{'Qt','Qbb','Qc','Qx'}`.

| Metric | Flux | Emc |
|---|---|---|
| NMSE | **1.6516** | **1.8578** |

- **CONFIRMED, not just assumed:** this closely reproduces the prior-art
  claim ("NMSE 1.65–1.90, worse than the mean") under controlled,
  documented conditions — Risk 4 is resolved, not carried over untested.
- **`qhull_opts` vs. no-opts: max abs diff = 0.** The custom Qhull options
  in `interpolant.m` make no difference to `'nearest'` predictions on this
  data. This also extends Milestone 5's exact-parity proof (fast predictor
  vs. raw `griddatan`) to the `qhull_opts` configuration, since both give
  identical predictions to the plain call already proven exact.

**Known risk this surfaces for Milestone 12:** griddatan underperforming
the mean at both tested dataset sizes (3000 and 10000) is a strong, two-way
confirmed result, not a fluke of one run — it should be stated plainly in
the final head-to-head report rather than softened.

**Adversarial review:** complete (independent agent, took ~19 minutes this
time — no repeat of Milestone 5's stall). Verdict: **PASS**. Every
numeric claim reproduced via a genuinely independent method, not a re-run
of the same code path: a from-scratch brute-force pairwise-distance search
(no `knnsearch`) matched Part A's NMSE to 6 decimal places, and a direct
raw `griddatan('nearest',...)` call (the real Delaunay/`dsearchn` path,
not the k-d tree) matched the brute-force result with max abs diff = 0.
Training-NMSE=0 was stress-tested with a full 3000×3000 pairwise distance
matrix (confirmed: every row's true nearest neighbor is itself, min
off-diagonal distance 0.009945, no near-duplicates found). The 1e-9
tolerance fix was confirmed legitimate (actual diff ≈1.654e-24, pure
floating-point noise). Part B's qhull_opts=0-diff claim was checked at 8
individual validation points, not just in aggregate — all 8 exact. One
non-blocking item: `git status` showed two unrelated-looking changes in
`network_validation/` timestamped the same minute as this milestone's
work; reviewed and attributed to coincidental concurrent activity, not
this script (no code path in `s03` touches that folder) — flagged for the
user's awareness, not a scope violation.

---

## Milestone 5 — Fast reusable predictor

**Status:** PASS WITH CAVEATS (final)

**What was done:** `s02_fast_predictor.m` builds a "fit once, query many
times" nearest-neighbor predictor (`createns`/`knnsearch`, a k-d tree) to
replace calling `griddatan('nearest',...)` fresh on every query. Design
choice follows directly from Milestone 4: since `'nearest'` doesn't use a
Delaunay/Qhull triangulation internally (no interior/exterior timing
asymmetry like `'linear'` showed), a k-d tree is the right reusable
structure — not the `delaunayn`+`tsearchn` wrapper originally planned in
`README.md`'s file-layout table, which was written before Milestone 4's
finding. `README.md` updated to match.

**Parity — exact, not approximate, in every test run:**
- Toy problem (interior queries): `griddatan` vs. fast predictor, max abs
  diff = **0**.
- Toy problem (exterior queries): max abs diff = **0**.
- Real data: the Milestone 2-verified 3000-row training corpus, full
  1000-row validation batch, flux target: max abs diff = **0**.

**Performance — measured, not estimated, for the whole-batch case:**
- One-time fit cost (3000 training rows): 0.023 s.
- `griddatan('nearest',...)` for the whole 1000-row validation batch:
  33.1 s.
- Fast predictor for the identical batch (after the one-time fit): 0.006 s.

**Performance for the workload that actually matters — repeated
single-point queries, the shape fminbnd and the joint-recovery grid sweep
actually use:** 500 repeated single-point queries via the fast predictor:
0.086 s total (0.0002 s/query). Calling `griddatan` fresh for each of those
500 queries would cost roughly 500 × 33.1 s ≈ 16,561 s (~4.6 hours) by
direct extrapolation from the measured single-call cost above. This is the
concrete number behind the Milestone-1 risk ("likely intractable in 6-D") —
now resolved, not just mitigated.

**Known limitation — now root-caused and confirmed, not just flagged:**
The adversarial reviewer (agent `a90d7cea58aa60db1`) pushed past "ties are
undocumented, so test for them" into actually reading `griddatan.m`'s
source and constructing a forced-tie stress test. Both findings were
independently reproduced here before being written down:

- **Mechanism** (confirmed by reading MATLAB's own `griddatan.m`, Milestone
  5 close-out): `griddatan` is **not** a plain nearest-neighbor search for
  either `'linear'` or `'nearest'`. Both methods always call `delaunayn`
  first; `'nearest'` then calls `dsearchn(x, tri, xi)`, which finds the
  nearest *vertex of the enclosing simplex*, not a global nearest neighbor.
  This also means `griddatan` pre-deduplicates exact-duplicate training
  rows by **averaging** their target values (`unique(x,'rows')` +
  `accumarray(...,@mean)`) before triangulating — the k-d tree predictor
  does no such thing; it treats every row as distinct.
- **Confirmed mismatch, reproduced independently:** a minimal 2-D forced
  exact-tie case (four points exactly equidistant from the query) gives
  `griddatan = 40`, k-d tree = `10` — a real disagreement, not a rounding
  artifact. The reviewer's own 6-D randomized forced-tie stress test found
  **101/200 (50.5%) mismatches**, 0 errors — i.e., when an exact tie is
  forced, the two methods have no reliable agreement and disagree roughly
  half the time. This is the expected consequence of simplex-vertex search
  (`dsearchn`) vs. global distance search (`knnsearch`), not a bug in
  either implementation.
- **Does this affect the data this project actually uses?** No, as far as
  tested: the real 3000-row training corpus has **zero exact-duplicate
  rows** and a minimum pairwise distance of 0.00994 (not near-zero) —
  consistent with the diff=0 parity results already measured on real data.
  This is a confirmed structural limitation of the fast predictor as a
  general substitute for `griddatan`, not a problem observed in this
  project's actual corpus. `fit_nearest_predictor` now warns if a future
  training set ever does contain duplicate rows, so this isn't a silent
  risk if the corpus changes later.
- **Feature-scale disparity — checked, not an issue:** the reviewer also
  tested whether the real data's large per-dimension scale spread
  (gamma1/lambda1 ~1e-3–1e-1 vs. delta/epsilon ~10–100) could make
  `createns`'s raw-Euclidean k-d tree disagree with `griddatan`'s internal
  distance metric. Confirmed diff=0 — both sides use unscaled Euclidean
  distance identically, so this is not a source of mismatch.
- **Timing stability — checked:** five repeated single-point
  `griddatan('nearest',...)` calls on the real data ranged 17.8–29.8 s
  (mean 22.6 s, std 5.3 s). The "~4.6 hours" extrapolation in the
  Performance section above is therefore a reasonable order-of-magnitude
  estimate from a single representative call, not a precise constant —
  call-to-call cost varies by roughly ±25%.

**Adversarial review:** complete (agent `a90d7cea58aa60db1`). The agent's
background MATLAB verification run finished but the agent itself did not
return a final written verdict after ~15 hours (last activity: idle,
waiting on its own background task notification). Per the user's
direction, its raw test output was read directly from its task transcript
and every consequential claim (the mechanism, the forced-tie mismatch rate,
the scale-disparity check, the timing stability) was independently
reproduced here before being committed to this document — the same
discipline used for Milestones 1, 2, and 4's reviews.

---

## Milestone 4 — Toy validation

**Status:** PASS WITH CAVEATS

**What was done:** `s01_toy_validation.m` tests `griddatan` on a known,
smooth, 6-D function (2000 training points, 500 test points — same order of
magnitude as the real problem), both `'nearest'` and `'linear'`, including
points deliberately outside the training bounding box.

**Results (accuracy/NaN numbers: confirmed bit-reproducible by the
adversarial review, re-run twice with identical `rng(7)` seeding):**
- `'nearest'`: NMSE = 0.117 (clearly beats predicting the mean), never
  produces `NaN`.
- `'linear'`, interior test points (sampled inside the training *bounding
  box*): NMSE = 0.00315 excluding `NaN`s — but **38.4% of queries returned
  `NaN` anyway**, even though every query point was inside the box.
- `'linear'`, exterior test points (deliberately outside the box): 100%
  `NaN`, confirming `interpolant.m`'s original finding.

**Caveat on the NMSE comparison (found by adversarial review, not by the
original milestone work):** the test function
`f = sin(x1) + 0.5cos(2x2) + x3*x4 - 0.3x5^2 + log1p(x6^2)` is not the
"fair, generic 6-D test" its original comment claimed. A Monte Carlo
variance decomposition (200k samples) shows `x3*x4` and `x5^2` alone
explain 97.9% of `f`'s variance — 3 of 6 input dimensions barely affect the
output. So the NMSE numbers above (0.117 vs. 0.00315) are softer evidence
than they look; they don't invalidate the milestone's actual purpose, since
**nearest-neighbor distance and convex-hull geometry — what the NaN-rate
finding depends on — use all 6 dimensions regardless of which ones drive
`f`'s value.** Also: linear's 0.00315 is computed only over the ~61.6% of
points that didn't return `NaN` (a non-random, likely-denser subset) — it
should not be read as "linear is ~37x more accurate than nearest in
general."

**Core finding, independently stress-tested and found robust (adversarial
review ran additional seeds/sizes beyond the script's own):** the convex
hull of a random point cloud in 6-D covers much less volume than its
bounding box (curse-of-dimensionality) — so "inside the box" does not mean
"inside the training convex hull." Interior-box NaN rate: 38.4% at
`n_train=2000` (this script), 40.6% at an alternate seed, 50.4% at
`n_train=1000`, 30.4% at `n_train=4000` — monotonically improving with more
training data, as expected, but every value comfortably clears the 10%
verdict threshold by 3x or more. Not a seed/size artifact. The review also
found the exterior-point test understates the fragility, if anything: even
a single input dimension pushed just 2% beyond the box edge (all other 5
dimensions solidly interior) was already 100% `NaN`.

**A claim in an earlier draft of this entry was wrong, and has been
corrected.** The original text attributed the 170.6s vs. 37.2s timing gap
between the two `'linear'` calls to "one-time Qhull warm-up cost,"
explicitly flagged at the time as untested ("not re-tested further here").
The adversarial review ran the direct, cheap test that claim invited:
calling `'linear'` on exterior points first, then interior points second
(reversed order). Result: exterior-first stayed fast, interior-second
stayed slow — the opposite of what warm-up predicts. Independently
re-confirmed in this session (exterior-first / interior-second,
single run): exterior 35.1s, interior 165.0s — same pattern. **Corrected
explanation: the cost tracks query content, not call order** — locating
which simplex contains an interior point requires a real search through the
triangulation (expensive), while rejecting an exterior point as outside the
convex hull short-circuits early (cheap). This doesn't change the
downstream conclusion (`griddatan` is expensive regardless), but the
written reasoning was factually wrong and is corrected here rather than
left standing. Also noted: timing numbers themselves are not
run-to-run reproducible (±20-25% system-load variance observed on rerun) —
unlike the accuracy/NaN numbers, they should be read as order-of-magnitude,
not precise.

**Decision:** `'nearest'` is the only method carried forward past this
milestone — matching `interpolant.m`'s original choice, now backed by
quantified reliability (no failure mode, robust across seeds/sizes) and
cost. A **hybrid** (try `'linear'`, fall back to `'nearest'` on `NaN`) was
not seriously evaluated here and is not ruled out — deferred to
Milestone 5/6: worth revisiting only if `'nearest'`-only accuracy proves
insufficient on the real validation set, given the NMSE-comparison caveat
above means linear's apparent accuracy edge isn't yet trustworthy enough to
justify the added complexity on its own.

**Other open item noted, not yet resolved:** `interpolant.m` calls
`griddatan` with explicit `qhull_opts = {'Qt','Qbb','Qc','Qx'}`; this
script used MATLAB's defaults. Convex-hull membership shouldn't depend on
this, but Milestone 6 should either match `interpolant.m`'s exact call
signature or document the deviation deliberately.

**Adversarial review:** complete — PASS WITH CAVEATS. Confirmed all
accuracy/NaN numbers bit-reproducible; found and refuted the warm-up
explanation with a direct test; found the test-function dimension-balance
issue; ran additional seed/size stress tests strengthening the core
finding; flagged the hybrid-fallback gap. No scope violations.

---

## Milestones 2+3 — Data and scaling parity

**Status:** PASS

**What was done:** wrote and ran `s00_data_and_scaling_parity.m`, which
loads the same training/validation CSVs the GP pipeline uses and proves
(not assumes) that a `rng(42); randperm(100000,3000)` reconstruction
matches the real GP training corpus, by comparing directly against an
already-trained, saved GP model object rather than trusting folder names or
saved variables (both were shown unreliable at Milestone 1).

**A real bug found and fixed, unrelated to the parity logic itself:**
MATLAB's `run()` fails on any `.m` filename starting with a digit, with a
misleading error ("Invalid text character... unsupported symbol") that
has nothing to do with the file's actual content. Confirmed by isolating it
with a trivial two-line test file: `00_test.m` (content: `disp(1+1)`) fails
with that exact error; `s00_test.m`, identical content, runs fine. This
invalidated the `00_`/`01_`/... naming convention from the Milestone 1
README — renamed the whole planned file layout to an `s##_` prefix
(`s00_`, `s01_`, ...) and fixed the one file that existed
(`00_data_and_scaling_parity.m` → `s00_data_and_scaling_parity.m`).

**A real data-provenance failure found and fully root-caused — the 1000-row
"full grid search" corpus (`saved_models/`) does not trace back to the
current training CSV; it predates it entirely.** Reconstructing a 1000-row
subset via the same `rng(42)`+`randperm` recipe and comparing to that
corpus's saved `model.X`:
- Per-row diffs ranged ~1.1 to ~88 (not floating-point noise — genuinely
  different values); none of the saved model's 1000 training rows appear
  anywhere in the current 100,000-row superset CSV, checked via exact row
  match against the *entire* superset.
- The adversarial review for this milestone pushed further rather than
  accepting "doesn't match, cause unknown": it tested alternate seeds (0-50,
  none matched), the other `data_bias_e5_*.csv` files (no match), and the
  project's older `data/data_baseline_models/e3_data/bias_train_e3.csv`
  (1000 rows, 8 columns) — **exact match, diff = 0**. Independently
  re-verified in this session. Corroborated via `git log`: `data_bias_e5_1.csv`
  (the 100k superset) was first added in commit `b61734c`, 2026-10-02
  ("Changed GP to run off of subset of 100k point dataset"); the
  `saved_models/` files are timestamped 2026-09-22, ten days earlier.
- Conclusion, now confirmed rather than inferred: those models were trained
  directly on the pre-100k-superset, e3-scale pipeline (`bias_train_e3.csv`),
  an earlier GP script version, not "stale data of unknown origin." `s00_data_and_scaling_parity.m`'s
  Check 1 updated to report this directly (diff = 0 against
  `bias_train_e3.csv`) rather than leaving it as an open question.

**The method itself is proven correct, not just assumed, via the 3000-row
corpus:** using the identical reconstruction approach, both the training
design matrix and the reproduced validation predictions (flux and emc)
matched the real `saved_models_superset_3000/` corpus **exactly — diff =
0**, not just within tolerance. This simultaneously proves: the training
CSV, the seed/subsetting recipe, the validation CSV, and the manual
target-standardization formula (mean/std of the training target column) are
all byte-identical to what the GP pipeline actually used.

**Consequence — simplified the project's parity plan:** since the 1000-row
corpus can't be verified and the 3000-row corpus is proven exact, this
project now uses `dataset_size = 3000` as its **sole** parity anchor for
both the NMSE-comparison stage and the parameter-recovery stage (dropping
the earlier two-size plan from Milestone 1). `README.md`'s "Dataset-size
parity" section updated accordingly. This means any Milestone 12 head-to-head
claim against the GP side's 1000-combination grid-search NMSE comparison
(`visualize_gp.m`) carries an explicit caveat: that specific GP corpus rests
on unverifiable training data. Reported as a limitation, not silently
worked around — regenerating it is the user's call, not assumed.

**Adversarial review (independent subagent, ran the script itself rather
than trusting the log):** confirmed all four diffs exactly as claimed;
attempted to break Check 2/3 (self-comparison, standardization artifacts,
row-order coincidence, degenerate columns) and could not — called them
"genuinely bulletproof." Found the digit-prefix `run()` bug was real but
characterized slightly too narrowly (it's not `run()`-specific — a bare
`matlab -batch "00_test"` fails identically, since the real cause is
MATLAB parsing any bare invocation as an identifier expression, and leading
digits aren't valid identifiers) — corrected in README. Found the Check-1
"stale/orphaned, cause unknown" conclusion was reached one step short of
fully root-causing it, then did that work itself (the `bias_train_e3.csv`
match and the git-log corroboration above) — independently re-verified by
me afterward, confirmed exact. Also suggested: document why checking one
`.mat` file per folder is sufficient (added as a comment), and add a
shape-assertion on the loaded CSV (added). Both applied. No scope
violations found.

---

## Milestone 1 — Scaffolding

**Status:** PASS (after corrections below)

**What was done:**
- Created `griddatan/` with `saved_models/`, `figures/`, `results/` subfolders,
  plus a scoped `.gitignore` (large artifacts in `saved_models/` and
  regenerated `figures/*.png` excluded; `results/` CSVs intentionally kept,
  matching the GP side's convention of committing comparison CSVs).
- Wrote `README.md`: purpose, scope boundary, reproducibility rules, known
  risks, planned file layout.
- No modeling code yet.

**Adversarial review (independent subagent, read-only, fact-checked every
claim against the actual files on disk):** found one central factual error
and several smaller gaps. Full findings below, with corrections applied.

**Error found and corrected — dataset-size parity claim.** The original
README asserted `dataset_size = 5000` for the NMSE-comparison stage as
settled fact. This was wrong, and traced back to me reading a stale version
of `gaussian_process.m` earlier in the session (the file had since changed
on disk, which I'd been warned about and didn't re-check before writing the
README — the exact overconfidence trap this review process exists to
catch). Follow-up: rather than just reframing this as an "open question" per
the review's suggestion, I resolved it directly by loading the actual saved
GP model objects in MATLAB (`stzd_flux_gprModel.NumObservations`), since
folder names and even the saved `dataset_size` variable turned out to be
unreliable (the oldest artifacts don't store it at all). Ground truth,
confirmed:

| Folder | Files | Confirmed training rows | Used by |
|---|---|---|---|
| `saved_models/` (unsuffixed) | 1000 (full grid) | 1000 | `visualize_gp.m` |
| `saved_models_superset_3000/` | 307 | 3000 | `synthetic_parameter_recovery_gp.m` |
| `saved_models_superset_5000/` | 10 (`best_models` shortlist) | 5000 | orphaned — nothing currently loads it |

README corrected: NMSE-comparison stage now mirrors `dataset_size = 1000`
(not 5000), parameter-recovery stage still mirrors `dataset_size = 3000`
(this half was already correct).

**Other review findings, applied:**
- Shared scaling file path was unspecific (`flux_scaling_params.mat` with no
  path) — corrected to `../parameter_recovery/flux_scaling_params.mat`
  (confirmed that's where it actually lives, not
  `gaussian_process-interpolant/`).
- `.gitignore` was missing for this folder's own large runtime artifacts —
  added, scoped to `griddatan/` only.
- Bootstrap parity (`n_boot = 1000`, `boot_size = 1000`, matching
  `gaussian_process.m`) was implicit — now stated explicitly as a commitment
  for Milestone 6.
- The "one old griddatan run" prior-art note was corrected to "two runs"
  (`FM_nearest.mat` and the mangled `FM_nearest_✐.mat`, different NMSE
  values each), and softened: their source script's own comment says "5000
  points" while its code sets `training_dataset_size = 10000` — an internal
  inconsistency, so the training size behind those specific NMSE numbers
  isn't fully trustworthy either. Both the NMSE values and the real training
  size get independently re-verified at Milestone 6.

**Prior-art findings carried forward (still to be independently re-verified
at Milestone 6, not assumed to hold):**
- The two old `'nearest'`-method runs report flux NMSE 1.65–1.69 and emc NMSE
  1.86–1.90 — worse than predicting the mean (NMSE = 1).
- `interpolant.m`'s `sprintf("FM_%s_%s", FitMethodList(i), training_dataset_size)`
  applies `%s` to a numeric value, which MATLAB converts via character code
  rather than text — confirmed on disk as the literal filename
  `FM_nearest_✐.mat` (`char(10000)` = `✐`). Will not be repeated here.
- `interpolant.m` found `'linear'` griddatan produces `NaN` outside the
  training convex hull and dropped it in favor of `'nearest'`.

**Scope check:** confirmed via `git status` at the `NTM_ML` root — only
`griddatan/` is untracked/new; nothing created or modified outside it.

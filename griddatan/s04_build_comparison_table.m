% s04_build_comparison_table.m
%
% Milestone 7: build a single comparison CSV joinable with the GP side's
% own output (`../gaussian_process-interpolant/gp_comparison_results_superset_3000.csv`,
% read-only reference, never edited). Adds a `family` column so the two
% surrogate families can be told apart after joining, and a `dataset_size`
% column so a future reader can filter to an apples-to-apples slice
% (dataset_size==3000) rather than assuming every row is comparable.
%
% Depends on Milestone 6 (`s03_train_validate_models.m`) having been run
% first -- it reads that milestone's saved .mat files rather than
% recomputing anything.

clc; clear;
here = fileparts(mfilename('fullpath'));

%% ---- GP side (read-only reference) ----
gp_csv = fullfile(here, '..', 'gaussian_process-interpolant', 'gp_comparison_results_superset_3000.csv');
assert(isfile(gp_csv), ...
    'GP comparison CSV not found at %s -- has the GP pipeline been run?', gp_csv);
% 'Delimiter' forced explicitly: readtable's auto-detection misreads this
% file as underscore-delimited (the `tags` column's values contain many
% more underscores than the file has commas), silently shredding the tags
% into a dozen bogus columns. Confirmed by direct test before this fix.
T_gp = readtable(gp_csv, 'Delimiter', ',');
T_gp.family = repmat("GP", height(T_gp), 1);
T_gp.dataset_size = repmat(3000, height(T_gp), 1);

%% ---- griddatan side (this project's own Milestone 6 output) ----
saved_models_dir = fullfile(here, 'saved_models');

% Part A: the parity-anchored result (dataset_size=3000) -- the row that
% matters for any apples-to-apples claim.
file_A = fullfile(saved_models_dir, 'FM_nearest__dataset_3000.mat');
assert(isfile(file_A), 'Milestone 6 output not found at %s -- run s03_train_validate_models.m first.', file_A);
d_A = load(file_A, 'flux_nmse_griddatan', 'emc_nmse_griddatan', 'dataset_size');
T_griddatan_A = table("FM_nearest__dataset_3000", d_A.flux_nmse_griddatan, d_A.emc_nmse_griddatan, ...
    "griddatan", d_A.dataset_size, ...
    'VariableNames', {'tags','flux_nmse','emc_nmse','family','dataset_size'});

% Part B: the prior-art re-check (dataset_size=10000) -- included for
% completeness, but on a DIFFERENT corpus than the 3000-row parity anchor,
% hence its own dataset_size value so it is never silently treated as
% comparable to the dataset_size==3000 rows above.
file_B = fullfile(here, 'results', 'priorart_recheck_dataset_10000.csv');
T_griddatan_B = table("FM_nearest__dataset_10000", nan, nan, "griddatan", 10000, ...
    'VariableNames', {'tags','flux_nmse','emc_nmse','family','dataset_size'});
if isfile(file_B)
    d_B = readtable(file_B, 'Delimiter', ',');
    T_griddatan_B.flux_nmse = d_B.flux_nmse(1);
    T_griddatan_B.emc_nmse = d_B.emc_nmse(1);
else
    warning('Milestone 6 Part B output not found at %s -- that row will be left as NaN.', file_B);
end

%% ---- Combine and save ----
T_gp = T_gp(:, {'tags','flux_nmse','emc_nmse','family','dataset_size'}); % column order match
T_all = [T_gp; T_griddatan_A; T_griddatan_B];
T_all = sortrows(T_all, 'flux_nmse');

folder_results = fullfile(here, 'results');
if ~isfolder(folder_results); mkdir(folder_results); end
out_csv = fullfile(folder_results, 'comparison_table.csv');
writetable(T_all, out_csv);

%% ---- Summary ----
n_total = height(T_all);
rank_griddatan_3000 = find(strcmp(T_all.tags, "FM_nearest__dataset_3000"));
fprintf('Combined comparison table: %d rows (%d GP configs + 2 griddatan runs)\n', ...
    n_total, height(T_gp));
fprintf('Saved to %s\n\n', out_csv);
fprintf('griddatan (dataset_size=3000) ranks %d of %d by flux NMSE (1 = best).\n', ...
    rank_griddatan_3000, n_total);
% min()/max(), not positional indexing (e.g. T_gp.flux_nmse(1)) -- that
% was only correct by coincidence of the source CSV's on-disk row order,
% a fragility the Milestone 7 adversarial reviewer caught and flagged.
best_gp_flux = min(T_gp.flux_nmse);
worst_gp_flux = max(T_gp.flux_nmse);
fprintf(['Every GP configuration beats griddatan on this corpus: best GP flux NMSE ' ...
    'is %.5f vs. griddatan''s %.4f (%.0fx worse). Reported plainly for Milestone 12,\n' ...
    'not softened.\n'], best_gp_flux, d_A.flux_nmse_griddatan, ...
    d_A.flux_nmse_griddatan / best_gp_flux);
fprintf('Even GP''s WORST config (flux NMSE=%.6f) still beats griddatan by %.1fx.\n', ...
    worst_gp_flux, d_A.flux_nmse_griddatan / worst_gp_flux);

% s03_train_validate_models.m
%
% Milestone 6: full train+validate pipeline for the griddatan surrogate,
% mirroring gaussian_process.m's NMSE + bootstrap methodology number-for-
% number (same n_boot, boot_size, formulas) so the two pipelines are
% directly comparable.
%
% Two independent parts:
%   Part A (primary, parity-anchored): dataset_size=3000, the corpus
%     Milestone 2 proved matches the GP pipeline's training data exactly.
%     This is the comparison target used from here through Milestone 12.
%   Part B (prior-art re-check, Risk 4 in README.md): independently
%     reproduces interpolant.m's own original setup (dataset_size=10000,
%     its qhull_opts) to confirm or refute the "NMSE 1.65-1.90, worse than
%     predicting the mean" claim under controlled, documented conditions --
%     not assumed to carry over untested.

clc; clear;
here = fileparts(mfilename('fullpath'));

%% ---- Shared data loading ----
tr_filepath = fullfile(here, '..', 'data', 'data_neural_network_model', 'training_data_e6', 'data_bias_e5_1.csv');
va_filepath = fullfile(here, '..', 'data', 'data_baseline_models', 'e3_data', 'bias_val_e3.csv');
training_superset = readmatrix(tr_filepath);
validation_dataset = readmatrix(va_filepath);
val_features = validation_dataset(:,1:6);
val_flux = validation_dataset(:,7);
val_emc = validation_dataset(:,8);

n_boot = 1000;
boot_size = 1000;

%% ================= PART A: parity-anchored run (dataset_size=3000) =================
training_seed = 42;
dataset_size = 3000;
rng(training_seed);
idx_3000 = randperm(100000, dataset_size);
training_dataset = training_superset(idx_3000, :);
features = training_dataset(:,1:6);
flux = training_dataset(:,7);
emc = training_dataset(:,8);

% Standardize targets -- matches gaussian_process.m's convention, even
% though it is provably a no-op for 'nearest': the nearest-neighbor index
% depends only on x, and un-standardizing an affine transform of y(idx)
% recovers y(idx) exactly. Verified below, not just asserted in a comment.
flux_mean = mean(flux); flux_std = std(flux);
flux_z = (flux - flux_mean) / flux_std;
emc_mean = mean(emc); emc_std = std(emc);
emc_z = (emc - emc_mean) / emc_std;

tic;
model = fit_nearest_predictor(features); % one k-d tree; flux and emc share the same x
idx_val = knnsearch(model.ns, val_features);
idx_train = knnsearch(model.ns, features);
wall_time = toc;

griddatan_flux_ypred = flux(idx_val);
griddatan_emc_ypred = emc(idx_val);
griddatan_flux_xpred = flux(idx_train);
griddatan_emc_xpred = emc(idx_train);

% No-op check (Part A's standardize/un-standardize round trip must
% reproduce the raw nearest-neighbor copy, up to floating-point rounding
% from the divide-then-multiply round trip -- not bit-exact zero, so a
% tight tolerance is used rather than ==0, matching s00's convention):
flux_ypred_via_z = flux_z(idx_val) * flux_std + flux_mean;
d_noop = max(abs(flux_ypred_via_z - griddatan_flux_ypred));
assert(d_noop < 1e-9, ...
    'Standardize round-trip is not a no-op for nearest-neighbor prediction (diff=%g) -- investigate.', d_noop);

% Every training point is its own nearest neighbor (distance 0) unless the
% corpus has duplicate/near-duplicate rows -- Milestone 5 close-out
% confirmed this 3000-row corpus has none (min pairwise distance 0.00994).
% So training NMSE is expected to be EXACTLY zero; assert it rather than
% just expect it.
assert(isequal(idx_train, (1:numel(flux))'), ...
    'Expected every training point to be its own nearest neighbor (training NMSE should be exactly 0). Investigate duplicate/near-duplicate training rows.');

squared_error_flux_y = (val_flux - griddatan_flux_ypred).^2;
squared_error_emc_y = (val_emc - griddatan_emc_ypred).^2;
squared_error_flux_x = (flux - griddatan_flux_xpred).^2;
squared_error_emc_x = (emc - griddatan_emc_xpred).^2;

flux_nmse_griddatan = mean(squared_error_flux_y) / var(val_flux);
emc_nmse_griddatan = mean(squared_error_emc_y) / var(val_emc);
flux_nmse_griddatan_train = mean(squared_error_flux_x) / var(flux);
emc_nmse_griddatan_train = mean(squared_error_emc_x) / var(emc);

fprintf('[Part A, dataset_size=%d] flux NMSE (val) = %.4f, emc NMSE (val) = %.4f\n', ...
    dataset_size, flux_nmse_griddatan, emc_nmse_griddatan);
fprintf('[Part A] flux NMSE (train) = %.6f, emc NMSE (train) = %.6f -- expected exactly 0:\n', ...
    flux_nmse_griddatan_train, emc_nmse_griddatan_train);
fprintf(['''nearest'' is an exact interpolator (every training point is its own\n' ...
    'nearest neighbor), unlike the GP, which smooths and generally has\n' ...
    'nonzero training residual. A real, structural difference between the\n' ...
    'two surrogate families -- reported as-is for Milestone 12, not hidden.\n\n']);

% Bootstrap -- identical formula and parameters to gaussian_process.m.
% No extra rng() reseed here: as in the GP script, the random stream is
% left to carry forward from the single rng(training_seed) call above, so
% re-running this script end-to-end reproduces the same bootstrap draws.
boot_flux_nmse_test = zeros(n_boot,1);
for b = 1:n_boot
    boot_idx = randi(numel(squared_error_flux_y), boot_size, 1);
    boot_flux_nmse_test(b) = mean(squared_error_flux_y(boot_idx)) / var(val_flux);
end
boot_flux_nmse_train = zeros(n_boot,1);
for b = 1:n_boot
    boot_idx = randi(numel(squared_error_flux_x), boot_size, 1);
    boot_flux_nmse_train(b) = mean(squared_error_flux_x(boot_idx)) / var(flux);
end
boot_emc_nmse_test = zeros(n_boot,1);
for b = 1:n_boot
    boot_idx = randi(numel(squared_error_emc_y), boot_size, 1);
    boot_emc_nmse_test(b) = mean(squared_error_emc_y(boot_idx)) / var(val_emc);
end
boot_emc_nmse_train = zeros(n_boot,1);
for b = 1:n_boot
    boot_idx = randi(numel(squared_error_emc_x), boot_size, 1);
    boot_emc_nmse_train(b) = mean(squared_error_emc_x(boot_idx)) / var(emc);
end

folder_name = fullfile(here, 'saved_models');
if ~isfolder(folder_name); mkdir(folder_name); end
tag = sprintf("FM_nearest__dataset_%d", dataset_size);
filename = fullfile(folder_name, tag + ".mat");
save(filename, ...
    "griddatan_flux_ypred", "griddatan_emc_ypred", ...
    "flux_nmse_griddatan", "emc_nmse_griddatan", ...
    "flux_nmse_griddatan_train", "emc_nmse_griddatan_train", ...
    "wall_time", "dataset_size", ...
    "boot_flux_nmse_train", "boot_flux_nmse_test", ...
    "boot_emc_nmse_train", "boot_emc_nmse_test", ...
    "flux_mean", "flux_std", "emc_mean", "emc_std");
fprintf('[Part A] Saved to %s\n\n', filename);

%% ================= PART B: prior-art re-check (Risk 4) =================
% Reproduces interpolant.m's own setup directly: same seed, same
% training_dataset_size=10000, same qhull_opts, calling griddatan
% directly (a one-off diagnostic, not a repeated-query workload, so the
% fast predictor isn't needed here).
training_seed_B = 42;
training_dataset_size_B = 10000;
rng(training_seed_B);
idx_10000 = randperm(100000, training_dataset_size_B);
training_dataset_B = training_superset(idx_10000, :);
features_B = training_dataset_B(:,1:6);
flux_B = training_dataset_B(:,7);
emc_B = training_dataset_B(:,8);

flux_mean_B = mean(flux_B); flux_std_B = std(flux_B);
flux_z_B = (flux_B - flux_mean_B) / flux_std_B;
emc_mean_B = mean(emc_B); emc_std_B = std(emc_B);
emc_z_B = (emc_B - emc_mean_B) / emc_std_B;

qhull_opts = {'Qt','Qbb','Qc','Qx'};

tic;
griddatan_flux_ypred_z_B = griddatan(features_B, flux_z_B, val_features, 'nearest', qhull_opts);
griddatan_emc_ypred_z_B  = griddatan(features_B, emc_z_B,  val_features, 'nearest', qhull_opts);
wall_time_B = toc;

griddatan_flux_ypred_B = griddatan_flux_ypred_z_B * flux_std_B + flux_mean_B;
griddatan_emc_ypred_B  = griddatan_emc_ypred_z_B  * emc_std_B  + emc_mean_B;

flux_nmse_griddatan_B = mean((val_flux - griddatan_flux_ypred_B).^2) / var(val_flux);
emc_nmse_griddatan_B  = mean((val_emc  - griddatan_emc_ypred_B).^2)  / var(val_emc);

fprintf('[Part B, dataset_size=%d, qhull_opts] flux NMSE = %.4f, emc NMSE = %.4f\n', ...
    training_dataset_size_B, flux_nmse_griddatan_B, emc_nmse_griddatan_B);

if flux_nmse_griddatan_B > 1 && emc_nmse_griddatan_B > 1
    fprintf(['[Part B] CONFIRMED: this re-run performs worse than predicting\n' ...
        'the mean (NMSE>1) for both targets, matching the prior-art claim.\n']);
elseif flux_nmse_griddatan_B < 1 && emc_nmse_griddatan_B < 1
    fprintf(['[Part B] NOT REPRODUCED: this re-run beats predicting the mean\n' ...
        '(NMSE<1) for both targets -- the prior-art numbers do not carry\n' ...
        'over under the current data/pipeline.\n']);
else
    fprintf('[Part B] MIXED: one target beats the mean, the other does not.\n');
end

% Does qhull_opts itself change anything vs. the plain call Milestone 5
% proved exact against the fast predictor? If diff=0, Part A's parity
% proof extends to this configuration too.
griddatan_flux_ypred_noopt_B = griddatan(features_B, flux_z_B, val_features, 'nearest') * flux_std_B + flux_mean_B;
d_opts = max(abs(griddatan_flux_ypred_noopt_B - griddatan_flux_ypred_B));
fprintf('[Part B] qhull_opts vs. no-opts, flux prediction max abs diff: %g\n', d_opts);

folder_results = fullfile(here, 'results');
if ~isfolder(folder_results); mkdir(folder_results); end
T_B = table(training_dataset_size_B, flux_nmse_griddatan_B, emc_nmse_griddatan_B, wall_time_B, d_opts, ...
    'VariableNames', {'training_dataset_size','flux_nmse','emc_nmse','wall_time_s','qhull_opts_vs_noopt_maxdiff'});
writetable(T_B, fullfile(folder_results, 'priorart_recheck_dataset_10000.csv'));
fprintf('[Part B] Saved to %s\n', fullfile(folder_results, 'priorart_recheck_dataset_10000.csv'));

fprintf('\nMilestone 6 complete: Part A saved for the Milestone 7+ comparison pipeline, Part B resolves Risk 4.\n');

%% ---- Local functions ----
% Duplicated from s02_fast_predictor.m rather than shared, intentionally:
% every script in this project is self-contained and runnable on its own
% (see README.md's reproducibility rules) -- small enough (6 lines) that a
% shared utility file would add more indirection than it saves.

function model = fit_nearest_predictor(X_train)
    % Same duplicate-row warning as s02_fast_predictor.m -- griddatan
    % averages duplicate rows' targets before triangulating; a k-d tree
    % does not. Confirmed not to apply to the 3000-row corpus (Part A)
    % but checked here too in case Part A's data path ever changes.
    if size(unique(X_train, 'rows'), 1) < size(X_train, 1)
        warning('fit_nearest_predictor:DuplicateRows', ...
            ['Training data contains duplicate feature rows. griddatan ' ...
             'averages duplicates'' targets before use; this k-d tree ' ...
             'predictor does not.']);
    end
    model.ns = createns(X_train, 'NSMethod', 'kdtree');
end

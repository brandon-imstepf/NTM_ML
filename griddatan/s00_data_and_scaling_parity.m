% s00_data_and_scaling_parity.m
%
% Milestones 2+3: prove this project loads byte-identical training and
% validation data, and identical target standardization, to the GP
% pipeline -- before any griddatan model is built on top of it.
%
% Verification strategy: folder names and even saved metadata turned out to
% be unreliable during scaffolding (see STATUS.md, Milestone 1), so this
% does not trust either. Instead it loads an *already-trained* GP model
% directly and compares its stored training design matrix (model.X) against
% our own reconstructed subset, row for row. Then it runs our own
% validation features through that same model and checks we reproduce its
% saved predictions exactly -- which proves validation-data parity and
% target-standardization parity in one shot, since both would have to be
% right for predictions to match.

clc; clear;
here = fileparts(mfilename('fullpath'));

%% Paths -- anchored to this script's own folder, never to pwd
gp_dir = fullfile(here, '..', 'gaussian_process-interpolant');
pr_dir = fullfile(here, '..', 'parameter_recovery');

tr_filepath = fullfile(here, '..', 'data', 'data_neural_network_model', 'training_data_e6', 'data_bias_e5_1.csv');
va_filepath = fullfile(here, '..', 'data', 'data_baseline_models', 'e3_data', 'bias_val_e3.csv');

training_superset  = readmatrix(tr_filepath);   % 100k rows
validation_dataset = readmatrix(va_filepath);
assert(size(training_superset,1) == 100000, ...
    'Expected a 100,000-row training superset -- CSV shape has changed.');

%% Reproduce GP's verified training subset
% GP resets rng(42) immediately before each randperm call, in an independent
% script run -- not inherited from a prior draw in the same session. The
% 3000-row subset below is the ONLY GP training corpus whose provenance
% could be independently proven (see Check 1 below); it is used as the sole
% parity anchor for this entire project, for both the NMSE comparison stage
% and the parameter-recovery stage.
training_seed = 42;
rng(training_seed); idx_3000 = randperm(100000, 3000);
subset_3000 = training_superset(idx_3000, :);

features_3000 = subset_3000(:,1:6);
flux_3000 = subset_3000(:,7);
emc_3000  = subset_3000(:,8);

tol = 1e-9;   % tight, but allows for benign floating-point summation order

%% Check 1 -- diagnostic only: root-cause confirmation that the unsuffixed
%% saved_models/ corpus (the "real" 1000-combination grid search) predates
%% the current 100k-superset pipeline entirely.
% NON-FATAL. It does NOT trace to data_bias_e5_1.csv (max per-row diff ~88;
% none of its 1000 training rows appear anywhere in the 100k-row superset,
% confirmed by exact-match search, not just within the expected subset).
% Root cause identified (not just ruled out): those models were trained
% directly on data/data_baseline_models/e3_data/bias_train_e3.csv (exact
% diff = 0, see below) -- the pre-100k-superset, e3-scale pipeline, which
% predates commit b61734c (2026-10-02, "Changed GP to run off of subset of
% 100k point dataset"). The saved_models/ files are timestamped 2026-09-22,
% ten days before that commit. So this is not an unexplained anomaly: it is
% an earlier pipeline version's output, confirmed, and correctly excluded
% from this project's comparison targets since it used different training
% data than everything else here.
rng(training_seed); idx_1000 = randperm(100000, 1000); %#ok<NASGU>
ref_1000 = load(fullfile(gp_dir, 'saved_models', ...
    'FM_exact__BF_constant__KF_ARDExponential__OPT_fmincon.mat'), 'stzd_flux_gprModel');
e3_train_filepath = fullfile(here, '..', 'data', 'data_baseline_models', 'e3_data', 'bias_train_e3.csv');
e3_train = readmatrix(e3_train_filepath);
d1 = max(abs(ref_1000.stzd_flux_gprModel.X(:) - reshape(e3_train(:,1:6), [], 1)));
fprintf(['[Check 1 -- DIAGNOSTIC, NON-FATAL] saved_models/ corpus vs. ' ...
    'bias_train_e3.csv (the pre-100k-superset pipeline''s training file), ' ...
    'max abs diff: %g (0 confirms root cause; that corpus is not used as ' ...
    'a parity target here either way)\n'], d1);

%% Check 2 -- 3000-row subset matches the parameter-recovery corpus (saved_models_superset_3000/)
% Only one file from this folder is checked. That's sufficient, not just
% convenient: gaussian_process.m computes features/flux_z/emc_z ONCE,
% outside its per-hyperparameter-combination loop, then reuses them for
% every fitrgp() call -- so every .mat file in this folder shares a
% bit-identical training matrix by construction, not by chance.
ref_3000 = load(fullfile(gp_dir, 'saved_models_superset_3000', ...
    'FM_exact__BF_constant__KF_ARDExponential__OPT_fmincon.mat'), ...
    'stzd_flux_gprModel', 'stzd_emc_gprModel', 'gpr_flux_ypred', 'gpr_emc_ypred');
d2 = max(abs(ref_3000.stzd_flux_gprModel.X(:) - features_3000(:)));
fprintf('[Check 2] 3000-row subset vs. saved GP model X, max abs diff: %g\n', d2);
assert(d2 < tol, 'Training-set reconstruction does NOT match the GP 3000-row corpus.');

%% Check 3 -- end-to-end: validation data + target standardization, flux and emc
% GP standardizes its *target* (flux/emc) manually, outside fitrgp, so it's
% not stored inside the model object. Reconstruct it the same way GP does
% and prove it's right by reproducing GP's own saved validation predictions.
flux_mean_3000 = mean(flux_3000); flux_std_3000 = std(flux_3000);
emc_mean_3000  = mean(emc_3000);  emc_std_3000  = std(emc_3000);

val_features = validation_dataset(:,1:6);

[flux_pred_z, ~] = predict(ref_3000.stzd_flux_gprModel, val_features);
our_flux_ypred = flux_pred_z * flux_std_3000 + flux_mean_3000;
d3_flux = max(abs(our_flux_ypred - ref_3000.gpr_flux_ypred));
fprintf('[Check 3a] reconstructed flux prediction vs. saved gpr_flux_ypred, max abs diff: %g\n', d3_flux);
assert(d3_flux < tol, 'Flux validation parity / target-standardization reconstruction failed.');

[emc_pred_z, ~] = predict(ref_3000.stzd_emc_gprModel, val_features);
our_emc_ypred = emc_pred_z * emc_std_3000 + emc_mean_3000;
d3_emc = max(abs(our_emc_ypred - ref_3000.gpr_emc_ypred));
fprintf('[Check 3b] reconstructed emc prediction vs. saved gpr_emc_ypred, max abs diff: %g\n', d3_emc);
assert(d3_emc < tol, 'EMC validation parity / target-standardization reconstruction failed.');

%% Shared scaling constants used later (parameter-recovery stage, Milestone 9+)
% Loaded here only to confirm they're reachable from this folder; not
% recomputed, not otherwise used until the recovery scripts.
flux_scaling = load(fullfile(pr_dir, 'flux_scaling_params.mat')); %#ok<NASGU>
emc_scaling  = load(fullfile(pr_dir, 'emc_scaling_params.mat'));  %#ok<NASGU>
fprintf('Shared flux/emc scaling constants loaded from parameter_recovery/ (used from Milestone 9 onward).\n');

fprintf('\nAll parity checks PASSED.\n');

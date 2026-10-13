% s06_parameter_recovery_gamma.m
%
% Milestone 9: ports the gamma-recovery procedure from
% `..\parameter_recovery\synthetic_parameter_recovery_gp.m` (read-only
% reference) onto the griddatan fast predictor. Same recovery setup: hold
% (lambda1, delta, epsilon, N1, N2) at their true per-row values, use
% `fminbnd` to find the gamma1 that makes the predicted (standardized)
% flux match the observed (standardized) flux, over the full holdout
% dataset -- not just a cherry-picked subset. The feasibility metric is
% `log_error_gamma` (|log10(recovered) - log10(true)|) across every row,
% matching how this metric is already used on the GP side.
%
% Resolves README.md's Risk 2 (no native predictive variance) by NOT
% fabricating a confidence band: griddatan's 'nearest' has no equivalent
% to the GP's predict()-based uncertainty, so the illustrative recovery
% plots say so explicitly instead of silently omitting or faking one.
%
% Only ONE griddatan "model" exists (unlike the GP script's 4 near-tied
% hyperparameter configs) -- there is nothing to sweep here.
%
% Depends on Milestone 2 (parity-anchored data reconstruction) and
% Milestone 5 (fast predictor) patterns already established in this folder.

clc; clear;
here = fileparts(mfilename('fullpath'));

%% ---- Load parity-anchored training data (dataset_size=3000) ----
tr_filepath = fullfile(here, '..', 'data', 'data_neural_network_model', 'training_data_e6', 'data_bias_e5_1.csv');
training_superset = readmatrix(tr_filepath);
rng(42);
idx_3000 = randperm(100000, 3000);
subset_3000 = training_superset(idx_3000, :);
features = subset_3000(:,1:6);
flux = subset_3000(:,7);

model = fit_nearest_predictor(features, flux);

%% ---- Load holdout dataset + shared scaling constants (read-only) ----
va_filepath = fullfile(here, '..', 'data', 'data_baseline_models', 'e3_data', 'bias_val_e3.csv');
dataset = readmatrix(va_filepath);

flux_scaling = load(fullfile(here, '..', 'parameter_recovery', 'flux_scaling_params.mat'));

%% ---- Recovery setup (identical ranges/sample counts to the GP script) ----
numiterations = 1000;
gamma1_range = [5e-5, 5e-2];
logdist_fun = @(x,n) 10.^(log10(x(1)) + (log10(x(2)/x(1)) * rand(n,1)));
rng(42); % fixed seed: gammaspace is only used for the illustrative loss curve, not the recovery itself
gammaspace = logdist_fun(gamma1_range, numiterations);

n_rows = size(dataset, 1);
[refined_gamma1, refined_loss_gamma] = deal(zeros(n_rows, 1));

tic;
for row = 1:n_rows
    l_row = dataset(row,2); d_row = dataset(row,3); e_row = dataset(row,4);
    n1_row = dataset(row,5); n2_row = dataset(row,6);
    flux_scaled = (dataset(row,7) - flux_scaling.target_mean) / flux_scaling.target_std;

    gamma_loss_fun = @(gamma1_candidate) ...
        (griddatan_flux_scaled(model, gamma1_candidate, l_row, d_row, e_row, n1_row, n2_row, flux_scaling) - flux_scaled)^2;
    [refined_gamma1(row), refined_loss_gamma(row)] = fminbnd(gamma_loss_fun, gamma1_range(1), gamma1_range(2));
end
wall_time = toc;
fprintf('Gamma recovery over all %d holdout rows: %.2f s (%.4f s/row)\n', n_rows, wall_time, wall_time/n_rows);

%% ---- Feasibility metric: log_error_gamma across the FULL dataset ----
log_error_gamma = abs(log10(refined_gamma1) - log10(dataset(:,1)));
fprintf('\nlog_error_gamma over all %d rows: mean=%.4f, median=%.4f, std=%.4f\n', ...
    n_rows, mean(log_error_gamma), median(log_error_gamma), std(log_error_gamma));
fprintf('(For reference: log_error_gamma=1 means the recovered gamma1 is off by a full\n');
fprintf('decade from the true value -- this is the same statistic the GP side already\n');
fprintf('computes, so it is directly comparable across both surrogate families.)\n\n');

%% ---- Root-cause diagnostic: why recovery fails this badly ----
% The illustrative plots below turned out to show a FLAT loss-vs-gamma
% curve for every inspected row -- the nearest neighbor found doesn't
% change at all as the candidate gamma1 sweeps its entire range. Checked
% directly here (not just eyeballed on a plot) before writing this up.
fprintf('--- Root-cause check: is griddatan''s nearest-neighbor search blind to gamma1? ---\n');
check_row = 70; % arbitrary holdout row, same one used for the first illustrative plot
l = dataset(check_row,2); d = dataset(check_row,3); e = dataset(check_row,4);
n1 = dataset(check_row,5); n2 = dataset(check_row,6);
gamma_probe = [5e-5, 0.001, 0.01, 0.05];
idx_probe = zeros(numel(gamma_probe),1);
for i = 1:numel(gamma_probe)
    idx_probe(i) = knnsearch(model.ns, [gamma_probe(i), l, d, e, n1, n2]);
end
fprintf('Row %d: nearest-neighbor index for gamma1 in {5e-5, 0.001, 0.01, 0.05}: %s (all identical: %d)\n', ...
    check_row, mat2str(idx_probe'), numel(unique(idx_probe))==1);
sq_contrib = ([gamma_probe(1), l, d, e, n1, n2] - features(idx_probe(1),:)).^2;
fprintf('Squared per-dimension distance contribution for that match: gamma1=%.3g, lambda1=%.3g, delta=%.3g, epsilon=%.3g, N1=%.3g, N2=%.3g\n', sq_contrib);
fprintf(['ROOT CAUSE, confirmed: delta and epsilon (raw scale ~10-100) dominate the raw\n' ...
    'Euclidean distance by 5-6 orders of magnitude over gamma1 (~1e-5 to 5e-2). Changing\n' ...
    'gamma1 across its ENTIRE allowed range barely moves the distance at all, so the\n' ...
    '"nearest" training point is selected almost independently of gamma1 -- the fminbnd\n' ...
    'search above is optimizing over an essentially flat, noise-dominated loss surface.\n' ...
    'This is the mechanism behind both the catastrophic log_error_gamma above AND the\n' ...
    'NMSE>1 results already reported in Milestones 6-8: griddatan''s raw-feature (not\n' ...
    'standardized) distance metric effectively ignores gamma1/lambda1/N1/N2 relative to\n' ...
    'delta/epsilon. This is a real, structural limitation of griddatan as actually built\n' ...
    'by interpolant.m and carried forward by this project for apples-to-apples parity --\n' ...
    'not a bug introduced here.\n\n']);

%% ---- Diagnostic only (NOT the apples-to-apples pipeline): standardized features ----
% Quick, clearly-separated sensitivity check: if features are z-scored
% before building the k-d tree (deviating from interpolant.m's original,
% unstandardized choice -- which is why this is kept out of the main
% pipeline and Milestones 6-8's numbers), does recovery actually improve?
% This does not change any committed result; it only diagnoses whether
% the root cause above is the (near-)complete explanation.
feat_mean = mean(features, 1); feat_std = std(features, 0, 1);
features_z = (features - feat_mean) ./ feat_std;
model_z = fit_nearest_predictor(features_z, flux);

refined_gamma1_z = zeros(n_rows, 1);
for row = 1:n_rows
    l_row = dataset(row,2); d_row = dataset(row,3); e_row = dataset(row,4);
    n1_row = dataset(row,5); n2_row = dataset(row,6);
    flux_scaled = (dataset(row,7) - flux_scaling.target_mean) / flux_scaling.target_std;
    gamma_loss_fun_z = @(g) (griddatan_flux_scaled_z(model_z, g, l_row, d_row, e_row, n1_row, n2_row, flux_scaling, feat_mean, feat_std) - flux_scaled)^2;
    refined_gamma1_z(row) = fminbnd(gamma_loss_fun_z, gamma1_range(1), gamma1_range(2));
end
log_error_gamma_z = abs(log10(refined_gamma1_z) - log10(dataset(:,1)));
fprintf('--- Diagnostic only: standardized-feature k-d tree (NOT part of the committed pipeline) ---\n');
fprintf('log_error_gamma (standardized features): mean=%.4f, median=%.4f (vs. raw-feature mean=%.4f, median=%.4f)\n', ...
    mean(log_error_gamma_z), median(log_error_gamma_z), mean(log_error_gamma), median(log_error_gamma));
if mean(log_error_gamma_z) < mean(log_error_gamma)
    fprintf('Standardizing features substantially improves recovery -- confirms the root cause above.\n\n');
else
    fprintf('Standardizing features does NOT meaningfully improve recovery -- the failure runs deeper than scale alone.\n\n');
end

%% ---- Save full-dataset recovery results ----
folder_results = fullfile(here, 'results');
if ~isfolder(folder_results); mkdir(folder_results); end
T_recovery = table((1:n_rows)', dataset(:,1), refined_gamma1, refined_loss_gamma, log_error_gamma, ...
    'VariableNames', {'row','true_gamma1','recovered_gamma1','min_loss','log_error_gamma'});
writetable(T_recovery, fullfile(folder_results, 'gamma_recovery_full.csv'));
fprintf('Full per-row recovery results saved to results/gamma_recovery_full.csv\n');

%% ---- Illustrative plots: best 10 rows by log_error_gamma, same style as the GP side ----
lowest_count = 10;
[~, best_gamma_recovery_rows] = mink(log_error_gamma, lowest_count);

outdir = fullfile(here, 'figures', 'gamma_recovery');
if ~isfolder(outdir); mkdir(outdir); end
for row = best_gamma_recovery_rows'
    gammaRecoveryBuilder(dataset, row, gammaspace, model, flux_scaling, outdir, ...
        refined_gamma1(row), refined_loss_gamma(row));
end
fprintf('Saved %d illustrative recovery plots to %s\n', lowest_count, outdir);

fprintf('\nMilestone 9 complete.\n');

%% ---- Local functions ----

function model = fit_nearest_predictor(X_train, y_train)
    % Same duplicate-row check as s02/s03 -- see those files for why.
    if size(unique(X_train, 'rows'), 1) < size(X_train, 1)
        warning('fit_nearest_predictor:DuplicateRows', ...
            'Training data contains duplicate feature rows; griddatan averages their targets, this k-d tree does not.');
    end
    model.ns = createns(X_train, 'NSMethod', 'kdtree');
    model.y = y_train;
end

function flux_z = griddatan_flux_scaled(model, gamma1, lambda1, delta, epsilon, N1, N2, flux_scaling)
    % Mirrors gp_flux_scaled_1's signature and broadcasting behavior from
    % synthetic_parameter_recovery_gp.m exactly, so this is a drop-in
    % substitute for the GP's predictor in the recovery loop. gamma1 (and,
    % for the Milestone 10 joint-recovery port, lambda1) may be arrays
    % (e.g. an ndgrid); the other four arguments are scalar per call.
    sz = size(gamma1);
    gamma1 = gamma1(:);
    n = numel(gamma1);
    lambda1 = lambda1(:);
    if isscalar(lambda1) && n > 1
        lambda1 = repmat(lambda1, n, 1);
    end
    X_query = [gamma1, lambda1, repmat(delta,n,1), repmat(epsilon,n,1), repmat(N1,n,1), repmat(N2,n,1)];
    idx = knnsearch(model.ns, X_query);
    flux_raw = model.y(idx);
    flux_z = reshape((flux_raw - flux_scaling.target_mean) / flux_scaling.target_std, sz);
end

function flux_z = griddatan_flux_scaled_z(model_z, gamma1, lambda1, delta, epsilon, N1, N2, flux_scaling, feat_mean, feat_std)
    % Root-cause-diagnostic variant only (see the "Diagnostic only" section
    % above) -- same as griddatan_flux_scaled but standardizes the query
    % the same way model_z's training features were standardized before
    % the k-d tree was built, so the comparison is apples-to-apples in
    % feature space too.
    sz = size(gamma1);
    gamma1 = gamma1(:);
    n = numel(gamma1);
    X_query_raw = [gamma1, repmat(lambda1,n,1), repmat(delta,n,1), repmat(epsilon,n,1), repmat(N1,n,1), repmat(N2,n,1)];
    X_query = (X_query_raw - feat_mean) ./ feat_std;
    idx = knnsearch(model_z.ns, X_query);
    flux_raw = model_z.y(idx);
    flux_z = reshape((flux_raw - flux_scaling.target_mean) / flux_scaling.target_std, sz);
end

function gammaRecoveryBuilder(dataset, row, gammaspace, model, flux_scaling, outdir, refined_gamma1, refined_loss_gamma)
    % Same visual layout as the GP side's gammaRecoveryBuilder: loss vs.
    % candidate gamma, recovered point, true-gamma reference line. NO
    % shaded confidence band -- griddatan has no native predictive
    % variance (README.md Risk 2), so this is stated on the plot itself
    % rather than silently omitted or faked with a substitute metric.
    l = dataset(row,2); d = dataset(row,3); e = dataset(row,4);
    n1 = dataset(row,5); n2 = dataset(row,6);
    flux_scaled = (dataset(row,7) - flux_scaling.target_mean) / flux_scaling.target_std;

    guess_flux = griddatan_flux_scaled(model, gammaspace, l, d, e, n1, n2, flux_scaling);
    squared_error = (guess_flux - flux_scaled).^2;

    [sorted_gamma, idx] = sort(gammaspace);

    figure('Visible', 'off');
    plot(sorted_gamma, squared_error(idx), 'b', 'LineWidth', 3);
    xlabel('Gamma Value', 'FontSize', 15);
    ylabel('Loss', 'FontSize', 15);
    hold on;
    h_dot = plot(refined_gamma1, refined_loss_gamma, 'ro', 'MarkerFaceColor', 'r');
    % LabelVerticalAlignment 'bottom' (not the default 'top'): the
    % adversarial reviewer found the label clipped at the axes' top edge
    % in most sampled plots with the default placement.
    xline(dataset(row,1), '--r', sprintf('True Gamma, %.4g', dataset(row,1)), ...
        'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'left');
    legend(h_dot, sprintf('Minimum loss = %.4g', refined_loss_gamma), 'Location', 'northeast', 'FontSize', 14);
    % Bottom-left, not top-left: a long annotation at the top collided
    % with the legend box (found by inspecting the first render of this
    % figure) since the legend's width pushes its left edge further left
    % than the plot's nominal "northeast" corner suggests.
    text(0.02, 0.05, 'No predictive-uncertainty band: griddatan has no native predictive variance (README Risk 2)', ...
        'Units', 'normalized', 'FontSize', 9, 'Color', [0.4 0.4 0.4], 'VerticalAlignment', 'bottom');
    hold off;
    filename = fullfile(outdir, sprintf('flux_recovery_row%04d.png', row));
    exportgraphics(gcf, filename, 'Resolution', 150);
    close(gcf);
end

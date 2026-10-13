% s07_parameter_recovery_joint.m
%
% Milestone 10: ports the joint (gamma1, lambda1) recovery procedure from
% `..\parameter_recovery\synthetic_parameter_recovery_gp.m` (read-only
% reference, its jointRecoveryFit/jointRecoveryPlot) onto the griddatan
% fast predictor: coarse grid search over a (gamma1, lambda1) loss
% surface, then an `fminsearch` polish in log-space.
%
% Expected outcome, stated up front rather than discovered as a surprise:
% Milestone 9 found griddatan's nearest-neighbor search is functionally
% blind to gamma1 because delta/epsilon's raw scale (~10-100) dominates
% Euclidean distance. Lambda1's raw scale (~0.01-0.1) is just as dwarfed,
% so the same mechanism predicts the 2-D (gamma1, lambda1) loss surface
% will also be a near-flat plateau set almost entirely by delta/epsilon,
% not by gamma1/lambda1. This milestone checks whether that prediction
% actually holds, rather than assuming it.
%
% Deliberate deviation from the GP script: `synthetic_parameter_recovery_gp.m`
% only runs joint recovery on the first `size(dataset,1)/10` rows (a
% performance shortcut for the GP, where each `predict()` call is
% comparatively expensive). The griddatan predictor is ~5 orders of
% magnitude faster per query (Milestone 5), so this project runs joint
% recovery over the FULL holdout set instead -- cheap here, and this
% project's established priority is full-dataset honesty over cherry-picked
% subsets.

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

%% ---- Grid setup (identical ranges/counts to the GP script) ----
gamma1_range = [5e-5, 5e-2];
lambda_range = [0.01, 0.1];
logdist_fun = @(x,n) 10.^(log10(x(1)) + (log10(x(2)/x(1)) * rand(n,1)));
gridN = 50;
rng(42);
gamma_grid = sort(logdist_fun(gamma1_range, gridN));
lambda_grid = sort(logdist_fun(lambda_range, gridN));
[GG, LL] = ndgrid(gamma_grid, lambda_grid);

%% ---- Joint recovery over the full holdout set ----
n_rows = size(dataset, 1);
[refined_gamma1_joint, refined_lambda_joint, refined_loss_joint] = deal(zeros(n_rows, 1));

tic;
for row = 1:n_rows
    d_row = dataset(row,3); e_row = dataset(row,4);
    n1_row = dataset(row,5); n2_row = dataset(row,6);
    flux_scaled = (dataset(row,7) - flux_scaling.target_mean) / flux_scaling.target_std;
    [refined_gamma1_joint(row), refined_lambda_joint(row), refined_loss_joint(row)] = ...
        jointRecoveryFit(gamma_grid, lambda_grid, GG, LL, gamma1_range, lambda_range, ...
            model, flux_scaling, flux_scaled, d_row, e_row, n1_row, n2_row);
end
wall_time = toc;
fprintf('Joint recovery over all %d holdout rows: %.2f s (%.4f s/row)\n', n_rows, wall_time, wall_time/n_rows);

%% ---- Feasibility metric: log_error_joint across the FULL dataset ----
log_error_joint = sqrt((log10(refined_gamma1_joint) - log10(dataset(:,1))).^2 + ...
    (log10(refined_lambda_joint) - log10(dataset(:,2))).^2);
fprintf('\nlog_error_joint over all %d rows: mean=%.4f, median=%.4f, std=%.4f\n', ...
    n_rows, mean(log_error_joint), median(log_error_joint), std(log_error_joint));

% Did the grid search actually land ON the grid's own boundary for most
% rows -- a telltale sign the loss surface is flat and the "minimum" is
% arbitrary (whatever corner the grid or fminsearch happened to settle
% on), not a real optimum. Checked directly rather than only inferred from
% the error statistic.
at_gamma_boundary = abs(refined_gamma1_joint - gamma1_range(1)) < 1e-6*gamma1_range(1) | ...
    abs(refined_gamma1_joint - gamma1_range(2)) < 1e-6*gamma1_range(2);
at_lambda_boundary = abs(refined_lambda_joint - lambda_range(1)) < 1e-6*lambda_range(1) | ...
    abs(refined_lambda_joint - lambda_range(2)) < 1e-6*lambda_range(2);
fprintf('Rows where recovered gamma1 landed on its range boundary: %d/%d (%.0f%%)\n', ...
    sum(at_gamma_boundary), n_rows, 100*mean(at_gamma_boundary));
fprintf('Rows where recovered lambda1 landed on its range boundary: %d/%d (%.0f%%)\n', ...
    sum(at_lambda_boundary), n_rows, 100*mean(at_lambda_boundary));

%% ---- Save full-dataset recovery results ----
folder_results = fullfile(here, 'results');
if ~isfolder(folder_results); mkdir(folder_results); end
T_recovery = table((1:n_rows)', dataset(:,1), dataset(:,2), refined_gamma1_joint, refined_lambda_joint, ...
    refined_loss_joint, log_error_joint, ...
    'VariableNames', {'row','true_gamma1','true_lambda1','recovered_gamma1','recovered_lambda1','min_loss','log_error_joint'});
writetable(T_recovery, fullfile(folder_results, 'joint_recovery_full.csv'));
fprintf('\nFull per-row joint recovery results saved to results/joint_recovery_full.csv\n');

%% ---- Illustrative plots: best 10 rows by log_error_joint ----
[~, best_joint_rows] = mink(log_error_joint, 10);
outdir = fullfile(here, 'figures', 'joint_recovery');
if ~isfolder(outdir); mkdir(outdir); end
for row = best_joint_rows'
    g_row = dataset(row,1); l_row = dataset(row,2);
    d_row = dataset(row,3); e_row = dataset(row,4);
    n1_row = dataset(row,5); n2_row = dataset(row,6);
    flux_scaled = (dataset(row,7) - flux_scaling.target_mean) / flux_scaling.target_std;
    jointRecoveryPlot(gamma_grid, lambda_grid, GG, LL, model, flux_scaling, flux_scaled, ...
        d_row, e_row, n1_row, n2_row, row, g_row, l_row, ...
        refined_gamma1_joint(row), refined_lambda_joint(row), outdir);
end
fprintf('Saved 10 illustrative joint-recovery plots to %s\n', outdir);

fprintf('\nMilestone 10 complete.\n');

%% ---- Local functions ----

function model = fit_nearest_predictor(X_train, y_train)
    if size(unique(X_train, 'rows'), 1) < size(X_train, 1)
        warning('fit_nearest_predictor:DuplicateRows', ...
            'Training data contains duplicate feature rows; griddatan averages their targets, this k-d tree does not.');
    end
    model.ns = createns(X_train, 'NSMethod', 'kdtree');
    model.y = y_train;
end

function flux_z = griddatan_flux_scaled(model, gamma1, lambda1, delta, epsilon, N1, N2, flux_scaling)
    % Same broadcasting convention as s06's version: gamma1 and/or lambda1
    % may be arrays (here, the full 50x50 ndgrid); the other four
    % arguments are scalar per call.
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

function [gamma_hat, lambda_hat, loss_hat] = jointRecoveryFit(gamma_grid, lambda_grid, GG, LL, ...
        gamma1_range, lambda_range, model, flux_scaling, flux_scaled, d_row, e_row, n1_row, n2_row)
    % Same two-stage approach as the GP script: coarse grid search for a
    % starting point, then an unconstrained fminsearch polish in
    % log-space, clamped back into range afterward.
    loss_surface = (griddatan_flux_scaled(model, GG, LL, d_row, e_row, n1_row, n2_row, flux_scaling) - flux_scaled).^2;
    [~, lin_idx] = min(loss_surface(:));
    [gi, li] = ind2sub(size(loss_surface), lin_idx);

    polish_fun = @(logp) (griddatan_flux_scaled(model, 10^logp(1), 10^logp(2), d_row, e_row, n1_row, n2_row, flux_scaling) - flux_scaled)^2;
    [logp_hat, loss_hat] = fminsearch(polish_fun, [log10(gamma_grid(gi)), log10(lambda_grid(li))]);
    logp_hat(1) = min(max(logp_hat(1), log10(gamma1_range(1))), log10(gamma1_range(2)));
    logp_hat(2) = min(max(logp_hat(2), log10(lambda_range(1))), log10(lambda_range(2)));
    gamma_hat = 10^logp_hat(1);
    lambda_hat = 10^logp_hat(2);
end

function jointRecoveryPlot(gamma_grid, lambda_grid, GG, LL, model, flux_scaling, flux_scaled, ...
        d_row, e_row, n1_row, n2_row, row, g_row, l_row, gamma_hat, lambda_hat, outdir)
    % Same visual layout as the GP side's jointRecoveryPlot: log10 loss
    % heatmap, mean/±1sigma contours, recovered point, true point, axes
    % expanded (not clipped) to always include the true and recovered
    % points -- matching the GP script's CURRENT live behavior (it also
    % has an older commented-out "skip row if truth is outside the grid"
    % variant that is not active there, so this is not active here either).
    loss_surface = (griddatan_flux_scaled(model, GG, LL, d_row, e_row, n1_row, n2_row, flux_scaling) - flux_scaled).^2;
    figure('Visible', 'off');
    imagesc(log10(lambda_grid), log10(gamma_grid), log10(loss_surface));
    set(gca, 'YDir', 'normal');
    xlabel('log_{10} \lambda_1', 'FontSize', 15); ylabel('log_{10} \gamma_1', 'FontSize', 15);
    colorbar; title(sprintf('griddatan joint loss surface, row %d', row));
    ylabel(colorbar, 'log_{10} squared error', 'FontSize', 15);
    log_loss = log10(loss_surface(:));
    mean_log_loss = mean(log_loss);
    std_log_loss = std(log_loss);
    xlim([min([log10(lambda_grid(1)), log10(l_row), log10(lambda_hat)]), ...
          max([log10(lambda_grid(end)), log10(l_row), log10(lambda_hat)])]);
    ylim([min([log10(gamma_grid(1)), log10(g_row), log10(gamma_hat)]), ...
          max([log10(gamma_grid(end)), log10(g_row), log10(gamma_hat)])]);
    hold on;
    % A relative tolerance, not std_log_loss>0: MATLAB's contour() treats
    % data as "constant" (and silently skips drawing, with a warning) under
    % its own near-equality tolerance, which a std that is merely tiny
    % floating-point jitter (not exactly 0) can still trigger -- found by
    % inspecting the first render of this figure, which produced "Contour
    % not rendered for constant ZData" warnings AND a legend with
    % dashed/dotted entries pointing at lines that were never drawn.
    is_flat = std_log_loss < 1e-6 * max(abs(mean_log_loss), 1);
    if ~is_flat
        contour(log10(lambda_grid), log10(gamma_grid), log10(loss_surface), ...
            [mean_log_loss, mean_log_loss], 'Color','k', 'LineWidth',2, 'LineStyle','--', ...
            'DisplayName','mean log-loss');
        contour(log10(lambda_grid), log10(gamma_grid), log10(loss_surface), ...
            [mean_log_loss-std_log_loss, mean_log_loss-std_log_loss], 'Color','k', 'LineWidth',1, 'LineStyle',':', ...
            'DisplayName','mean \pm 1\sigma log-loss');
        contour(log10(lambda_grid), log10(gamma_grid), log10(loss_surface), ...
            [mean_log_loss+std_log_loss, mean_log_loss+std_log_loss], 'Color','w', 'LineWidth',1, 'LineStyle',':', ...
            'HandleVisibility','off');
    else
        % Surface is flat to within floating-point noise -- contour() has
        % nothing meaningful to draw. Stated on the plot and skipped
        % entirely, rather than left to silently fail with a warning and a
        % legend entry for a line that was never drawn.
        text(0.02, 0.02, 'Loss surface is flat (no gamma1/lambda1 sensitivity): contours omitted', ...
            'Units', 'normalized', 'FontSize', 9, 'Color', [1 1 1], 'VerticalAlignment', 'bottom');
    end
    plot(log10(lambda_hat), log10(gamma_hat), 'o', 'MarkerFaceColor','r', ...
        'MarkerEdgeColor','k', 'MarkerSize', 18, 'LineWidth', 2, 'DisplayName', 'recovered');
    plot(log10(l_row), log10(g_row), 'x', 'Color','m', ...
        'MarkerSize', 18, 'LineWidth', 3, 'DisplayName', 'true');
    legend('Location', 'northeast', 'FontSize', 12);
    hold off;
    filename = fullfile(outdir, sprintf('joint_recovery_row%04d.png', row));
    exportgraphics(gcf, filename, 'Resolution', 150);
    close(gcf);
end

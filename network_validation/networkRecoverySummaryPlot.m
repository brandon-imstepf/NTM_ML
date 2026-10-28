function T = networkRecoverySummaryPlot(rec_dir)
% One-figure summary of the single-parameter recoveries on the network, for
% both lambda2 conventions (rows: l2f = lambda2 0.01, l2eq = lambda2 = lambda1).
%   left:  models x parameters, median recovery error over the truths
%   right: models x parameters, share of candidate runs that failed (non-finite N)
%
% Error = |recovered - true| in search coordinates (log10 for gamma1 and
% lambda1, linear for delta and epsilon), as a fraction of the training range
% in those coordinates (the same measure as parameter_recovery/recoverySummaryPlot.m).
% A "b" after the number marks cells where at least half of the recovered values
% lie on a bound of the search range.
%
% Reads <rec_dir>/analysis_recovery_single_<conv>/recovery_network.mat
% (analyze_network_recovery.m). Writes <rec_dir>/recovery_summary_network.png
% and recovery_summary_network.csv (one row per convention, model and parameter).
%
%   networkRecoverySummaryPlot('outputs_recovery')

    convs = ["l2f", "l2eq"];
    conv_labels = ["\lambda_2 = 0.01", "\lambda_2 = \lambda_1"];
    param_range = [5e-5 5e-2; 0.01 0.1; 10 100; 10 100];
    convs_found = arrayfun(@(c) isfile(fullfile(rec_dir, "analysis_recovery_single_" + c, 'recovery_network.mat')), convs);
    convs = convs(convs_found);  conv_labels = conv_labels(convs_found);
    c_max = 0.3; % colour scale top for the error (errors above it share the darkest colour)

    T = table();
    fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 34 9*numel(convs) + 2]);
    tl = tiledlayout(fig, numel(convs), 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, 'Network single-parameter recovery (8 interior truths)', 'FontSize', 14);
    for c = 1:numel(convs)
        S = load(fullfile(rec_dir, "analysis_recovery_single_" + convs(c), 'recovery_network.mat'));
        rn = S.rn;
        n_models = numel(rn.model_labels);
        n_params = numel(rn.param_names);
        n_truths = numel(rn.truth_ids);

        [med, worst, bad_share, depth, best_nmse] = deal(nan(n_models, n_params));
        bound = false(n_models, n_params);
        for m = 1:n_models
            for i = 1:n_params
                u = @(x) x;
                if rn.param_log(i), u = @(x) log10(x); end
                lo = u(param_range(i,1));  hi = u(param_range(i,2));
                err = abs(u(squeeze(rn.hat(i,m,:))) - u(rn.true_params(:,i))) / (hi - lo);
                on_bound = abs(u(squeeze(rn.hat(i,m,:))) - lo) < 1e-3*(hi - lo) | ...
                    abs(u(squeeze(rn.hat(i,m,:))) - hi) < 1e-3*(hi - lo);
                loss = cell2mat(reshape(rn.loss(i,m,:), [], 1));
                med(m,i) = median(err, 'omitnan');
                worst(m,i) = max(err);
                bound(m,i) = mean(on_bound) >= 0.5;
                bad_share(m,i) = mean(~isfinite(loss));
                depth(m,i) = median(arrayfun(@(t) max(rn.loss{i,m,t}) / min(rn.loss{i,m,t}), 1:n_truths));
                best_nmse(m,i) = median(squeeze(rn.min_loss(i,m,:)), 'omitnan');
                T = [T; table(convs(c), rn.model_labels(m), rn.param_names(i), med(m,i), worst(m,i), ...
                    sum(on_bound), bad_share(m,i), depth(m,i), best_nmse(m,i), n_truths, 'VariableNames', ...
                    {'convention', 'model', 'param', 'median_error', 'worst_error', 'n_on_bound', ...
                    'share_failed', 'median_depth_max_over_min', 'median_best_nmse', 'n_truths'})]; %#ok<AGROW>
            end
        end

        % median error
        ax = nexttile(tl, 2*(c-1) + 1);
        imagesc(ax, min(med, c_max));
        set(ax, 'CLim', [0 c_max], 'TickLength', [0 0], 'XTick', 1:n_params, 'XTickLabel', rn.param_labels, ...
            'YTick', 1:n_models, 'YTickLabel', rn.model_labels, 'FontSize', 11);
        colormap(ax, sequentialMap([0.60 0.16 0.05]));
        cellText(ax, med, c_max, bound);
        title(ax, sprintf('Median error (fraction of range), %s', conv_labels(c)), 'FontWeight', 'normal');
        cb = colorbar(ax);
        ylabel(cb, sprintf('\\geq %.1f same colour', c_max));

        % failed or far runs
        ax = nexttile(tl, 2*(c-1) + 2);
        imagesc(ax, bad_share);
        set(ax, 'CLim', [0 0.25], 'TickLength', [0 0], 'XTick', 1:n_params, 'XTickLabel', rn.param_labels, ...
            'YTick', 1:n_models, 'YTickLabel', rn.model_labels, 'FontSize', 11);
        colormap(ax, sequentialMap([0.25 0.25 0.25]));
        cellText(ax, bad_share, 0.25, false(n_models, n_params));
        title(ax, sprintf('Share of runs failed (non-finite), %s', conv_labels(c)), 'FontWeight', 'normal');
        colorbar(ax);
    end
    exportgraphics(fig, fullfile(rec_dir, 'recovery_summary_network.png'), 'Resolution', 200);
    close(fig);
    writetable(T, fullfile(rec_dir, 'recovery_summary_network.csv'));
end

function cellText(ax, values, c_top, bound)
% value of every cell written on it; white ink on dark cells
    for m = 1:size(values, 1)
        for i = 1:size(values, 2)
            txt = sprintf('%.2f', values(m,i));
            if bound(m,i), txt = [txt ' b']; end
            ink = [0 0 0];
            if min(values(m,i), c_top) > 0.6*c_top, ink = [1 1 1]; end
            text(ax, i, m, txt, 'HorizontalAlignment', 'center', 'Color', ink, 'FontSize', 10);
        end
    end
end

function cmap = sequentialMap(dark)
% one hue, light (small) to dark (large)
    light = [0.97 0.97 0.97];
    t = linspace(0, 1, 256)';
    cmap = light + t .* (dark - light);
end

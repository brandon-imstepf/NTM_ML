function T = recoverySummaryPlot(outdir)
% One-figure summary of the single-parameter recoveries, for each output the
% loss was computed on separately (rows: flux, w1) and every group of locked sets
% (columns: interior, edge, ood). Each panel is models (rows) x parameters
% (columns); a cell is the median recovery error over the sets of the group.
%
% Error = |recovered - true| in search coordinates (log10 for gamma1 and
% lambda1, linear for delta and epsilon), as a fraction of the training range
% in those coordinates. 0 = exact, 1 = the whole range. For out-of-range sets
% the true value lies outside the search range, so the error is measured to the
% nearest bound of the range (the closest value the search can return).
% Cells where the parameter is flat (not identifiable) in every set of the
% group are grey and marked "n.i.". A "b" after the number marks cells where
% most recovered values lie on a bound of the search range: for edge and
% out-of-range sets that is where the truth (or its nearest value) is, so a
% small error there can come from the search stopping at the bound.
%
% Reads <outdir>/<output>/recovery_curves.mat (synthetic_parameter_recovery_all_models.m).
% Writes <outdir>/recovery_summary_single.png and recovery_summary_single.csv.
%
%   recoverySummaryPlot('recovery_all_models')

    outputs = ["flux", "w1"];   % each output on its own (recovery_all_models/both/ is not shown)
    groups = ["interior", "edge", "ood"];
    group_labels = ["Interior", "Edge", "Out-of-range"];
    outputs = outputs(arrayfun(@(o) isfile(fullfile(outdir, o, 'recovery_curves.mat')), outputs));
    c_max = 0.3; % colour scale top (errors above it share the darkest colour)

    T = table();
    fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 38 9.5*numel(outputs) + 2]);
    tl = tiledlayout(fig, numel(outputs), numel(groups), 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, 'Single-parameter recovery: median error (fraction of training range)', 'FontSize', 14);
    for o = 1:numel(outputs)
        S = load(fullfile(outdir, outputs(o), 'recovery_curves.mat'));
        rc = S.rc;
        n_models = numel(rc.model_labels);
        n_params = numel(rc.params);
        set_group = string(regexp(rc.set_ids, '^[^_]+', 'match', 'once'));

        % error of every recovery, n_params x n_models x n_sets
        err = nan(n_params, n_models, numel(rc.set_ids));
        on_bound = false(n_params, n_models, numel(rc.set_ids));
        for i = 1:n_params
            p = rc.params(i);
            u = @(x) x;
            if p.log, u = @(x) log10(x); end
            lo = u(p.range(1));  hi = u(p.range(2));
            u_true = min(max(u(rc.true_params(:,i)), lo), hi);   % truth clamped into the range
            u_hat = u(squeeze(rc.single_hat(i,:,:)));             % n_models x n_sets
            err(i,:,:) = abs(u_hat - u_true') / (hi - lo);
            on_bound(i,:,:) = abs(u_hat - lo) < 1e-3*(hi - lo) | abs(u_hat - hi) < 1e-3*(hi - lo);
        end

        for g = 1:numel(groups)
            in_group = set_group == groups(g);
            med = median(err(:,:,in_group), 3, 'omitnan')';                 % n_models x n_params
            flat = all(rc.single_flat(:,:,in_group), 3)';
            bound = (mean(on_bound(:,:,in_group), 3) > 0.5)';
            ax = nexttile(tl, (o-1)*numel(groups) + g);
            shown = min(med, c_max);
            shown(flat) = NaN;
            imagesc(ax, shown, 'AlphaData', ~isnan(shown));
            set(ax, 'Color', [0.85 0.85 0.85], 'CLim', [0 c_max], 'TickLength', [0 0], ...
                'XTick', 1:n_params, 'XTickLabel', {rc.params.label}, ...
                'YTick', 1:n_models, 'YTickLabel', rc.model_labels, 'FontSize', 11);
            colormap(ax, sequentialMap());
            for m = 1:n_models
                for i = 1:n_params
                    if flat(m,i)
                        txt = 'n.i.';  ink = [0.3 0.3 0.3];
                    else
                        txt = sprintf('%.2f', med(m,i));
                        if bound(m,i), txt = [txt ' b']; end
                        ink = [0 0 0];
                        if shown(m,i) > 0.6*c_max, ink = [1 1 1]; end
                    end
                    text(ax, i, m, txt, 'HorizontalAlignment', 'center', 'Color', ink, 'FontSize', 10);
                end
            end
            title(ax, sprintf('%s, loss on %s (%d sets)', group_labels(g), rc.output_label, sum(in_group)), ...
                'FontWeight', 'normal');
            for m = 1:n_models
                for i = 1:n_params
                    T = [T; table(outputs(o), groups(g), rc.model_labels(m), string(rc.params(i).name), ...
                        med(m,i), flat(m,i), bound(m,i), sum(in_group), 'VariableNames', ...
                        {'output', 'group', 'model', 'param', 'median_error', 'flat', 'on_bound', 'n_sets'})]; %#ok<AGROW>
                end
            end
        end
    end
    cb = colorbar(ax);
    cb.Layout.Tile = 'east';
    ylabel(cb, sprintf('median |error| / range  (\\geq %.1f same colour)', c_max));
    exportgraphics(fig, fullfile(outdir, 'recovery_summary_single.png'), 'Resolution', 200);
    close(fig);
    writetable(T, fullfile(outdir, 'recovery_summary_single.csv'));
end

function cmap = sequentialMap()
% one hue, light (small error) to dark (large error)
    light = [0.97 0.97 0.97];
    dark = [0.60 0.16 0.05];
    t = linspace(0, 1, 256)';
    cmap = light + t .* (dark - light);
end

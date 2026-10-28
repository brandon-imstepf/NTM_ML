function recoveryGroupMeanPlot(rc, group, outdir)
% Mean recovery figures for one group of locked sets (e.g. "interior", "edge",
% "ood"), every model side by side (columns). A set belongs to the group named
% by the start of its set_id, up to the first underscore.
%
% The sets in a group have different true values, so the curves are first
% moved onto axes relative to the truth, then averaged:
%     log-searched p (gamma1, lambda1):     log10(p / true p)
%     linear-searched p (delta, epsilon):   p - true p
% 0 is the truth, so a model that recovers well has its minimum at 0.
%   by_group/single_mean_<group>.png  one row per parameter: each set's loss curve
%                                     thin, the group mean thick, recovered dots
%   by_group/joint_mean_<group>.png   one row per pair: mean log10 loss surface,
%                                     each set's recovered point (dots), truth (x)
% The search grid covers the training range only, so a set does not cover the
% whole relative axis (an edge set whose truth is on the lower bound has nothing
% below 0). The mean is drawn only where every set in the group has a value;
% elsewhere only the thin per-set curves or grey show. Flat (not identifiable)
% recoveries have no dots; a pair flat along both parameters is left blank.
% rc is the struct saved in recovery_curves.mat.

    set_group = string(regexp(rc.set_ids, '^[^_]+', 'match', 'once'));
    in_group = find(set_group == group);
    if numel(in_group) < 2
        return % a mean of one set is the set's own figure
    end
    n_models = numel(rc.model_labels);
    n_sets = numel(in_group);
    n_params = numel(rc.params);
    n_pairs = size(rc.pairs, 1);
    true_params = rc.true_params(in_group,:);
    group_title = groupLabel(group);
    if isfield(rc, 'output_label')
        group_title = sprintf('%s  [loss on %s]', group_title, rc.output_label);
    end
    group_outdir = fullfile(outdir, 'by_group');
    if ~exist(group_outdir, 'dir'), mkdir(group_outdir); end

    % relative axis for each parameter, as wide as its search range on either side of 0
    rel = cell(n_params, 1);
    for i = 1:n_params
        p = rc.params(i);
        if p.log
            rel{i}.fun = @(x, x_true) log10(x ./ x_true);
            rel{i}.label = sprintf('log_{10}(%s / true)', p.label);
            span = log10(p.range(2) / p.range(1));
        else
            rel{i}.fun = @(x, x_true) x - x_true;
            rel{i}.label = sprintf('%s - true', p.label);
            span = diff(p.range);
        end
        rel{i}.curve_axis = linspace(-span, span, 241)';
        rel{i}.grid_axis = linspace(-span, span, 61)';
    end

    %% single-parameter figure
    fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 40 8.5*n_params]);
    tl = tiledlayout(fig, n_params, n_models, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, sprintf('%s: mean over %d sets (single-parameter)', group_title, n_sets), 'FontSize', 14);
    for i = 1:n_params
        p = rc.params(i);
        x_axis = rel{i}.curve_axis;
        for m = 1:n_models
            ax = nexttile(tl, (i-1)*n_models + m); hold(ax, 'on');
            rgb = hex2rgb(rc.model_colors(m));
            curves = nan(numel(x_axis), n_sets);
            for k = 1:n_sets
                s = in_group(k);
                curves(:,k) = interp1(rel{i}.fun(p.curve_grid, true_params(k,i)), ...
                    log10(rc.single_curves(:,i,m,s)), x_axis, 'linear', NaN);
                plot(ax, x_axis, curves(:,k), 'Color', [rgb 0.3], 'LineWidth', 0.75);
            end
            % no 'omitnan': a point any set does not cover stays NaN (not drawn)
            plot(ax, x_axis, mean(curves, 2), 'Color', rgb, 'LineWidth', 2.5);
            xline(ax, 0, '--k', 'LineWidth', 1);
            hats = squeeze(rc.single_hat(i,m,in_group));
            losses = squeeze(rc.single_loss(i,m,in_group));
            plot(ax, rel{i}.fun(hats, true_params(:,i)), log10(losses), 'o', ...
                'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'MarkerSize', 5);   % NaN (flat) not drawn
            if all(rc.single_flat(i,m,in_group))
                title(ax, sprintf('%s:  %s flat (not identifiable)', rc.model_labels(m), p.label), 'FontWeight', 'normal');
            else
                title(ax, sprintf('%s:  median |error| = %.2g', rc.model_labels(m), ...
                    median(abs(rel{i}.fun(hats, true_params(:,i))), 'omitnan')), 'FontWeight', 'normal');
            end
            set(ax, 'Box', 'off', 'TickDir', 'out');
            xlim(ax, x_axis([1 end]));
            grid(ax, 'on');
            xlabel(ax, rel{i}.label);
            if m == 1, ylabel(ax, 'log_{10} loss'); end
            hold(ax, 'off');
        end
    end
    ax1 = nexttile(tl, 1); hold(ax1, 'on');
    h = [plot(ax1, NaN, NaN, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.75), ...
         plot(ax1, NaN, NaN, '-', 'Color', [0.3 0.3 0.3], 'LineWidth', 2.5), ...
         plot(ax1, NaN, NaN, 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k'), ...
         plot(ax1, NaN, NaN, '--k')];
    lg = legend(h, {'each set', 'group mean', 'recovered, each set', 'truth'}, ...
        'Orientation', 'horizontal', 'Box', 'off', 'FontSize', 11);
    lg.Layout.Tile = 'south';
    hold(ax1, 'off');
    exportgraphics(fig, fullfile(group_outdir, sprintf('single_mean_%s.png', group)), 'Resolution', 130);
    close(fig);

    %% joint figure
    % mean surfaces on the relative axes, first, for one shared colour scale
    mean_surfaces = cell(n_pairs, n_models);
    for k = 1:n_pairs
        a = rc.pairs(k,1);  b = rc.pairs(k,2);
        [DA, DB] = ndgrid(rel{a}.grid_axis, rel{b}.grid_axis);
        for m = 1:n_models
            surfaces = nan([size(DA), n_sets]);
            for j = 1:n_sets
                s = in_group(j);
                F = griddedInterpolant({rel{a}.fun(rc.params(a).joint_grid, true_params(j,a)), ...
                    rel{b}.fun(rc.params(b).joint_grid, true_params(j,b))}, ...
                    log10(rc.joint_surfaces(:,:,k,m,s)), 'linear', 'none');
                surfaces(:,:,j) = F(DA, DB);
            end
            mean_surfaces{k,m} = mean(surfaces, 3);
        end
    end
    all_means = cell2mat(cellfun(@(x) x(:), mean_surfaces(:), 'UniformOutput', false));
    clim_all = [min(all_means), max(all_means)];
    if ~all(isfinite(clim_all)) || clim_all(1) == clim_all(2), clim_all = [-1 1]; end
    % contour levels, the same in every panel so a shade means the same loss everywhere
    n_levels = 20;
    levels = linspace(clim_all(1), clim_all(2), n_levels + 1);

    fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 40 8.5*n_pairs]);
    tl = tiledlayout(fig, n_pairs, n_models, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, sprintf('%s: mean over %d sets (joint)', group_title, n_sets), 'FontSize', 14);
    for k = 1:n_pairs
        a = rc.pairs(k,1);  b = rc.pairs(k,2);   % a on the y axis, b on the x axis
        for m = 1:n_models
            ax = nexttile(tl, (k-1)*n_models + m); hold(ax, 'on');
            flat = squeeze(rc.joint_flat(k,:,m,in_group));   % 2 x n_sets
            if all(flat, 'all')
                axis(ax, 'off');
                title(ax, sprintf('%s:  %s, %s flat (not identifiable)', rc.model_labels(m), ...
                    rc.params(a).label, rc.params(b).label), 'FontWeight', 'normal');
                hold(ax, 'off');
                continue
            end
            % filled contours of mean log10 loss (shading); NaN (not covered by every set) stays grey
            [X, Y] = meshgrid(rel{b}.grid_axis, rel{a}.grid_axis);
            contourf(ax, X, Y, mean_surfaces{k,m}, levels, 'LineColor', [0.25 0.25 0.25], 'LineWidth', 0.3);
            img_ax = ax; % any surface panel carries the shared colour scale
            set(ax, 'Box', 'off', 'TickDir', 'out', 'CLim', clim_all, 'Color', [0.85 0.85 0.85], 'Layer', 'top');
            axis(ax, 'tight');
            hats = squeeze(rc.joint_hat(k,:,m,in_group));   % 2 x n_sets, NaN where flat
            plot(ax, rel{b}.fun(hats(2,:)', true_params(:,b)), rel{a}.fun(hats(1,:)', true_params(:,a)), 'o', ...
                'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'MarkerSize', 6);
            plot(ax, 0, 0, 'x', 'Color', 'm', 'MarkerSize', 14, 'LineWidth', 2.5);
            if all(flat(1,:))
                title(ax, sprintf('%s:  %s flat', rc.model_labels(m), rc.params(a).label), 'FontWeight', 'normal');
            elseif all(flat(2,:))
                title(ax, sprintf('%s:  %s flat', rc.model_labels(m), rc.params(b).label), 'FontWeight', 'normal');
            else
                title(ax, rc.model_labels(m), 'FontWeight', 'normal');
            end
            xlabel(ax, rel{b}.label);
            if m == 1, ylabel(ax, rel{a}.label); end
            hold(ax, 'off');
        end
    end
    cb = colorbar(img_ax);
    cb.Layout.Tile = 'east';
    ylabel(cb, 'mean log_{10} loss');
    ax1 = nexttile(tl, 1); hold(ax1, 'on');
    h = [plot(ax1, NaN, NaN, 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k'), ...
         plot(ax1, NaN, NaN, 'x', 'Color', 'm', 'LineWidth', 2)];
    lg = legend(h, {'recovered, each set', 'truth'}, 'Orientation', 'horizontal', 'Box', 'off', 'FontSize', 11);
    lg.Layout.Tile = 'south';
    hold(ax1, 'off');
    exportgraphics(fig, fullfile(group_outdir, sprintf('joint_mean_%s.png', group)), 'Resolution', 130);
    close(fig);
end

function label = groupLabel(group)
% name shown in figures for a set_id prefix (the same names as the network-validation table)
    switch group
        case "interior", label = 'Interior';
        case "edge",     label = 'Edge';
        case "ood",      label = 'Out-of-range';
        otherwise,       label = char(group);
    end
end

function rgb = hex2rgb(hex)
    h = char(hex);
    rgb = hex2dec({h(2:3), h(4:5), h(6:7)})' / 255;
end

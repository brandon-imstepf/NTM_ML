function plotNelderMeadRecovery(nm_dir, outdir, label, extra, truth_file)
% Figures for the Nelder-Mead network recoveries in nm_dir (the nm_*.mat files
% of nelder_mead_network_recovery.m, one per parameter, all for one truth):
%   nm_search_<label>.png       loss vs parameter value: every simulation of the
%                               search, in the order fminsearch ran them
%   nm_convergence_<label>.png  parameter value and loss vs simulation number
%   nm_tau_<label>.png          NMSE(t) of node tau against the truth over time,
%                               for the start and the recovered simulation
%   nm_regions_<label>_<param>.png  tau over time in each region: truth, start
%                               and recovered, one figure per searched parameter
%   nm_region_loss_<label>.png  loss of each region alone, start vs recovered
% extra (optional) is a struct array with fields param, value, N: further
% surrogate runs at the truth except for one parameter, drawn in the search
% figure as a loss curve (the loss is computed here against the same truth).
% A parameter with extra runs but no search gets a panel with the curve only.
% truth_file (optional) is the truth run the searches were fitted to; by default
% the path saved in the search results (result.truth_file).
%
%   plotNelderMeadRecovery('outputs_nm', 'outputs_nm', 'mlp_vs_numeric')

    if nargin < 4 || isempty(extra), extra = struct('param', {}, 'value', {}, 'N', {}); end
    if ~isfolder(outdir), mkdir(outdir); end
    param_names = ["gamma1", "lambda1", "delta", "epsilon"];
    param_labels = ["\gamma_1", "\lambda_1", "\delta", "\epsilon"];
    param_range = [5e-5 5e-2; 0.01 0.1; 10 100; 10 100];
    param_log = [true, true, false, false];

    % searches, one per parameter
    runs = cell(1, numel(param_names));
    for f = dir(fullfile(nm_dir, 'nm_*.mat'))'
        S = load(fullfile(f.folder, f.name), 'trace', 'result');
        if ~isfield(S, 'result') || numel(S.result.param) > 1, continue, end   % not finished, or several parameters (plotNelderMeadRecovery2D)
        runs{param_names == S.result.param} = S;
    end
    has_run = ~cellfun(@isempty, runs);
    shown = find(has_run | ismember(param_names, string({extra.param})));
    assert(any(has_run), 'no finished searches in %s', nm_dir);

    % truth: true values and node tau
    r1 = runs{find(has_run, 1)}.result;
    here = fileparts(mfilename('fullpath'));
    T = readtable(fullfile(here, 'manifests', 'recovery_truths.csv'), 'TextType', 'string');
    T = T(T.convention == r1.conv & T.truth_id == r1.truth_id, :);
    theta_true = [T.gamma1, T.lambda1, T.delta, T.epsilon];
    if nargin < 5 || isempty(truth_file), truth_file = r1.truth_file; end
    truth = load(truth_file, 'N', 'trange', 'meta', 'region_names');
    days = truth.trange * 180;
    nmse = @(N) sum((N(:) - truth.N(:)).^2) / sum(truth.N(:).^2);
    nmse_t = @(N) sum((N - truth.N).^2, 1) ./ sum(truth.N.^2, 1);
    truth_name = upper(string(truth.meta.method));
    if truth_name == "NUMERIC", truth_name = "numeric"; end
    % two lines, so the title fits figures of one panel
    fig_title = {sprintf('Nelder-Mead recovery, %s vs %s truth, %g days', upper(r1.method), truth_name, round(days(end), 1)), ...
        sprintf('true values:  \\gamma_1 = %g,  \\lambda_1 = %g,  \\delta = %g,  \\epsilon = %g', theta_true)};
    n_max = max(cellfun(@(S) numel(S.trace), runs(has_run)));   % one colour scale for search order

    order_map = sequentialMap();
    path_grey = [0.6 0.6 0.6];

    %% 1. Search: loss vs parameter value
    n_rows = 1 + (numel(shown) > 3);   % four panels as 2 x 2, so the figure fits a page width
    n_cols = ceil(numel(shown) / n_rows);
    fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 max(14*n_cols, 24) 11*n_rows + 3]);
    tl = tiledlayout(fig, n_rows, n_cols, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, fig_title, 'FontSize', 13);
    legend_h = gobjects(0);  legend_names = {};
    for k = 1:numel(shown)
        i = shown(k);
        ax = nexttile(tl); hold(ax, 'on');
        h = gobjects(0);  names = {};
        ex = extra(string({extra.param}) == param_names(i));
        if ~isempty(ex)
            [v, o] = sort([ex.value]);
            L = arrayfun(@(e) nmse(e.N), ex(o));
            h(end+1) = plot(ax, v, L, '-s', 'Color', path_grey, 'MarkerFaceColor', path_grey, ...
                'MarkerSize', 5, 'LineWidth', 1.5);
            names{end+1} = 'loss at fixed values (not part of the search)';
        end
        if has_run(i)
            tr = runs{i}.trace;  res = runs{i}.result;
            v = [tr.value];  L = [tr.nmse];  n = numel(v);
            plot(ax, v, L, '-', 'Color', [path_grey 0.7], 'LineWidth', 0.75);   % order of the search
            h(end+1) = scatter(ax, v, L, 60, 1:n, 'filled', 'MarkerEdgeColor', 'k');
            names{end+1} = 'simulation of the search (colour = order)';
            text(ax, v(1), L(1), '  start', 'FontSize', 10, 'VerticalAlignment', 'top');
            [~, b] = min(L);
            h(end+1) = plot(ax, v(b), L(b), 'o', 'MarkerSize', 15, 'Color', [0.8 0.1 0.1], 'LineWidth', 2);
            names{end+1} = 'recovered (lowest loss)';
            colormap(ax, order_map);  clim(ax, [1 n_max]);
            order_ax = ax;
            title(ax, sprintf('%s: true %.3g, recovered %.3g  (%d simulations)', ...
                param_labels(i), res.true_value, res.recovered_value, n), 'FontWeight', 'normal');
        else
            title(ax, sprintf('%s: true %.3g  (no search run)', param_labels(i), theta_true(i)), 'FontWeight', 'normal');
        end
        h(end+1) = xline(ax, theta_true(i), '--k', 'LineWidth', 1.5);
        names{end+1} = 'true value';
        if param_log(i), set(ax, 'XScale', 'log'); end
        set(ax, 'YScale', 'log', 'Box', 'off', 'TickDir', 'out', 'FontSize', 11);
        shown_values = theta_true(i);
        if has_run(i), shown_values = [shown_values, runs{i}.trace.value]; end
        xlim(ax, axisLimits([shown_values, ex.value], param_range(i,:), param_log(i)));
        grid(ax, 'on');
        xlabel(ax, param_labels(i));
        if mod(k - 1, n_cols) == 0, ylabel(ax, 'loss: NMSE of node \tau'); end
        % legend: every kind of mark used in any panel, once
        for j = find(~ismember(names, legend_names))
            legend_h(end+1) = h(j);  legend_names{end+1} = names{j}; %#ok<AGROW>
        end
        hold(ax, 'off');
    end
    lg = legend(legend_h, legend_names, 'Orientation', 'horizontal', 'Box', 'off', 'FontSize', 11);
    lg.Layout.Tile = 'south';
    cb = colorbar(order_ax);
    cb.Layout.Tile = 'east';
    ylabel(cb, 'simulation number (order of the search)');
    exportgraphics(fig, fullfile(outdir, sprintf('nm_search_%s.png', label)), 'Resolution', 150);
    close(fig);

    %% 2. Convergence: value and loss vs simulation number
    searched = find(has_run);
    fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 max(12*numel(searched), 20) 19]);
    tl = tiledlayout(fig, 2, numel(searched), 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, fig_title, 'FontSize', 13);
    for k = 1:numel(searched)
        i = searched(k);
        tr = runs{i}.trace;  n = numel(tr);
        ax = nexttile(tl, k); hold(ax, 'on');
        plot(ax, 1:n, [tr.value], '-o', 'Color', [0.15 0.3 0.6], 'MarkerFaceColor', [0.15 0.3 0.6], 'LineWidth', 1.5);
        yline(ax, theta_true(i), '--k', 'true', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left');
        if param_log(i), set(ax, 'YScale', 'log'); end
        ylim(ax, axisLimits([theta_true(i), tr.value], param_range(i,:), param_log(i)));
        set(ax, 'Box', 'off', 'TickDir', 'out', 'FontSize', 11);
        grid(ax, 'on');
        title(ax, param_labels(i), 'FontWeight', 'normal');
        ylabel(ax, sprintf('%s simulated', param_labels(i)));
        hold(ax, 'off');

        ax = nexttile(tl, numel(searched) + k); hold(ax, 'on');
        L = [tr.nmse];
        h = [plot(ax, 1:n, L, 'o', 'Color', [0.15 0.3 0.6], 'MarkerFaceColor', [0.15 0.3 0.6]), ...
             stairs(ax, 1:n, cummin(L), '-k', 'LineWidth', 1.5)];
        set(ax, 'YScale', 'log', 'Box', 'off', 'TickDir', 'out', 'FontSize', 11);
        grid(ax, 'on');
        xlabel(ax, 'simulation number');
        if k == 1, ylabel(ax, 'loss: NMSE of node \tau'); end
        hold(ax, 'off');
    end
    lg = legend(h, {'loss of that simulation', 'lowest loss so far'}, 'Orientation', 'horizontal', 'Box', 'off', 'FontSize', 11);
    lg.Layout.Tile = 'south';
    exportgraphics(fig, fullfile(outdir, sprintf('nm_convergence_%s.png', label)), 'Resolution', 150);
    close(fig);

    %% 3. Node tau error over time: start vs recovered
    keep = 2:numel(days);   % skip t = 0 (both start from the same seeding)
    fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 max(12*numel(searched), 20) 11]);
    tl = tiledlayout(fig, 1, numel(searched), 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, fig_title, 'FontSize', 13);
    for k = 1:numel(searched)
        i = searched(k);
        tr = runs{i}.trace;
        [~, b] = min([tr.nmse]);
        ax = nexttile(tl); hold(ax, 'on');
        e_start = nmse_t(tr(1).N);  e_best = nmse_t(tr(b).N);
        h = [plot(ax, days(keep), e_start(keep), '-', 'Color', path_grey, 'LineWidth', 2), ...
             plot(ax, days(keep), e_best(keep), '-', 'Color', [0.8 0.1 0.1], 'LineWidth', 2)];
        set(ax, 'YScale', 'log', 'Box', 'off', 'TickDir', 'out', 'FontSize', 11);
        grid(ax, 'on');
        title(ax, sprintf('%s:  start %.3g,  recovered %.3g', param_labels(i), tr(1).value, tr(b).value), 'FontWeight', 'normal');
        xlabel(ax, 'days');
        if k == 1, ylabel(ax, 'NMSE(t) of node \tau'); end
        hold(ax, 'off');
    end
    lg = legend(h, {'start of the search (middle of the range)', 'recovered'}, 'Orientation', 'horizontal', 'Box', 'off', 'FontSize', 11);
    lg.Layout.Tile = 'south';
    exportgraphics(fig, fullfile(outdir, sprintf('nm_tau_%s.png', label)), 'Resolution', 150);
    close(fig);

    %% 4. Tau per region: truth, start and recovered, one panel per region
    n_reg = size(truth.N, 1);
    reg_cols = ceil(sqrt(n_reg * 1.2));
    reg_rows = ceil(n_reg / reg_cols);
    reg_names = regionNames(truth, n_reg);
    nmse_r = @(N) sum((N - truth.N).^2, 2) ./ sum(truth.N.^2, 2);   % one loss per region
    for k = 1:numel(searched)
        i = searched(k);
        tr = runs{i}.trace;
        [~, b] = min([tr.nmse]);
        e_reg = nmse_r(tr(b).N);
        fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 6*reg_cols 5*reg_rows + 3]);
        tl = tiledlayout(fig, reg_rows, reg_cols, 'TileSpacing', 'compact', 'Padding', 'compact');
        title(tl, [fig_title, {sprintf('%s search: start %.3g, recovered %.3g', param_labels(i), tr(1).value, tr(b).value)}], 'FontSize', 13);
        for r = 1:n_reg
            ax = nexttile(tl); hold(ax, 'on');
            h = [plot(ax, days, truth.N(r,:), '-k', 'LineWidth', 2), ...
                 plot(ax, days, tr(1).N(r,:), ':', 'Color', path_grey, 'LineWidth', 1.75), ...
                 plot(ax, days, tr(b).N(r,:), '--', 'Color', [0.8 0.1 0.1], 'LineWidth', 1.75)];
            set(ax, 'Box', 'off', 'TickDir', 'out', 'FontSize', 8);
            grid(ax, 'on');
            title(ax, [wrapName(reg_names(r)), sprintf('loss %.2g', e_reg(r))], 'FontWeight', 'normal', 'FontSize', 9, 'Interpreter', 'none');
            if r > n_reg - reg_cols, xlabel(ax, 'days'); end
            if mod(r - 1, reg_cols) == 0, ylabel(ax, '\tau'); end
            hold(ax, 'off');
        end
        lg = legend(h, {sprintf('truth (%s)', truth_name), 'start of the search', 'recovered'}, ...
            'Orientation', 'horizontal', 'Box', 'off', 'FontSize', 11);
        lg.Layout.Tile = 'south';
        exportgraphics(fig, fullfile(outdir, sprintf('nm_regions_%s_%s.png', label, param_names(i))), 'Resolution', 150);
        close(fig);
    end

    %% 5. Loss per region: start vs recovered, one row per searched parameter
    fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 34 8*numel(searched) + 4]);
    tl = tiledlayout(fig, numel(searched), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, fig_title, 'FontSize', 13);
    for k = 1:numel(searched)
        i = searched(k);
        tr = runs{i}.trace;
        [~, b] = min([tr.nmse]);
        ax = nexttile(tl); hold(ax, 'on');
        h = [plot(ax, 1:n_reg, nmse_r(tr(1).N), 'o', 'Color', path_grey, 'MarkerFaceColor', path_grey, 'MarkerSize', 7), ...
             plot(ax, 1:n_reg, nmse_r(tr(b).N), 'o', 'Color', [0.8 0.1 0.1], 'MarkerFaceColor', [0.8 0.1 0.1], 'MarkerSize', 7), ...
             yline(ax, nmse(tr(b).N), '--', 'Color', [0.8 0.1 0.1], 'LineWidth', 1)];
        set(ax, 'YScale', 'log', 'Box', 'off', 'TickDir', 'out', 'FontSize', 10, ...
            'XTick', 1:n_reg, 'XTickLabel', [], 'TickLabelInterpreter', 'none');
        xlim(ax, [0.5 n_reg + 0.5]);
        grid(ax, 'on');
        title(ax, sprintf('%s:  start %.3g,  recovered %.3g', param_labels(i), tr(1).value, tr(b).value), 'FontWeight', 'normal');
        ylabel(ax, 'loss of the region');
        hold(ax, 'off');
    end
    set(ax, 'XTickLabel', reg_names, 'XTickLabelRotation', 60);   % region names on the bottom row only
    lg = legend(h, {'start of the search', 'recovered', 'recovered, all regions together'}, ...
        'Orientation', 'horizontal', 'Box', 'off', 'FontSize', 11);
    lg.Layout.Tile = 'south';
    exportgraphics(fig, fullfile(outdir, sprintf('nm_region_loss_%s.png', label)), 'Resolution', 150);
    close(fig);
end

function lim = axisLimits(values, range, is_log)
% the whole range, or around the values when they cover less than a tenth of it
% (searches started close to the truth)
    f = @(x) x;  g = @(u) u;
    if is_log, f = @log10;  g = @(u) 10.^u; end
    u = f(values);  span = diff(f(range));
    lim = range;
    if max(u) - min(u) < 0.1 * span
        pad = max(0.15 * (max(u) - min(u)), 1e-3 * span);
        lim = g([min(u) - pad, max(u) + pad]);
    end
end

function names = regionNames(truth, n_reg)
% region names saved with the run, or numbers if the run has none
    if isfield(truth, 'region_names')
        names = string(truth.region_names(:))';
    else
        names = "region " + string(1:n_reg);
    end
end

function cmap = sequentialMap()
% one hue, light (early simulations) to dark (late)
    light = [0.78 0.86 0.97];
    dark = [0.05 0.17 0.45];
    t = linspace(0, 1, 256)';
    cmap = light + t .* (dark - light);
end

function lines = wrapName(name)
% long region names on two lines, split at the first comma
    lines = name;
    k = strfind(name, ", ");
    if strlength(name) > 26 && ~isempty(k)
        lines = [extractBefore(name, k(1) + 1), extractAfter(name, k(1) + 1)];
    end
end

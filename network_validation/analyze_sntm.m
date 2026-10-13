function analyze_sntm(out_dir)
% Turn the aggregate_sntm summaries into the numbers and figures for the
% network-wide validation section. Run aggregate_sntm(out_dir) first.
%
% Writes to <out_dir>/summary/:
%   paper_tables.tex          rows for tab:validation-accuracy and tab:validation-coverage,
%                             plus runtime and non-finite counts
%   fig_validation_dist.pdf/png   per-simulation pooled R^2, pooled NMSE, median b0 by group
%   fig_validation_time.pdf/png   NMSE(t) and b0(t) (median across settings, IQR band),
%                                 and the numeric runs' regional CV(t)
% Pooled = all regions and time points of a simulation at once. R^2 at a single
% time point is not plotted: it is unstable once N is nearly equal across regions.

if nargin < 1, cfg = sntm_config(); out_dir = cfg.out_dir; end
sum_dir = fullfile(out_dir, 'summary');
C = readtable(fullfile(sum_dir, 'accuracy_by_case.csv'), 'TextType', 'string');
T = readtable(fullfile(sum_dir, 'accuracy_by_time.csv'), 'TextType', 'string');
V = readtable(fullfile(sum_dir, 'coverage.csv'), 'TextType', 'string');

methods = ["mlp", "sr", "gp", "interp", "lr"];
mlabels = ["MLP", "SR", "GP", "NN-int.", "LR"];
groups = ["interior", "edge", "ood", "interior_r0.002", "fig8", "fig8_r0.01"];
glabels = ["Interior", "Edge", "Out-of-range", "Interior (init 0.002)", "Fig. 8 (init 0.002)", "Fig. 8 (init 0.01)"];
colors = hex2rgb(["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]);
markers = {'o', 's', 'd', '^', 'v'};
groups = groups(ismember(groups, unique(C.group)));
glabels = glabels(ismember(["interior", "edge", "ood", "interior_r0.002", "fig8", "fig8_r0.01"], groups));

%% Tables
fid = fopen(fullfile(sum_dir, 'paper_tables.tex'), 'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%% tab:validation-accuracy  median [min, max] across simulations\n');
fprintf(fid, '%% columns: R^2 all regions+times | NMSE all regions+times | median R^2(t) from day 120 | median b0\n');
for gi = 1:numel(groups)
    for mi = 1:numel(methods)
        r = C(C.group == groups(gi) & C.method == methods(mi), :);
        if isempty(r), continue, end
        lead = ""; if mi == 1, lead = glabels(gi); end
        fprintf(fid, '%-12s & %-7s & %s & %s & %s & %s \\\\  %% n=%d\n', lead, mlabels(mi), ...
            mrange(r.r2_pooled), mrange(r.nmse_pooled), mrange(r.r2_ss_median), mrange(r.slope0_median), height(r));
    end
    if gi < numel(groups), fprintf(fid, '\\midrule\n'); end
end

fprintf(fid, '\n%% tab:validation-coverage\n');
for gi = 1:numel(groups)
    r = V(V.group == groups(gi), :);
    w = r.n_edge_inputs / sum(r.n_edge_inputs); % pool edge inputs across the group's simulations
    fprintf(fid, '%-12s & %s & %s & %s & %s \\\\  %% n=%d\n', glabels(gi), num(max(r.N_max)), ...
        num(w' * r.frac_above_train_q99), num(w' * r.frac_above_train_q999), ...
        num(w' * r.frac_above_train_max), height(r));
end

fprintf(fid, '\n%% runtime, mean +- s.e.m. [s]\n');
num_rt = unique(C(:, {'setting_id', 'numeric_runtime_s'}), 'rows').numeric_runtime_s;
fprintf(fid, '%% numeric: %s +- %s (n=%d)\n', num(mean(num_rt)), num(sem(num_rt)), numel(num_rt));
for mi = 1:numel(methods)
    rt = C.runtime_s(C.method == methods(mi));
    if isempty(rt), continue, end
    fprintf(fid, '%% %s: %s +- %s (n=%d)\n', mlabels(mi), num(mean(rt)), num(sem(rt)), numel(rt));
end

fprintf(fid, '\n%% simulations with non-finite N: %d of %d\n', sum(C.nonfinite), height(C));
for k = find(C.nonfinite)'
    fprintf(fid, '%%   %s\n', C.case_id(k));
end
fprintf('wrote %s\n', fullfile(sum_dir, 'paper_tables.tex'));

%% Figure: per-simulation distributions
R2_FLOOR = -1; % pooled R^2 below this is drawn at the floor as an open marker
f = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 26 8.5]);
tl = tiledlayout(f, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
panels = {'r2_pooled', 'R^2_{1:1} (all regions and times)', 'linear'; ...
          'nmse_pooled', 'NMSE (all regions and times)', 'log'; ...
          'slope0_median', 'median b_0', 'log'};
rng(0);
for p = 1:size(panels, 1)
    ax = nexttile(tl); hold(ax, 'on');
    for gi = 1:numel(groups)
        for mi = 1:numel(methods)
            y = C.(panels{p, 1})(C.group == groups(gi) & C.method == methods(mi));
            if isempty(y), continue, end
            x = gi + (mi - 3) * 0.15 + (rand(size(y)) - 0.5) * 0.08;
            clipped = strcmp(panels{p, 1}, 'r2_pooled') & y < R2_FLOOR;
            scatter(ax, x(~clipped), y(~clipped), 36, colors(mi, :), 'filled', markers{mi}, ...
                'MarkerEdgeColor', 'w', 'LineWidth', 0.5, 'HandleVisibility', 'off');
            scatter(ax, x(clipped), R2_FLOOR + 0 * y(clipped), 36, colors(mi, :), markers{mi}, ...
                'LineWidth', 1, 'HandleVisibility', 'off');
            plot(ax, gi + (mi - 3) * 0.15 + [-0.06 0.06], median(y) * [1 1], '-', ...
                'Color', colors(mi, :), 'LineWidth', 2, 'HandleVisibility', 'off');
        end
    end
    set(ax, 'YScale', panels{p, 3}, 'XTick', 1:numel(groups), 'XTickLabel', glabels, ...
        'XLim', [0.4 numel(groups) + 0.6], 'TickDir', 'out', 'Box', 'off', 'FontSize', 9);
    ax.YGrid = 'on'; ax.GridAlpha = 0.12;
    ylabel(ax, panels{p, 2});
    switch panels{p, 1}
        case 'r2_pooled', yline(ax, 1, ':', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
            ylim(ax, [R2_FLOOR - 0.05, 1.02]);
        case 'slope0_median', yline(ax, 1, ':', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
    end
end
for mi = 1:numel(methods) % legend entries only
    scatter(ax, nan, nan, 36, colors(mi, :), 'filled', markers{mi}, 'DisplayName', mlabels(mi));
end
lg = legend(ax, 'Orientation', 'horizontal', 'Box', 'off'); lg.Layout.Tile = 'north';
save_fig(f, fullfile(sum_dir, 'fig_validation_dist'));

%% Figure: accuracy over simulated time
f = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 26 5.5 * numel(groups)]);
tl = tiledlayout(f, numel(groups), 3, 'TileSpacing', 'compact', 'Padding', 'compact');
for gi = 1:numel(groups)
    for p = 1:3
        ax = nexttile(tl); hold(ax, 'on');
        if p < 3
            metric = ["nmse", "slope0"]; metric = metric(p);
            for mi = 1:numel(methods)
                r = T(T.group == groups(gi) & T.method == methods(mi), :);
                if isempty(r), continue, end
                [days, q] = time_quantiles(r.days, r.(metric));
                fill(ax, [days; flipud(days)], [q(:, 1); flipud(q(:, 3))], colors(mi, :), ...
                    'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                plot(ax, days, q(:, 2), '-', 'Color', colors(mi, :), 'LineWidth', 1.5, 'DisplayName', mlabels(mi));
            end
        else % spatial variation of the numeric solution, one value per setting and time point
            r = T(T.group == groups(gi), :);
            [~, first] = unique(r.setting_id + "|" + string(r.days));
            r = r(first, :);
            [days, q] = time_quantiles(r.days, r.num_cv);
            fill(ax, [days; flipud(days)], [q(:, 1); flipud(q(:, 3))], [0.4 0.4 0.4], ...
                'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
            plot(ax, days, q(:, 2), '-', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5, 'HandleVisibility', 'off');
        end
        switch p
            case 1, ylabel(ax, 'NMSE');
                title(ax, glabels(gi), 'FontWeight', 'normal', 'HorizontalAlignment', 'left', 'Units', 'normalized', 'Position', [0 1.02]);
            case 2, yline(ax, 1, ':', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off'); ylabel(ax, 'b_0');
            case 3, ylabel(ax, 'CV of numeric N'); title(ax, 'regional spread', 'FontWeight', 'normal', 'FontSize', 9);
        end
        set(ax, 'YScale', 'log', 'TickDir', 'out', 'Box', 'off', 'FontSize', 9, 'XLim', [0 180]);
        ax.YGrid = 'on'; ax.GridAlpha = 0.12;
        if gi == numel(groups), xlabel(ax, 'time (days)'); end
    end
end
lg = legend(nexttile(tl, 1), 'Orientation', 'horizontal', 'Box', 'off'); lg.Layout.Tile = 'north';
save_fig(f, fullfile(sum_dir, 'fig_validation_time'));
end

function [days, q] = time_quantiles(t, v)
% 25th, 50th and 75th percentiles of v across simulations at each time point
[days, ~, k] = unique(t);
q = [accumarray(k, v, [], @(x) prctile_(x, 25)), accumarray(k, v, [], @(x) prctile_(x, 50)), ...
     accumarray(k, v, [], @(x) prctile_(x, 75))];
end

function s = mrange(v)
% Over finite values; simulations with non-finite N are listed separately in the file
v = v(isfinite(v));
if isempty(v), s = 'n/a'; return, end
s = sprintf('%s [%s, %s]', num(median(v)), num(min(v)), num(max(v)));
end

function s = num(x)
if ~isfinite(x), s = 'n/a'; elseif x == 0, s = '0';
elseif abs(x) >= 1e4 || abs(x) < 1e-3, s = sprintf('$%.2g$', x); s = regexprep(s, 'e([+-])0*(\d+)', '\\cdot10^{$1$2}'); s = strrep(s, '{+', '{');
else, s = sprintf('%.3g', x);
end
end

function s = sem(v)
s = std(v) / sqrt(numel(v));
end

function q = prctile_(v, p) % sorted-index percentile; no Statistics Toolbox needed
v = sort(v(isfinite(v)));
if isempty(v), q = nan; else, q = v(max(1, ceil(p / 100 * numel(v)))); end
end

function rgb = hex2rgb(hex)
rgb = zeros(numel(hex), 3);
for k = 1:numel(hex)
    h = char(hex(k)); rgb(k, :) = hex2dec({h(2:3), h(4:5), h(6:7)})' / 255;
end
end

function save_fig(f, base)
exportgraphics(f, base + ".pdf", 'ContentType', 'vector');
exportgraphics(f, base + ".png", 'Resolution', 200);
fprintf('wrote %s.pdf/.png\n', base);
end

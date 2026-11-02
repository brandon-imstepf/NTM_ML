function plot_t60_fits(grid, search_dir)
% Tau concentration over time for the 60-day recoveries (run_nm_t60.slurm):
%   truths_t60.png         every region of each numeric truth, default grid vs
%                          0.225 days (is the truth smooth?)
%   truth_wiggle.png       the region with the most changes of direction in each
%                          truth, every grid, numeric and MLP
%   fit_<param>.png        per truth: the numeric truth, the MLP at the true
%                          values and the MLP at the recovered value, summed
%                          over regions
%   fit_nodes_<param>.png  the same at day 59.4, region by region
% grid        truth grid of the fits: base, dt0.45 or dt0.225 (default dt0.225)
% search_dir  folder of the numeric-truth searches (default
%             outputs_nm/t60_<grid>/numeric_truth); figures go to <search_dir>/../figures
%
%   plot_t60_fits()
%   plot_t60_fits('dt0.225', 'outputs_nm/t60_dt0.225_near/numeric_truth')

here = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(grid), grid = "dt0.225"; end
grid = string(grid);
if nargin < 2 || isempty(search_dir)
    search_dir = fullfile(here, 'outputs_nm', "t60_" + grid, 'numeric_truth');
end
fig_dir = fullfile(fileparts(search_dir), 'figures');
if ~isfolder(fig_dir), mkdir(fig_dir); end
T = readtable(fullfile(here, 'manifests', 'recovery_truths.csv'), 'TextType', 'string');
T = T(T.convention == "l2f", :);
params = ["gamma1", "lambda1", "delta", "epsilon"];
param_tex = ["\gamma_1", "\lambda_1", "\delta", "\epsilon"];
ink = [0.23 0.23 0.22];  gray = [0.55 0.55 0.53];
blue = double([0x2a 0x78 0xd6]) / 255;  orange = double([0xeb 0x68 0x34]) / 255;
truth_file = @(g, m, s) fullfile(here, "outputs_t60_" + g, "rec_l2f__" + m + "__ECl__" + s + ".mat");

%% Truths: every region, default grid vs 0.225 days
f = figure('Visible', 'off', 'Color', 'w', 'Position', [0 0 1400 640]);
tl = tiledlayout(f, 2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
for s = 1:height(T)
    ax = clean_axes(nexttile(tl));  hold(ax, 'on');
    b = load(truth_file("base", "numeric", T.truth_id(s)), 'N', 'trange');
    r = load(truth_file("dt0.225", "numeric", T.truth_id(s)), 'N', 'trange');
    hb = plot(ax, b.trange * 180, 1e3 * b.N', '-', 'Color', orange, 'LineWidth', 1);
    hr = plot(ax, r.trange * 180, 1e3 * r.N', '-', 'Color', ink, 'LineWidth', 0.75);
    title(ax, short_label(T.truth_id(s)), 'FontWeight', 'normal', 'FontSize', 9);
    xlim(ax, [0 60]);
    if mod(s, 4) == 1, ylabel(ax, 'tau concentration N (	imes10^{-3})'); end
    if s > 4, xlabel(ax, 'day'); end
end
lg = legend([hb(1), hr(1)], ["numeric truth, default grid", "numeric truth, 0.225-day steps"], ...
    'Orientation', 'horizontal', 'Box', 'off');
lg.Layout.Tile = 'south';
title(tl, 'Numeric truths, all 30 regions, 60 days', 'FontSize', 13);
exportgraphics(f, fullfile(fig_dir, 'truths_t60.png'), 'Resolution', 150);
close(f);

%% The most wiggly region of each truth: every grid, numeric and MLP
grids = ["base", "dt0.45", "dt0.225"];
gcol = [blue; orange; double([0x1b 0xaf 0x7a]) / 255];
f = figure('Visible', 'off', 'Color', 'w', 'Position', [0 0 1400 640]);
tl = tiledlayout(f, 2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
for s = 1:height(T)
    ax = clean_axes(nexttile(tl));  hold(ax, 'on');
    r = load(truth_file("base", "numeric", T.truth_id(s)), 'N', 'region_names');
    d = diff(r.N, 1, 2);
    [n_changes, i] = max(sum(d(:, 1:end-1) .* d(:, 2:end) < 0, 2));
    h = gobjects(0);
    for g = 1:numel(grids)
        x = load(truth_file(grids(g), "numeric", T.truth_id(s)), 'N', 'trange');
        h(g) = plot(ax, x.trange * 180, 1e4 * x.N(i, :), '-', 'Color', gcol(g, :), 'LineWidth', 1.25); %#ok<AGROW>
    end
    m = load(truth_file("dt0.225", "mlp", T.truth_id(s)), 'N', 'trange');
    hm = plot(ax, m.trange * 180, 1e4 * m.N(i, :), '--', 'Color', gray, 'LineWidth', 1.25);
    title(ax, {short_label(T.truth_id(s)), sprintf('%s (%d changes, default grid)', r.region_names{i}, n_changes)}, ...
        'FontWeight', 'normal', 'FontSize', 8);
    xlim(ax, [0 60]);
    if mod(s, 4) == 1, ylabel(ax, 'N (	imes10^{-4})'); end
    if s > 4, xlabel(ax, 'day'); end
end
lg = legend([h, hm], ["numeric, default grid", "numeric, 0.45 d", "numeric, 0.225 d", "MLP at the true values, 0.225 d"], ...
    'Orientation', 'horizontal', 'Box', 'off');
lg.Layout.Tile = 'south';
title(tl, 'Region with the most changes of direction in each numeric truth (default grid)', 'FontSize', 13);
exportgraphics(f, fullfile(fig_dir, 'truth_wiggle.png'), 'Resolution', 150);
close(f);

%% Fits: truth, MLP at the true values, MLP at the recovered value
for p = 1:numel(params)
    f1 = figure('Visible', 'off', 'Color', 'w', 'Position', [0 0 1400 640]);
    tl1 = tiledlayout(f1, 2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
    f2 = figure('Visible', 'off', 'Color', 'w', 'Position', [0 0 1400 700]);
    tl2 = tiledlayout(f2, 2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
    for s = 1:height(T)
        num = load(truth_file(grid, "numeric", T.truth_id(s)), 'N', 'trange');
        mlp = load(truth_file(grid, "mlp", T.truth_id(s)), 'N');
        S = load(fullfile(search_dir, "nm_l2f__mlp__" + T.truth_id(s) + "__" + params(p) + ".mat"), 'best', 'result');
        day = num.trange * 180;
        nmse = @(x) sum((x(:) - num.N(:)).^2) / sum(num.N(:).^2);
        lbl = [sprintf("numeric truth"), ...
            sprintf("MLP at true %s = %.3g (NMSE %.2g)", params(p), S.result.true_value, nmse(mlp.N)), ...
            sprintf("MLP at recovered %.3g (NMSE %.2g, exitflag %d)", S.result.recovered_value, S.best.nmse, S.result.exitflag)];

        ax = clean_axes(nexttile(tl1));  hold(ax, 'on');
        h = [plot(ax, day, sum(num.N, 1), '-', 'Color', ink, 'LineWidth', 2), ...
             plot(ax, day, sum(mlp.N, 1), '--', 'Color', gray, 'LineWidth', 1.5), ...
             plot(ax, day, sum(S.best.N, 1), '-', 'Color', blue, 'LineWidth', 1.5)];
        legend(ax, h, lbl, 'Location', 'southeast', 'Box', 'off', 'FontSize', 7);
        title(ax, short_label(T.truth_id(s)), 'FontWeight', 'normal', 'FontSize', 9);
        xlim(ax, [0 60]);
        if mod(s, 4) == 1, ylabel(ax, 'N summed over regions'); end
        if s > 4, xlabel(ax, 'day'); end

        ax = clean_axes(nexttile(tl2));  hold(ax, 'on');
        x = num.N(:, end);
        lim = [min(x(x > 0)) / 2, max(x) * 2];
        plot(ax, lim, lim, '-', 'Color', ink, 'HandleVisibility', 'off');
        h = [plot(ax, x, mlp.N(:, end), 'o', 'Color', gray, 'MarkerSize', 6, 'LineWidth', 1), ...
             plot(ax, x, S.best.N(:, end), 'o', 'MarkerFaceColor', blue, 'MarkerEdgeColor', 'w', 'MarkerSize', 7)];
        set(ax, 'XScale', 'log', 'YScale', 'log', 'XLim', lim, 'YLim', lim);
        legend(ax, h, lbl(2:3), 'Location', 'southeast', 'Box', 'off', 'FontSize', 7);
        title(ax, short_label(T.truth_id(s)), 'FontWeight', 'normal', 'FontSize', 9);
        if mod(s, 4) == 1, ylabel(ax, 'MLP N, day 59.4'); end
        if s > 4, xlabel(ax, 'numeric truth N, day 59.4'); end
    end
    title(tl1, sprintf('%s searched (others at their true values), %s truth grid: tau summed over the 30 regions', ...
        param_tex(p), grid), 'FontSize', 13);
    title(tl2, sprintf('%s searched, %s truth grid: each region at day 59.4 (diagonal = exact)', ...
        param_tex(p), grid), 'FontSize', 13);
    exportgraphics(f1, fullfile(fig_dir, "fit_" + params(p) + ".png"), 'Resolution', 150);
    exportgraphics(f2, fullfile(fig_dir, "fit_nodes_" + params(p) + ".png"), 'Resolution', 150);
    close(f1);  close(f2);
end
fprintf('figures in %s\n', fig_dir);
end

function ax = clean_axes(ax)
set(ax, 'Box', 'off', 'TickDir', 'out', 'FontSize', 9, 'XColor', [0.32 0.32 0.31], ...
    'YColor', [0.32 0.32 0.31], 'GridColor', [0.85 0.85 0.84], 'GridAlpha', 1);
grid(ax, 'on');
end

function label = short_label(id)
% g0.0005_l0.02_d80_e20 -> \gamma5e-4 \lambda0.02 \delta80 \epsilon20
t = split(id, "_");
label = sprintf('\\gamma%s \\lambda%s \\delta%s \\epsilon%s', ...
    regexprep(sprintf('%.0e', str2double(extractAfter(t(1), 1))), 'e-0', 'e-'), ...
    extractAfter(t(2), 1), extractAfter(t(3), 1), extractAfter(t(4), 1));
end

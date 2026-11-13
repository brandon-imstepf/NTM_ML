function T = plot_fig8_gp(out_file, src_fmt, fig_title, num_file)
% Fig. 8 of the manuscript (with the GP, 11/10): regional tau of each sNTM
% against the numeric NTM at months 1-5, at the Fig. 8 parameters (gamma1 5e-5,
% lambda 0.01, delta 10, epsilon 50, seed ECl, init_rescale 0.01; the runs of
% manifests/fig8_gp.csv in outputs_fig8/, 180 days).
% One point per region (30). Red line: least-squares fit with intercept; dashed
% line: y = x. Both axes of a panel have the same range (square panels).
% R^2 of that fit (squared correlation) and the NMSE of Eq. 7 over the regions
% at that month (nmse_eq7.m).
% Concentrations in nM (N is stored in uM).
%
%   plot_fig8_gp();                 % figures/fig8_gp.png and .pdf
%   plot_fig8_gp('somewhere/fig8');
%   plot_fig8_gp(out, 'outputs_consistent/rec_l2f__%s__ECl__g0.0005_l0.02_d20_e80.mat', title)
%       same figure for another run set: src_fmt is relative to network_validation/
%       with %s for the method (numeric, mlp, gp, sr, lr); see plot_fig8_interior.m
%   plot_fig8_gp(out, src_fmt, title, num_file)
%       numeric run from num_file (relative to network_validation/) instead of
%       src_fmt; see plot_fig8_sets.m
%   T: table of method, month, R2, NMSE, slope

here = fileparts(mfilename('fullpath'));
nv = fileparts(here);
if nargin < 1, out_file = fullfile(here, 'figures', 'fig8_gp'); end
if ~isfolder(fileparts(out_file)), mkdir(fileparts(out_file)); end
if nargin < 2 || isempty(src_fmt), src_fmt = 'outputs_fig8/fig8__%s__ECl__g5e-05_l0.01_d10_e50.mat'; end
if nargin < 3, fig_title = 'Network wide NTM simulation regional tau accuracy'; end
file = @(m) fullfile(nv, sprintf(src_fmt, m));
if nargin >= 4 && ~isempty(num_file), num_path = fullfile(nv, num_file); else, num_path = file("numeric"); end
methods = ["mlp", "gp", "sr", "lr", "interp"];
labels = ["MLP", "GP", "Symbolic Regression", "Linear Regression", "Nearest-Neighbor"];
% the nearest-neighbor (interp) row only where that run exists (so far the Primary
% Parameter Set, outputs_consistent/rec_l2f__interp__CA1__g0.005_l0.02_d20_e80.mat)
has = arrayfun(@(m) isfile(file(m)), methods);  has(1:4) = true;
methods = methods(has);  labels = labels(has);
months = 1:5;

num = load(num_path, 'N', 'trange');
T = table('Size', [0 5], 'VariableTypes', ["string", "double", "double", "double", "double"], ...
    'VariableNames', ["method", "month", "R2", "NMSE", "slope"]);
day = num.trange * 180;
fig = figure('Units', 'centimeters', 'Position', [1 1 26 5.5 * numel(methods)], 'Color', 'w');
tl = tiledlayout(fig, numel(methods), numel(months), 'TileSpacing', 'compact', 'Padding', 'compact');
for r = 1:numel(methods)
    s = load(file(methods(r)), 'N');
    for c = 1:numel(months)
        [~, h] = min(abs(day - 30 * months(c)));
        x = num.N(:, h) * 1e3;  y = s.N(:, h) * 1e3;
        ax = nexttile(tl);  hold(ax, 'on');
        if ~all(isfinite(y))
            % the sNTM simulation diverged (NaN/Inf) before this month
            d_fail = day(find(any(~isfinite(s.N), 1), 1));
            text(ax, 0.5, 0.5, sprintf('simulation diverged\n(day %.1f)', d_fail), 'Units', 'normalized', ...
                'HorizontalAlignment', 'center', 'FontSize', 9, 'FontName', 'Times', 'Color', [0.6 0 0]);
            set(ax, 'FontSize', 8, 'FontName', 'Times', 'Box', 'on', 'XTick', [], 'YTick', []);
            if r == 1, title(ax, sprintf('Month %d', months(c)), 'FontSize', 12, 'FontName', 'Times'); end
            if c == 1, ylabel(ax, labels(r), 'FontSize', 11, 'FontWeight', 'bold', 'FontName', 'Times'); end
            fprintf('%-4s month %d: diverged at day %.1f\n', methods(r), months(c), d_fail);
            T(end+1, :) = {methods(r), months(c), NaN, NaN, NaN}; %#ok<AGROW>
            continue
        end
        % same range on both axes (square panel), so a slope of 1 looks like one
        lim = [min([x; y]), max([x; y])];  lim = lim + [-1 1] * 0.06 * diff(lim);
        h_id = plot(ax, lim, lim, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.8);    % y = x
        scatter(ax, x, y, 12, [0.12 0.47 0.71], 'filled');
        p = polyfit(x, y, 1);
        xl = [min(x), max(x)];
        h_fit = plot(ax, xl, polyval(p, xl), 'r-', 'LineWidth', 1);
        xlim(ax, lim);  ylim(ax, lim);  axis(ax, 'square');
        if r == 1 && c == 1
            lg = legend(ax, [h_id, h_fit], ["y = x", "Least-squares fit"], 'Orientation', 'horizontal', ...
                'FontName', 'Times', 'FontSize', 11, 'Box', 'off');
            lg.Layout.Tile = 'south';
        end
        R = corrcoef(x, y);
        nmse = nmse_eq7(y, x);
        set(ax, 'FontSize', 8, 'FontName', 'Times', 'Box', 'on', 'TickDir', 'out');
        % R^2 and NMSE above the panel (inside it they would cover points far from y = x)
        stats = sprintf('R^2 = %.3f,  NMSE = %.3g', R(1, 2)^2, nmse);
        if r == 1
            title(ax, sprintf('Month %d', months(c)), 'FontSize', 12, 'FontName', 'Times');
            subtitle(ax, stats, 'FontSize', 9, 'FontName', 'Times');
        else
            title(ax, stats, 'FontSize', 9, 'FontName', 'Times', 'FontWeight', 'normal');
        end
        if c == 1, ylabel(ax, labels(r), 'FontSize', 11, 'FontWeight', 'bold', 'FontName', 'Times'); end
        fprintf('%-4s month %d: R^2 %.4f, NMSE %.3g, slope %.3f\n', methods(r), months(c), R(1, 2)^2, nmse, p(1));
        T(end+1, :) = {methods(r), months(c), R(1, 2)^2, nmse, p(1)}; %#ok<AGROW>
    end
end
xlabel(tl, '[Tau] (nM) - Numeric NTM', 'FontName', 'Times', 'FontSize', 12);
ylabel(tl, '[Tau] (nM) - sNTM', 'FontName', 'Times', 'FontSize', 12);
title(tl, fig_title, 'FontName', 'Times', 'FontSize', 15, 'FontWeight', 'bold');
exportgraphics(fig, out_file + ".png", 'Resolution', 300);
exportgraphics(fig, out_file + ".pdf", 'ContentType', 'vector');
fprintf('wrote %s.png/.pdf\n', out_file);
end

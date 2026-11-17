function T = plot_fig8_gp(out_file, src_fmt, fig_title, num_file)
% Fig. 8 of the manuscript (with the GP, 11/10): regional tau of each sNTM
% against the numeric NTM at months 1-5, at the Fig. 8 parameters (gamma1 5e-5,
% lambda 0.01, delta 10, epsilon 50, seed ECl, init_rescale 0.01; the runs of
% manifests/fig8_gp.csv in outputs_fig8/, 180 days).
% One point per region (30). Red line (11/17, advisor): least-squares fit through
% the origin, y = b x with slope b = sum(x y) / sum(x^2), so that b - 1 is the
% proportional bias of the sNTM; dashed line: y = x. Both axes of a panel have the
% same range (square panels). Each panel gives b, the R^2 of that fit,
% 1 - sum (y - b x)^2 / sum (y - mean y)^2, and the NMSE of Eq. 7 over the regions
% at that month (nmse_eq7.m). (Before 11/17: fit with intercept, R^2 = squared
% correlation, intercept marked by a diamond/triangle.)
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
%   T: table of method, month, R2 (of the fit through the origin), NMSE, slope b

here = fileparts(mfilename('fullpath'));
nv = fileparts(here);
if nargin < 1, out_file = fullfile(here, 'figures', 'fig8_gp'); end
if ~isfolder(fileparts(out_file)), mkdir(fileparts(out_file)); end
if nargin < 2 || isempty(src_fmt), src_fmt = 'outputs_fig8/fig8__%s__ECl__g5e-05_l0.01_d10_e50.mat'; end
if nargin < 3, fig_title = 'Network wide NTM simulation regional tau accuracy'; end
file = @(m) fullfile(nv, sprintf(src_fmt, m));
if nargin >= 4 && ~isempty(num_file), num_path = fullfile(nv, num_file); else, num_path = file("numeric"); end
% rows in the paper's model order (11/17): LR, SR, 1-NN, GP, MLP
methods = ["lr", "sr", "interp", "gp", "mlp"];
labels = ["Linear Regression", "Symbolic Regression", "1-NN Interpolation", "Gaussian Process", "MLP"];
% the nearest-neighbor (interp) row only where that run exists (so far the Primary
% Parameter Set, outputs_consistent/rec_l2f__interp__CA1__g0.005_l0.02_d20_e80.mat)
has = arrayfun(@(m) isfile(file(m)), methods);  has(methods ~= "interp") = true;
methods = methods(has);  labels = labels(has);
months = 1:5;
fs = 1.25;          % all font sizes x1.25 (11/16)

num = load(num_path, 'N', 'trange');
T = table('Size', [0 5], 'VariableTypes', ["string", "double", "double", "double", "double"], ...
    'VariableNames', ["method", "month", "R2", "NMSE", "slope"]);
day = num.trange * 180;
% one axis range for every panel (11/13): from 0 to the largest finite
% concentration of any model at any of the months shown, so that panels can be
% compared
hm = arrayfun(@(m) find(abs(day - 30 * m) == min(abs(day - 30 * m)), 1), months);
S = cell(size(methods));  vmax = max(num.N(:, hm), [], 'all');
for r = 1:numel(methods)
    S{r} = load(file(methods(r)), 'N');
    v = S{r}.N(:, hm);  v = v(isfinite(v));
    if ~isempty(v), vmax = max(vmax, max(v)); end
end
lim = [0, 1.05 * vmax * 1e3];
fig = figure('Units', 'centimeters', 'Position', [1 1 29 5.6 * numel(methods)], 'Color', 'w');
tl = tiledlayout(fig, numel(methods), numel(months), 'TileSpacing', 'compact', 'Padding', 'compact');
for r = 1:numel(methods)
    s = S{r};
    for c = 1:numel(months)
        h = hm(c);
        x = num.N(:, h) * 1e3;  y = s.N(:, h) * 1e3;
        ax = nexttile(tl);  hold(ax, 'on');
        if ~all(isfinite(y))
            % the sNTM simulation diverged (NaN/Inf) before this month
            d_fail = day(find(any(~isfinite(s.N), 1), 1));
            text(ax, 0.5, 0.5, sprintf('simulation diverged\n(day %.1f)', d_fail), 'Units', 'normalized', ...
                'HorizontalAlignment', 'center', 'FontSize', 9 * fs, 'FontName', 'Times', 'Color', [0.6 0 0]);
            set(ax, 'FontSize', 8 * fs, 'FontName', 'Times', 'Box', 'on', 'XTick', [], 'YTick', []);
            if r == 1, title(ax, sprintf('Month %d', months(c)), 'FontSize', 12 * fs, 'FontName', 'Times'); end
            if c == 1, ylabel(ax, labels(r), 'FontSize', 11 * fs, 'FontWeight', 'bold', 'FontName', 'Times'); end
            fprintf('%-4s month %d: diverged at day %.1f\n', methods(r), months(c), d_fail);
            T(end+1, :) = {methods(r), months(c), NaN, NaN, NaN}; %#ok<AGROW>
            continue
        end
        h_id = plot(ax, lim, lim, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.8);    % y = x
        scatter(ax, x, y, 12, [0.12 0.47 0.71], 'filled');
        b = sum(x .* y) / sum(x .^ 2);                      % least squares through the origin
        xl = [0, max(x)];
        h_fit = plot(ax, xl, b * xl, 'r-', 'LineWidth', 1);
        r2 = 1 - sum((y - b * x) .^ 2) / sum((y - mean(y)) .^ 2);
        xlim(ax, lim);  ylim(ax, lim);  axis(ax, 'square');
        % x ticks every 0.4 (every 0.2 do not fit side by side at this font size)
        set(ax, 'XTick', 0:0.4:lim(2), 'YTick', 0:0.2:lim(2), 'XTickLabelRotation', 0);
        if r == 1 && c == 1
            lg = legend(ax, [h_id, h_fit], ["y = x", "Least-squares fit through the origin, y = bx"], ...
                'Orientation', 'horizontal', 'FontName', 'Times', 'FontSize', 11 * fs, 'Box', 'off');
        end
        nmse = nmse_eq7(y, x);
        set(ax, 'FontSize', 8 * fs, 'FontName', 'Times', 'Box', 'on', 'TickDir', 'out');
        % slope, R^2 and NMSE in a small box at the bottom right
        text(ax, 0.97, 0.04, {sprintf('b = %.3f', b), sprintf('R^2 = %.3f', r2), sprintf('NMSE = %.3g', nmse)}, ...
            'Units', 'normalized', 'HorizontalAlignment', 'right', ...
            'VerticalAlignment', 'bottom', 'FontName', 'Times', 'FontSize', 7.5 * fs, 'BackgroundColor', 'w', ...
            'EdgeColor', [0.6 0.6 0.6], 'Margin', 1.5);
        if r == 1, title(ax, sprintf('Month %d', months(c)), 'FontSize', 12 * fs, 'FontName', 'Times'); end
        if c == 1, ylabel(ax, labels(r), 'FontSize', 11 * fs, 'FontWeight', 'bold', 'FontName', 'Times'); end
        fprintf('%-6s month %d: b %.3f, R^2 %.4f, NMSE %.3g\n', methods(r), months(c), b, r2, nmse);
        T(end+1, :) = {methods(r), months(c), r2, nmse, b}; %#ok<AGROW>
    end
end
% x-label and legend (11/16): the legend is one row in the layout's bottom slot,
% its items spread across the width, and the x-label is the legend's title, so
% that it sits directly above the legend (as xlabel(tl) it came below it)
lg.String = lg.String + "        ";
lg.ItemTokenSize = [30 18];
lg.Layout.Tile = 'south';
title(lg, '[Tau] (nM) - Numeric NTM', 'FontName', 'Times', 'FontSize', 12 * fs, 'FontWeight', 'normal');
ylabel(tl, '[Tau] (nM) - sNTM', 'FontName', 'Times', 'FontSize', 12 * fs);
title(tl, fig_title, 'FontName', 'Times', 'FontSize', 15 * fs, 'FontWeight', 'bold');
exportgraphics(fig, out_file + ".png", 'Resolution', 300);
exportgraphics(fig, out_file + ".pdf", 'ContentType', 'vector');
fprintf('wrote %s.png/.pdf\n', out_file);
end

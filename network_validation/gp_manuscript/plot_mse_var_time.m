function T = plot_mse_var_time(out_file, src_fmt, fig_title)
% Why the NMSE of the network simulations grows over time (11/18, Suzanne):
% the NMSE over the regions at time t (Eq. 7, nmse_eq7.m) is MSE(t) / Var(t), the
% mean squared error of the sNTM regional concentrations divided by the variance
% of the numeric regional concentrations. Per model (panels, paper order LR, SR,
% 1-NN, GP, MLP from left to right): MSE(t) of the sNTM, Var(t) of the numeric NTM
% (the same in every panel) and their ratio, the NMSE, against time, on a
% logarithmic axis. Concentrations in nM.
%
%   plot_mse_var_time();   % Primary Parameter Set -> figures/mse_var_primary.png/.pdf
%   T: method, day, mse, var_num, nmse

here = fileparts(mfilename('fullpath'));
nv = fileparts(here);
addpath(nv);
if nargin < 1, out_file = fullfile(here, 'figures', 'mse_var_primary'); end
if nargin < 2, src_fmt = 'outputs_consistent/rec_l2f__%s__CA1__g0.005_l0.02_d20_e80.mat'; end
if nargin < 3, fig_title = 'Error and variance of the regional concentrations over time'; end
methods = ["lr", "sr", "interp", "gp", "mlp"];
labels = ["LR", "SR", "1-NN", "GP", "MLP"];
num = load(fullfile(nv, sprintf(src_fmt, 'numeric')), 'N', 'trange');
day = num.trange * 180;
X = num.N * 1e3;
v = mean((X - mean(X, 1)) .^ 2, 1);              % variance over the regions at each time
k = day > 0;                                     % t = 0: all tau in the seed, both runs equal
T = table();
W = 13.1;  H = 9.5;
fig = figure('Units', 'centimeters', 'Position', [1 1 W H], 'Color', 'w', ...
    'PaperUnits', 'centimeters', 'PaperSize', [W H], 'PaperPosition', [0 0 W H]);
tl = tiledlayout(fig, 2, numel(methods), 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:numel(methods)
    s = load(fullfile(nv, sprintf(src_fmt, methods(i))), 'N');
    Y = s.N * 1e3;
    mse = mean((Y - X) .^ 2, 1);
    % top row: MSE of the sNTM and variance of the numeric NTM (nM^2)
    ax = nexttile(tl, i);  hold(ax, 'on');
    h1 = plot(ax, day(k), mse(k), '-', 'Color', [0.85 0.1 0.1], 'LineWidth', 1.2);
    h2 = plot(ax, day(k), v(k), '-', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.2);
    set(ax, 'YScale', 'log', 'FontName', 'Times', 'FontSize', 9.5, 'Box', 'on', 'TickDir', 'out', ...
        'XLim', [0 180], 'XTick', 0:90:180, 'XTickLabel', [], 'YLim', [1e-6 1], 'YTick', 10.^(-6:2:0));
    if i > 1, ax.YTickLabel = []; else, ylabel(ax, 'nM^2', 'FontName', 'Times', 'FontSize', 11); end
    title(ax, labels(i), 'FontName', 'Times', 'FontSize', 12);
    % bottom row: their ratio, the NMSE (Eq. 7)
    ax2 = nexttile(tl, numel(methods) + i);
    h3 = plot(ax2, day(k), mse(k) ./ v(k), '-', 'Color', [0.12 0.47 0.71], 'LineWidth', 1.2);
    set(ax2, 'YScale', 'log', 'FontName', 'Times', 'FontSize', 9.5, 'Box', 'on', 'TickDir', 'out', ...
        'XLim', [0 180], 'XTick', 0:90:180, 'XTickLabelRotation', 0, 'YLim', [1e-3 1e3], 'YTick', 10.^(-3:2:3));
    if i > 1, ax2.YTickLabel = []; else, ylabel(ax2, 'NMSE', 'FontName', 'Times', 'FontSize', 11); end
    T = [T; table(repmat(methods(i), nnz(k), 1), day(k)', mse(k)', v(k)', (mse(k) ./ v(k))', ...
        'VariableNames', ["method", "day", "mse", "var_num", "nmse"])]; %#ok<AGROW>
    fprintf('%-6s day 30/90/150: MSE %s, Var %s, NMSE %s\n', methods(i), ...
        mat2str(interp1(day, mse, [30 90 150]), 2), mat2str(interp1(day, v, [30 90 150]), 2), ...
        mat2str(interp1(day, mse ./ max(v, eps), [30 90 150]), 2));
end
if fig_title ~= "", title(tl, fig_title, 'FontName', 'Times', 'FontSize', 12, 'FontWeight', 'bold'); end
lg = legend([h1 h2 h3], ["MSE of the sNTM", "Variance of the numeric NTM", "NMSE = MSE / variance"], ...
    'NumColumns', 2, 'FontName', 'Times', 'FontSize', 10, 'Box', 'off');
lg.Layout.Tile = 'south';
% x-label above the legend, as its title (as in plot_fig8_gp)
title(lg, 'Time (days)', 'FontName', 'Times', 'FontSize', 12, 'FontWeight', 'normal');
writetable(T, fullfile(nv, 'results', 'mse_var_time_primary.csv'));
print(fig, out_file + ".png", '-dpng', '-r300');
print(fig, out_file + ".pdf", '-dpdf', '-vector');
fprintf('wrote %s.png/.pdf\n', out_file);
end

function plot_nn_k_curve(out_file)
% Candidate Results figure for nearest-neighbor interpolation with a complexity
% axis (as the number of features for linear regression, the complexity for
% symbolic regression and the epochs for the MLP): training and held-out NMSE
% (Eq. 7) of k-nearest-neighbor averaging against k, from results/nn_k_curve.csv
% (nn_k_curve.m). k decreases to the right (flexibility increases), so that the
% axis reads like the others; the network model (k = 1) is shown in red.
% Points: mean over 10^3 bootstrapped datasets of 10^3 samples; error bars: s.d.
%
%   plot_nn_k_curve();     % figures/results_variants/gd_results_k.png/.pdf

here = fileparts(mfilename('fullpath'));
nv = fileparts(here);
if nargin < 1, out_file = fullfile(here, 'figures', 'results_variants', 'gd_results_k'); end
T = readtable(fullfile(nv, 'results', 'nn_k_curve.csv'), 'TextType', 'string');
target = ["flux", "emc"];  tlab = ["\phi_\tau", "\Delta_\tau"];
c_train = [0 0.447 0.741];  c_test = [0.466 0.674 0.188];  c_sel = [1 0 0];
font = 'Times New Roman';
fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 20 9], 'Color', 'w');
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for t = 1:2
    a = sortrows(T(T.target == target(t), :), 'k');
    ax = nexttile(tl);  hold(ax, 'on');  box(ax, 'on');
    h1 = errorbar(ax, a.k, a.train_nmse_mean, a.train_nmse_sd, '-o', 'Color', c_train, 'MarkerFaceColor', c_train, ...
        'MarkerSize', 4, 'LineWidth', 1.4, 'CapSize', 3);
    h2 = errorbar(ax, a.k, a.heldout_nmse_mean, a.heldout_nmse_sd, '-o', 'Color', c_test, 'MarkerFaceColor', c_test, ...
        'MarkerSize', 4, 'LineWidth', 1.4, 'CapSize', 3);
    k1 = a.k == 1;
    h3 = plot(ax, [1 1], [a.train_nmse_mean(k1), a.heldout_nmse_mean(k1)], 'o', 'MarkerEdgeColor', c_sel, ...
        'MarkerFaceColor', c_sel, 'MarkerSize', 6);
    [~, kb] = min(a.heldout_nmse_mean);
    xline(ax, a.k(kb), ':', 'Color', [0.4 0.4 0.4], 'LineWidth', 1);
    set(ax, 'XScale', 'log', 'XDir', 'reverse', 'XLim', [0.8 160], 'YLim', [0 0.05], 'XTick', [1 2 5 10 20 50 100], 'XTickLabel', ["1", "2", "5", "10", "20", "50", "100"], 'XMinorTick', 'off', 'FontName', font, ...
        'FontSize', 12, 'TickDir', 'out');
    title(ax, tlab(t), 'FontName', font, 'FontSize', 16, 'FontWeight', 'normal');
    xlabel(ax, 'Number of neighbors averaged, k', 'FontName', font, 'FontSize', 12);
    if t == 1, ylabel(ax, 'NMSE', 'FontName', font, 'FontSize', 13); else, ax.YTickLabel = []; end
    if t == 2
        legend(ax, [h1 h2 h3(1)], {'Training', 'Held-out test', 'Network model (k = 1)'}, 'Location', 'northwest', ...
            'FontName', font, 'FontSize', 10);
    end
end
title(tl, 'Nearest-neighbor interpolation: number of neighbors averaged', 'FontName', font, 'FontSize', 16, 'FontWeight', 'bold');
exportgraphics(fig, out_file + ".png", 'Resolution', 300, 'BackgroundColor', 'w');
exportgraphics(fig, out_file + ".pdf", 'ContentType', 'vector', 'BackgroundColor', 'w');
close(fig);
fprintf('wrote %s\n', out_file);
end

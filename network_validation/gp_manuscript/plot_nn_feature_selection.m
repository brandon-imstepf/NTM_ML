function plot_nn_feature_selection(out_file)
% Results figure for the nearest-neighbor surrogate with the number of inputs as
% the complexity axis, as the forward feature selection of linear regression
% (Fig. 3): leave-one-out (training) and held-out NMSE (Eq. 7) of 1-nearest-
% neighbor interpolation against the number of inputs, the inputs added in the
% order of forward selection (results/nn_feature_selection.csv,
% nn_feature_selection.m). The network model (all six inputs) is shown in red.
% Points: mean over 10^3 bootstrapped datasets of 10^3 samples; error bars: s.d.
%
%   plot_nn_feature_selection();   % figures/results_variants/gd_results_inputs.png/.pdf

here = fileparts(mfilename('fullpath'));
nv = fileparts(here);
if nargin < 1, out_file = fullfile(here, 'figures', 'results_variants', 'gd_results_inputs'); end
T = readtable(fullfile(nv, 'results', 'nn_feature_selection.csv'), 'TextType', 'string');
target = ["flux", "emc"];  tlab = ["\phi_\tau", "\Delta_\tau"];
sym = containers.Map(["gamma1", "lambda1", "delta", "epsilon", "N1", "N2"], ...
    ["\gamma", "\lambda", "\delta", "\epsilon", "N_1", "N_2"]);
c_loo = [0 0.447 0.741];  c_test = [0.466 0.674 0.188];  c_sel = [1 0 0];
font = 'Times New Roman';
fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 20 9.5], 'Color', 'w');
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for t = 1:2
    a = sortrows(T(T.target == target(t), :), 'n_inputs');
    ax = nexttile(tl);  hold(ax, 'on');  box(ax, 'on');
    h1 = errorbar(ax, a.n_inputs, a.loo_nmse_mean, a.loo_nmse_sd, '-o', 'Color', c_loo, 'MarkerFaceColor', c_loo, ...
        'MarkerSize', 4, 'LineWidth', 1.4, 'CapSize', 3);
    h2 = errorbar(ax, a.n_inputs, a.heldout_nmse_mean, a.heldout_nmse_sd, '-o', 'Color', c_test, ...
        'MarkerFaceColor', c_test, 'MarkerSize', 4, 'LineWidth', 1.4, 'CapSize', 3);
    h3 = plot(ax, [6 6], [a.loo_nmse_mean(end), a.heldout_nmse_mean(end)], 'o', 'MarkerEdgeColor', c_sel, ...
        'MarkerFaceColor', c_sel, 'MarkerSize', 6);
    lab = strings(1, 6);
    for k = 1:6, lab(k) = string(sym(a.added(k)));  if k > 1, lab(k) = "+" + lab(k); end, end
    set(ax, 'YScale', 'log', 'XLim', [0.6 6.4], 'YLim', [0.01 3], 'XTick', 1:6, 'XTickLabel', lab, ...
        'FontName', font, 'FontSize', 12, 'TickDir', 'out', 'TickLabelInterpreter', 'tex');
    title(ax, tlab(t), 'FontName', font, 'FontSize', 16, 'FontWeight', 'normal');
    xlabel(ax, 'Inputs, added in order of selection (1 to 6)', 'FontName', font, 'FontSize', 12);
    if t == 1, ylabel(ax, 'NMSE', 'FontName', font, 'FontSize', 13); else, ax.YTickLabel = []; end
    if t == 2
        legend(ax, [h1 h2 h3(1)], {'Leave-one-out (training)', 'Held-out test', 'Network model (six inputs)'}, ...
            'Location', 'southwest', 'FontName', font, 'FontSize', 10);
    end
end
title(tl, '1-Nearest-neighbor interpolation: input selection', 'FontName', font, 'FontSize', 16, 'FontWeight', 'bold');
exportgraphics(fig, out_file + ".png", 'Resolution', 300, 'BackgroundColor', 'w');
exportgraphics(fig, out_file + ".pdf", 'ContentType', 'vector', 'BackgroundColor', 'w');
close(fig);
fprintf('wrote %s\n', out_file);
end

function plot_gp_sigma_curve(out_file)
% Candidate Results figure for the GP with a complexity axis: training and
% held-out NMSE (Eq. 7) of the paper's GP configuration (5,000 samples) refitted
% with the noise standard deviation fixed at each value, from
% results/gp_sigma_curve.csv (gp_sigma_curve.m). The noise decreases to the right
% (flexibility increases), so that the axis reads like the number of features,
% complexity and epochs of the other surrogates; the paper's GP (noise 0.01, the
% lower bound at which its fitted noise ended) is shown in red.
% Points: mean over 10^3 bootstrapped datasets of 10^3 samples; error bars: s.d.
%
%   plot_gp_sigma_curve();     % figures/results_variants/gp_results_sigma.png/.pdf

here = fileparts(mfilename('fullpath'));
nv = fileparts(here);
if nargin < 1, out_file = fullfile(here, 'figures', 'results_variants', 'gp_results_sigma'); end
T = sortrows(readtable(fullfile(nv, 'results', 'gp_sigma_curve.csv')), 'sigma');
target = ["flux", "emc"];  tlab = ["\phi_\tau", "\Delta_\tau"];
c_train = [0 0.447 0.741];  c_test = [0.466 0.674 0.188];  c_sel = [1 0 0];
font = 'Times New Roman';
fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 20 9], 'Color', 'w');
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for t = 1:2
    ax = nexttile(tl);  hold(ax, 'on');  box(ax, 'on');
    h1 = errorbar(ax, T.sigma, T.("train_" + target(t) + "_mean"), T.("train_" + target(t) + "_sd"), '-o', ...
        'Color', c_train, 'MarkerFaceColor', c_train, 'MarkerSize', 4, 'LineWidth', 1.4, 'CapSize', 3);
    h2 = errorbar(ax, T.sigma, T.("heldout_" + target(t) + "_mean"), T.("heldout_" + target(t) + "_sd"), '-o', ...
        'Color', c_test, 'MarkerFaceColor', c_test, 'MarkerSize', 4, 'LineWidth', 1.4, 'CapSize', 3);
    k = abs(T.sigma - 0.01) < 1e-12;
    h3 = plot(ax, [0.01 0.01], [T.("train_" + target(t) + "_mean")(k), T.("heldout_" + target(t) + "_mean")(k)], 'o', ...
        'MarkerEdgeColor', c_sel, 'MarkerFaceColor', c_sel, 'MarkerSize', 6);
    set(ax, 'XScale', 'log', 'YScale', 'log', 'XDir', 'reverse', 'XLim', [5e-5 2], 'FontName', font, ...
        'FontSize', 12, 'TickDir', 'out');
    title(ax, tlab(t), 'FontName', font, 'FontSize', 16, 'FontWeight', 'normal');
    xlabel(ax, 'Noise standard deviation (fixed)', 'FontName', font, 'FontSize', 12);
    if t == 1, ylabel(ax, 'NMSE', 'FontName', font, 'FontSize', 13); end
    if t == 2
        legend(ax, [h1 h2 h3(1)], {'Training', 'Held-out test', 'Selected model'}, 'Location', 'best', ...
            'FontName', font, 'FontSize', 10);
    end
end
title(tl, 'Gaussian process: noise level', 'FontName', font, 'FontSize', 16, 'FontWeight', 'bold');
exportgraphics(fig, out_file + ".png", 'Resolution', 300, 'BackgroundColor', 'w');
exportgraphics(fig, out_file + ".pdf", 'ContentType', 'vector', 'BackgroundColor', 'w');
close(fig);
fprintf('wrote %s\n', out_file);
end

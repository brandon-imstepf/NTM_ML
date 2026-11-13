function plot_gp_results_fig(out_file)
% Results figure for the GP (main text), in the style of the linear regression
% feature-selection figure (images/Figure4.png): training and held-out NMSE
% (Eq. 7) of the GP against the number of training samples, for phi_tau and
% Delta_tau; the selected model (n = 5000) in red. GP only (the comparison with the other surrogates and griddatan is
% Supplementary Fig. S4, plot_supp_gp_size.m).
% Data: results/supp_gp_griddatan.csv (supp_gp_griddatan.m): the paper's GP
% configuration (exact, linear basis, ARD Matern 3/2, fminunc) fit on n rows of
% data_bias_e5_1.csv; n = 5000 is the selected model. Training NMSE on the n
% training rows; held-out NMSE on data_bias_e5_11.csv (10^5 samples never used
% for fitting or selection). Points: mean over 10^3 bootstrap resamples of 10^3
% samples; error bars: s.d.
%
%   plot_gp_results_fig();      % gp_manuscript/figures/gp_results.png/.pdf

here = fileparts(mfilename('fullpath'));
nv = fileparts(here);
if nargin < 1, out_file = fullfile(here, 'figures', 'gp_results'); end
T = readtable(fullfile(nv, 'results', 'supp_gp_griddatan.csv'), 'TextType', 'string');
gp = sortrows(T(startsWith(T.model, "gp_") & T.model ~= "gp_lrsr_1000", :), 'n_train');
n_sel = 5000;

target = ["flux", "emc"];  tlab = ["\phi_\tau", "\Delta_\tau"];
c_train = [0 0.447 0.741];  c_test = [0.466 0.674 0.188];  c_sel = [1 0 0];
font = 'Times New Roman';
fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 20 9], 'Color', 'w');
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for t = 1:2
    ax = nexttile(tl);  hold(ax, 'on');  box(ax, 'on');
    h1 = errorbar(ax, gp.n_train, gp.("train_" + target(t) + "_nmse_mean"), gp.("train_" + target(t) + "_nmse_sd"), ...
        '-o', 'Color', c_train, 'MarkerFaceColor', c_train, 'MarkerSize', 4, 'LineWidth', 1.4, 'CapSize', 3);
    h2 = errorbar(ax, gp.n_train, gp.("heldout_" + target(t) + "_nmse_mean"), gp.("heldout_" + target(t) + "_nmse_sd"), ...
        '-o', 'Color', c_test, 'MarkerFaceColor', c_test, 'MarkerSize', 4, 'LineWidth', 1.4, 'CapSize', 3);
    % selected model: its two points drawn in red
    k = gp.n_train == n_sel;
    h3 = plot(ax, gp.n_train(k), gp.("train_" + target(t) + "_nmse_mean")(k), 'o', ...
        'MarkerEdgeColor', c_sel, 'MarkerFaceColor', c_sel, 'MarkerSize', 6);
    plot(ax, gp.n_train(k), gp.("heldout_" + target(t) + "_nmse_mean")(k), 'o', ...
        'MarkerEdgeColor', c_sel, 'MarkerFaceColor', c_sel, 'MarkerSize', 6);
    set(ax, 'YScale', 'log', 'XLim', [0 5400], 'YLim', [1e-8 1], 'XTick', 0:1000:5000, ...
        'YTick', 10.^(-8:2:0), 'FontName', font, 'FontSize', 12, 'TickDir', 'out');
    title(ax, tlab(t), 'FontName', font, 'FontSize', 18, 'FontWeight', 'normal');
    if t == 1, ylabel(ax, 'NMSE', 'FontName', font, 'FontSize', 14); else, ax.YTickLabel = []; end
    if t == 2
        legend(ax, [h1 h2 h3], {'Training', 'Held-out test', 'Selected model'}, 'Location', 'east', ...
            'FontName', font, 'FontSize', 11);
    end
end
xlabel(tl, 'Number of GP Training Samples', 'FontName', font, 'FontSize', 14);
title(tl, 'Gaussian process training set size', 'FontName', font, 'FontSize', 18, 'FontWeight', 'bold');
exportgraphics(fig, out_file + ".png", 'Resolution', 300, 'BackgroundColor', 'w');
exportgraphics(fig, out_file + ".pdf", 'ContentType', 'vector', 'BackgroundColor', 'w');
close(fig);
fprintf('wrote %s.png/.pdf\n', out_file);
end

function fig8_swap_plot()
% Swap test figure (FIG8_OFFSET.md, step 3): total soluble tau sum_i N_i of each
% run divided by the numeric run, over the 180-day simulation (log time axis so
% days 0-10, where the offset forms, are visible).
%   full MLP sNTM, A (numeric phi_tau + MLP Delta_tau), B (MLP phi_tau + numeric
%   Delta_tau), and the sum of the deficits of A and B (dashed) for comparison
%   with the full MLP deficit.
% Writes outputs_fig8/figures/swap_test.png and .pdf.
%
%   fig8_swap_plot();

here = fileparts(mfilename('fullpath'));
f = @(m) fullfile(here, 'outputs_fig8', "fig8__" + m + "__ECl__g5e-05_l0.01_d10_e50.mat");
num = load(f("numeric"), 'N', 'trange');
day = num.trange(:)' * 180;
ratio = @(m) sum(getfield(load(f(m), 'N'), 'N'), 1) ./ sum(num.N, 1);
r_mlp = ratio("mlp");
r_A = ratio("swapA_numflux_mlpw1");
r_B = ratio("swapB_mlpflux_numw1");
r_sum = 1 - ((1 - r_A) + (1 - r_B));
keep = day > 0;

col = {[42 120 214]/255, [235 104 52]/255, [27 175 122]/255};   % blue, orange, aqua
ink = [0.25 0.25 0.25];  grid_c = [0.85 0.85 0.85];
fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'centimeters', 'Position', [0 0 17 12]);
ax = axes(fig);  hold(ax, 'on');
yline(ax, 1, '-', 'Color', [0.55 0.55 0.55], 'LineWidth', 1, 'HandleVisibility', 'off');
h(1) = plot(ax, day(keep), r_mlp(keep), '-', 'Color', col{1}, 'LineWidth', 2);
h(2) = plot(ax, day(keep), r_A(keep), '-', 'Color', col{2}, 'LineWidth', 2);
h(3) = plot(ax, day(keep), r_B(keep), '-', 'Color', col{3}, 'LineWidth', 2);
h(4) = plot(ax, day(keep), r_sum(keep), '--', 'Color', ink, 'LineWidth', 1.2);
set(ax, 'XScale', 'log', 'XLim', [min(day(keep)) 180], 'YLim', [0.86 1.02], ...
    'XTick', [0.5 1 2 5 10 30 60 150], 'TickDir', 'out', 'Box', 'off', ...
    'XColor', ink, 'YColor', ink, 'GridColor', grid_c, 'GridAlpha', 1, 'FontSize', 10);
grid(ax, 'on');  ax.XMinorGrid = 'off';
xlabel(ax, 'Simulated day (log scale)');
ylabel(ax, '\Sigma_i N_i  /  numeric \Sigma_i N_i');
title(ax, 'Fig. 8 offset: swap test at \gamma=5\times10^{-5}, \lambda=0.01, \delta=10, \epsilon=50', ...
    'FontWeight', 'normal', 'FontSize', 10, 'Color', ink);
legend(ax, h, {'Full MLP sNTM (MLP \phi_\tau, MLP \Delta_\tau)', ...
    'A: numeric \phi_\tau, MLP \Delta_\tau', 'B: MLP \phi_\tau, numeric \Delta_\tau', ...
    'A and B deficits added'}, 'Location', 'southoutside', 'NumColumns', 2, 'Box', 'off', 'TextColor', ink);
% direct end labels
x_end = 180;
text(ax, x_end, r_mlp(end), sprintf('  %.3f', r_mlp(end)), 'Color', ink, 'FontSize', 9);
text(ax, x_end, r_A(end), sprintf('  %.3f', r_A(end)), 'Color', ink, 'FontSize', 9);
text(ax, x_end, r_B(end), sprintf('  %.3f', r_B(end)), 'Color', ink, 'FontSize', 9);
% (end labels sit just right of the axis)

out = fullfile(here, 'outputs_fig8', 'figures');
if ~isfolder(out), mkdir(out); end
exportgraphics(fig, fullfile(out, 'swap_test.png'), 'Resolution', 200);
exportgraphics(fig, fullfile(out, 'swap_test.pdf'), 'ContentType', 'vector');
close(fig);
fprintf('wrote %s\n', fullfile(out, 'swap_test.png'));
end

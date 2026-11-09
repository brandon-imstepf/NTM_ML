function plot_fig8_gp(out_file)
% Fig. 8 of GP_manuscript.tex: regional tau of each sNTM against the numeric
% NTM at months 1-6, at the Fig. 8 parameters (gamma1 5e-5, lambda 0.01,
% delta 10, epsilon 50, seed ECl; the Table 3 timing runs, repeat 1, 180 days).
% One point per region (30). Red line: least-squares fit with intercept.
% R^2 of that fit (squared correlation) and
%   NMSE = sum_i (s_i - n_i)^2 / sum_i n_i^2 over the regions at that month.
% Concentrations in nM (N is stored in uM).
%
%   plot_fig8_gp();                 % figures/fig8_gp.png and .pdf
%   plot_fig8_gp('somewhere/fig8');

here = fileparts(mfilename('fullpath'));
nv = fileparts(here);
if nargin < 1, out_file = fullfile(here, 'figures', 'fig8_gp'); end
if ~isfolder(fileparts(out_file)), mkdir(fileparts(out_file)); end
file = @(m) fullfile(nv, 'outputs_timing', "timing__" + m + "__ECl__g5e-05_l0.01_d10_e50__rep1.mat");
methods = ["mlp", "gp", "sr", "lr"];
labels = ["MLP", "GP", "Symbolic Regression", "Linear Regression"];
months = 1:6;

num = load(file("numeric"), 'N', 'trange');
day = num.trange * 180;
fig = figure('Units', 'centimeters', 'Position', [1 1 30 22], 'Color', 'w');
tl = tiledlayout(fig, numel(methods), numel(months), 'TileSpacing', 'compact', 'Padding', 'compact');
for r = 1:numel(methods)
    s = load(file(methods(r)), 'N');
    for c = 1:numel(months)
        [~, h] = min(abs(day - 30 * months(c)));
        x = num.N(:, h) * 1e3;  y = s.N(:, h) * 1e3;
        ax = nexttile(tl);  hold(ax, 'on');
        scatter(ax, x, y, 12, [0.12 0.47 0.71], 'filled');
        p = polyfit(x, y, 1);
        xl = [min(x), max(x)];
        plot(ax, xl, polyval(p, xl), 'r-', 'LineWidth', 1);
        R = corrcoef(x, y);
        nmse = sum((y - x).^2) / sum(x.^2);
        text(ax, 0.05, 0.93, sprintf('R^2 = %.3f', R(1, 2)^2), 'Units', 'normalized', 'FontSize', 9, 'FontName', 'Times');
        text(ax, 0.95, 0.07, sprintf('NMSE = %.3g', nmse), 'Units', 'normalized', 'FontSize', 9, ...
            'FontName', 'Times', 'HorizontalAlignment', 'right');
        set(ax, 'FontSize', 8, 'FontName', 'Times', 'Box', 'on', 'TickDir', 'out');
        if r == 1, title(ax, sprintf('Month %d', months(c)), 'FontSize', 12, 'FontName', 'Times'); end
        if c == 1, ylabel(ax, labels(r), 'FontSize', 11, 'FontWeight', 'bold', 'FontName', 'Times'); end
        fprintf('%-4s month %d: R^2 %.4f, NMSE %.3g, slope %.3f\n', methods(r), months(c), R(1, 2)^2, nmse, p(1));
    end
end
xlabel(tl, '[Tau] (nM) - Numeric NTM', 'FontName', 'Times', 'FontSize', 12);
ylabel(tl, '[Tau] (nM) - sNTM', 'FontName', 'Times', 'FontSize', 12);
title(tl, 'Network wide NTM simulation regional tau accuracy', 'FontName', 'Times', 'FontSize', 15, 'FontWeight', 'bold');
exportgraphics(fig, out_file + ".png", 'Resolution', 300);
exportgraphics(fig, out_file + ".pdf", 'ContentType', 'vector');
fprintf('wrote %s.png/.pdf\n', out_file);
end

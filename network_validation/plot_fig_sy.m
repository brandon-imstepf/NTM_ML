function plot_fig_sy(out_file, run_set)
% Supplementary Fig. SY (R1.m1): network error of each sNTM against the numeric
% NTM over the 8 l2f parameter settings and the seeds in
% results/network_checks_l2f.csv (run network_checks_l2f first), or, with
% run_set = "consistent", in results/network_checks_l2f_consistent.csv.
% (a) NMSE = sum (s - n)^2 / sum n^2 over all regions and time points, one point
%     per simulation, log scale; (b) R^2 over the same points.
% Non-finite simulations are drawn at the top of (a) as crosses; in (b) the
% axis stops at R^2 = -1 (a diverged GP simulation reaches -2e3).
%
%   plot_fig_sy();                          % outputs_nm/fig_sy.png and .pdf
%   plot_fig_sy('figs/fig_sy');
%   plot_fig_sy([], "consistent");             % outputs_nm/fig_sy_consistent.png and .pdf

here = fileparts(mfilename('fullpath'));
if nargin < 2, run_set = "original"; end
[~, tag] = l2f_runs(run_set);
if nargin < 1 || isempty(out_file), out_file = fullfile(here, 'outputs_nm', "fig_sy" + tag); end
A = readtable(fullfile(here, 'results', "network_checks_l2f" + tag + ".csv"), 'TextType', 'string');
methods = ["mlp", "gp", "sr", "lr"];
labels = ["MLP", "GP", "SR", "LR"];
seeds = unique(A.seed, 'stable');
colors = [0.17 0.40 0.69; 0.85 0.45 0.10; 0.30 0.62 0.30];   % one per seed
markers = ['o', 's', '^'];
rng(0);   % horizontal jitter only

fig = figure('Units', 'centimeters', 'Position', [2 2 17 7], 'Color', 'w');
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
metrics = {'nmse', 'r2'};
ylabels = {'NMSE vs numeric NTM', 'R^2 vs numeric NTM'};
for p = 1:2
    ax = nexttile(tl);  hold(ax, 'on');
    y_all = A.(metrics{p})(A.method ~= "numeric");
    y_top = max(y_all(isfinite(y_all)));
    for k = 1:numel(methods)
        for s = 1:numel(seeds)
            y = A.(metrics{p})(A.method == methods(k) & A.seed == seeds(s));
            x = k + (s - 2) * 0.18 + 0.05 * (rand(size(y)) - 0.5);
            ok = isfinite(y);
            scatter(ax, x(ok), y(ok), 22, colors(s, :), markers(s), 'filled', 'MarkerFaceAlpha', 0.8, ...
                'DisplayName', char(seeds(s)));
            if p == 1 && any(~ok)
                plot(ax, x(~ok), repmat(y_top * 3, nnz(~ok), 1), 'x', 'Color', colors(s, :), ...
                    'MarkerSize', 7, 'LineWidth', 1.2, 'HandleVisibility', 'off');
            end
        end
        y = A.(metrics{p})(A.method == methods(k));
        plot(ax, k + [-0.35 0.35], median(y, 'omitnan') * [1 1], 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    end
    set(ax, 'XTick', 1:numel(methods), 'XTickLabel', labels, 'XLim', [0.4, numel(methods) + 0.6], ...
        'FontSize', 9, 'TickDir', 'out', 'Box', 'off');
    ax.YGrid = 'on';  ax.GridAlpha = 0.15;
    ylabel(ax, ylabels{p});
    if p == 1
        set(ax, 'YScale', 'log', 'YMinorGrid', 'off', 'YMinorTick', 'off');
        title(ax, '(a)', 'FontWeight', 'normal');
    else
        ylim(ax, [max(-1, min(A.r2(isfinite(A.r2)))) - 0.05, 1.02]);
        title(ax, '(b)', 'FontWeight', 'normal');
        legend(ax, compose('seed %s', seeds), 'Location', 'southwest', 'Box', 'off', 'FontSize', 8);
    end
end
exportgraphics(fig, out_file + ".png", 'Resolution', 300);
exportgraphics(fig, out_file + ".pdf", 'ContentType', 'vector');
fprintf('wrote %s.png/.pdf (%d simulations per method)\n', out_file, sum(A.method == "mlp"));
end

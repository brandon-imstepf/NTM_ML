function plot_shapley(models)
% Shapley-value figures from results/shapley_<model>_<target>.csv
% (shapley_surrogates.m). One figure per model (figures/shapley/shapley_<model>):
% rows phi_tau and Delta_tau; left, the mean |Shapley value| of each input as a
% share of the sum over the inputs (global importance); right, the Shapley value
% of each input at every query point, coloured by the value of that input
% (rank within the query points: low blue, high red).
% Also figures/shapley/shapley_importance_gp_vs_interp: the shares of both
% models side by side.
%
%   plot_shapley();                 % gp and interp

if nargin < 1, models = ["gp", "interp"]; end
here = fileparts(mfilename('fullpath'));
nv = fileparts(here);
fig_dir = fullfile(here, 'figures', 'shapley');
if ~isfolder(fig_dir), mkdir(fig_dir); end
font = 'Times New Roman';
inputs = ["gamma1", "lambda1", "delta", "epsilon", "N1", "N2"];
ilab = ["\gamma", "\lambda", "\delta", "\epsilon", "N_1", "N_2"];
targets = ["flux", "emc"];  tlab = ["\phi_\tau", "\Delta_\tau"];
mlab = struct('gp', "Gaussian process", 'interp', "Nearest-neighbor interpolation", ...
    'lr', "Linear regression", 'sr', "Symbolic regression");
cmap = interp1([0 0.5 1], [0.12 0.40 0.75; 0.85 0.85 0.85; 0.80 0.15 0.15], linspace(0, 1, 256));
share = struct();

for m = string(models)
    fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 24 16], 'Color', 'w');
    tl = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    for t = 1:2
        T = readtable(fullfile(nv, 'results', "shapley_" + m + "_" + targets(t) + ".csv"));
        phi = T{:, "phi_" + inputs};  X = T{:, inputs};
        imp = mean(abs(phi));  s = imp / sum(imp);  share.(m)(t, :) = s;
        [~, ord] = sort(s, 'descend');

        ax = nexttile(tl);  hold(ax, 'on');  box(ax, 'on');
        barh(ax, 1:6, 100 * s(ord), 0.6, 'FaceColor', [0.30 0.45 0.65], 'EdgeColor', 'none');
        set(ax, 'YDir', 'reverse', 'YTick', 1:6, 'YTickLabel', ilab(ord), 'FontName', font, 'FontSize', 12, ...
            'TickDir', 'out', 'XLim', [0 100 * max(s) * 1.25]);
        for k = 1:6
            text(ax, 100 * s(ord(k)) + 1, k, sprintf('%.0f%%', 100 * s(ord(k))), 'FontName', font, 'FontSize', 10);
        end
        xlabel(ax, 'Share of mean |Shapley value| (%)', 'FontName', font, 'FontSize', 11);
        ylabel(ax, tlab(t), 'FontName', font, 'FontSize', 18, 'Rotation', 0, 'HorizontalAlignment', 'right');

        ax = nexttile(tl, [1 2]);  hold(ax, 'on');  box(ax, 'on');
        sc = 10^floor(log10(max(abs(phi(:)))));
        rng(3);
        for k = 1:6
            i = ord(k);
            c = tiedrank(X(:, i)) / size(X, 1);
            y = k + 0.32 * (rand(size(phi, 1), 1) - 0.5);
            scatter(ax, phi(:, i) / sc, y, 8, c, 'filled', 'MarkerFaceAlpha', 0.7);
        end
        xline(ax, 0, '-', 'Color', [0.5 0.5 0.5]);
        colormap(ax, cmap);  clim(ax, [0 1]);
        set(ax, 'YDir', 'reverse', 'YTick', 1:6, 'YTickLabel', ilab(ord), 'YLim', [0.5 6.5], 'FontName', font, ...
            'FontSize', 12, 'TickDir', 'out');
        xlabel(ax, sprintf('Shapley value (\\times10^{%d})', log10(sc)), 'FontName', font, 'FontSize', 11);
        cb = colorbar(ax, 'Ticks', [0 1], 'TickLabels', {'low', 'high'});
        cb.Label.String = 'Input value';  cb.FontName = font;  cb.FontSize = 10;
    end
    title(tl, mlab.(m) + ": Shapley values of the inputs (" + height(T) + " held-out samples)", 'FontName', font, ...
        'FontSize', 15, 'FontWeight', 'bold');
    save_fig(fig, fullfile(fig_dir, "shapley_" + m));
end

if all(isfield(share, ["gp", "interp"]))
    fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 20 8], 'Color', 'w');
    tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    for t = 1:2
        ax = nexttile(tl);
        b = bar(ax, 100 * [share.gp(t, :); share.interp(t, :)]', 0.8, 'EdgeColor', 'none');
        b(1).FaceColor = [0.30 0.45 0.65];  b(2).FaceColor = [0.85 0.55 0.20];
        set(ax, 'XTickLabel', ilab, 'FontName', font, 'FontSize', 12, 'TickDir', 'out', 'Box', 'on');
        title(ax, tlab(t), 'FontName', font, 'FontSize', 16, 'FontWeight', 'normal');
        if t == 1, ylabel(ax, 'Share of mean |Shapley value| (%)', 'FontName', font, 'FontSize', 11); end
        if t == 2, legend(ax, ["Gaussian process", "Nearest-neighbor interpolation"], 'FontName', font, 'FontSize', 10, 'Location', 'northeast'); end
    end
    save_fig(fig, fullfile(fig_dir, "shapley_importance_gp_vs_interp"));
end
end

function save_fig(fig, out)
exportgraphics(fig, out + ".png", 'Resolution', 300, 'BackgroundColor', 'w');
exportgraphics(fig, out + ".pdf", 'ContentType', 'vector', 'BackgroundColor', 'w');
close(fig);
fprintf('wrote %s\n', out);
end

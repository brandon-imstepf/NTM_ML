function plot_results_variants(which)
% Candidate "Results" figures (main text) for the GP and for griddatan, each
% surrogate on its own (no comparison with the other surrogates):
%   gp_results_v1   = plot_gp_results_fig.m (training and held-out NMSE vs n)
%   gp_results_v2   held-out predicted vs true (density), phi_tau and Delta_tau
%   gp_results_v3   v1 (top) and v2 (bottom) together
%   gp_results_v4   held-out error across the input space: MSE in each decile of
%                   each input, divided by the variance of the target over the
%                   whole held-out set (so that the bins are comparable)
%   gd_results_v1   griddatan 'nearest': held-out NMSE vs n for three input
%                   scalings (griddatan_size_curve.m); the network model in red
%   gd_results_v2   as gp_results_v2, for the network model (log gamma, lambda; 10^6)
%   gd_results_v3   gd v1 (top) and v2 (bottom) together
%   gd_results_v4   as gp_results_v4, for the network model
% Held out: data_bias_e5_11.csv (10^5; never used for fitting or selection).
% GP: the paper's GP (n = 5000). griddatan network model: nearest neighbour over
% 10^6 rows on log10 gamma1, log10 lambda1 and the other inputs, all z-scored
% (the "interp" sNTM, make_surrogate.m).
% Predictions on the held-out set are cached in figures/cache/heldout_pred.mat.
%
%   plot_results_variants();                 % all
%   plot_results_variants("gp");             % GP only

if nargin < 1, which = ["gp", "gd"]; end
here = fileparts(mfilename('fullpath'));
nv = fileparts(here);
addpath(nv);
fig_dir = fullfile(here, 'figures', 'results_variants');
if ~isfolder(fig_dir), mkdir(fig_dir); end
P = heldout_predictions(here, nv);
font = 'Times New Roman';
tlab = ["\phi_\tau", "\Delta_\tau"];

for w = string(which)
    if w == "gp"
        pred = P.gp;  name = "Gaussian process";  stem = "gp_results";
        if ~isfile(fullfile(fig_dir, stem + "_v1.png"))
            plot_gp_results_fig(fullfile(fig_dir, stem + "_v1"));
        end
    else
        pred = P.interp;  name = "Nearest-neighbor interpolation";  stem = "gd_results";
        fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 20 9], 'Color', 'w');
        tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
        gd_size_panels(tl, nv, font, tlab);
        title(tl, 'Nearest-neighbor interpolation training set size', 'FontName', font, 'FontSize', 18, 'FontWeight', 'bold');
        save_fig(fig, fullfile(fig_dir, stem + "_v1"));
    end

    % v2: parity
    fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 20 9.5], 'Color', 'w');
    tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    parity_panels(tl, P.y, pred, font, tlab);
    title(tl, name + ": held-out predictions", 'FontName', font, 'FontSize', 18, 'FontWeight', 'bold');
    save_fig(fig, fullfile(fig_dir, stem + "_v2"));

    % v3: size curve + parity
    fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 20 18], 'Color', 'w');
    tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    if w == "gp", gp_size_panels(tl, nv, font, tlab); else, gd_size_panels(tl, nv, font, tlab); end
    parity_panels(tl, P.y, pred, font, tlab);
    title(tl, name, 'FontName', font, 'FontSize', 18, 'FontWeight', 'bold');
    for pl = [struct('s', 'a', 'y', 0.90), struct('s', 'b', 'y', 0.44)]     % panel letters
        annotation(fig, 'textbox', [0.005 pl.y 0.05 0.05], 'String', pl.s, 'FontName', font, 'FontSize', 20, ...
            'FontWeight', 'bold', 'EdgeColor', 'none');
    end
    save_fig(fig, fullfile(fig_dir, stem + "_v3"));

    % v4: error across the input space
    fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 24 13], 'Color', 'w');
    tl = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    decile_panels(tl, P.X, P.y, pred, font, tlab);
    title(tl, name + ": held-out error across the input space", 'FontName', font, 'FontSize', 16, 'FontWeight', 'bold');
    save_fig(fig, fullfile(fig_dir, stem + "_v4"));
end
end

function P = heldout_predictions(here, nv)
cache = fullfile(here, 'figures', 'cache', 'heldout_pred.mat');
if isfile(cache), P = load(cache); return; end
cfg = sntm_config();
H = readmatrix(fullfile(cfg.ntm_ml_root, 'data', 'data_neural_network_model', 'training_data_e6', 'data_bias_e5_11.csv'));
P.X = H(:, 1:6);  P.y = H(:, 7:8);
G = load(cfg.gp_model_file);
P.gp = [predict(G.stzd_flux_gprModel, P.X) * G.flux_std + G.flux_mean, ...
        predict(G.stzd_emc_gprModel, P.X) * G.emc_std + G.emc_mean];
S = make_surrogate('interp', cfg);
P.interp = [S.predict_flux(P.X), S.predict_w1(P.X)];
if ~isfolder(fileparts(cache)), mkdir(fileparts(cache)); end
save(cache, '-struct', 'P');
end

function gp_size_panels(tl, nv, font, tlab)
T = readtable(fullfile(nv, 'results', 'supp_gp_griddatan.csv'), 'TextType', 'string');
gp = sortrows(T(startsWith(T.model, "gp_") & T.model ~= "gp_lrsr_1000", :), 'n_train');
target = ["flux", "emc"];
c_train = [0 0.447 0.741];  c_test = [0.466 0.674 0.188];  c_sel = [1 0 0];
for t = 1:2
    ax = nexttile(tl);  hold(ax, 'on');  box(ax, 'on');
    h1 = errorbar(ax, gp.n_train, gp.("train_" + target(t) + "_nmse_mean"), gp.("train_" + target(t) + "_nmse_sd"), ...
        '-o', 'Color', c_train, 'MarkerFaceColor', c_train, 'MarkerSize', 4, 'LineWidth', 1.4, 'CapSize', 3);
    h2 = errorbar(ax, gp.n_train, gp.("heldout_" + target(t) + "_nmse_mean"), gp.("heldout_" + target(t) + "_nmse_sd"), ...
        '-o', 'Color', c_test, 'MarkerFaceColor', c_test, 'MarkerSize', 4, 'LineWidth', 1.4, 'CapSize', 3);
    k = gp.n_train == 5000;
    h3 = plot(ax, [gp.n_train(k) gp.n_train(k)], [gp.("train_" + target(t) + "_nmse_mean")(k), ...
        gp.("heldout_" + target(t) + "_nmse_mean")(k)], 'o', 'MarkerEdgeColor', c_sel, 'MarkerFaceColor', c_sel, 'MarkerSize', 6);
    set(ax, 'YScale', 'log', 'XLim', [0 5400], 'YLim', [1e-8 1], 'XTick', 0:1000:5000, 'XTickLabelRotation', 0, ...
        'XTickLabel', ["0", "", "2000", "", "4000", ""], 'YTick', 10.^(-8:2:0), 'FontName', font, 'FontSize', 12, 'TickDir', 'out');
    title(ax, tlab(t), 'FontName', font, 'FontSize', 16, 'FontWeight', 'normal');
    xlabel(ax, 'Number of training samples', 'FontName', font, 'FontSize', 12);
    if t == 1, ylabel(ax, 'NMSE', 'FontName', font, 'FontSize', 13); else, ax.YTickLabel = []; end
    if t == 2
        legend(ax, [h1 h2 h3(1)], {'Training', 'Held-out test', 'Selected model'}, 'Location', 'east', ...
            'FontName', font, 'FontSize', 10);
    end
end
end

function gd_size_panels(tl, nv, font, tlab)
T = readtable(fullfile(nv, 'results', 'griddatan_size_curve.csv'), 'TextType', 'string');
sc = ["raw", "zscore", "logz"];
lab = ["Raw inputs", "Z-scored inputs", "log_{10} \gamma, \lambda, then z-scored"];
col = [0.55 0.55 0.55; 0 0.447 0.741; 0.466 0.674 0.188];
target = ["flux", "emc"];
for t = 1:2
    ax = nexttile(tl);  hold(ax, 'on');  box(ax, 'on');
    h = gobjects(1, 4);
    for s = 1:3
        r = sortrows(T(T.scaling == sc(s), :), 'n_train');
        h(s) = errorbar(ax, r.n_train, r.("heldout_" + target(t) + "_nmse_mean"), r.("heldout_" + target(t) + "_nmse_sd"), ...
            '-o', 'Color', col(s, :), 'MarkerFaceColor', col(s, :), 'MarkerSize', 4, 'LineWidth', 1.4, 'CapSize', 3);
    end
    k = T.scaling == "logz" & T.n_train == 1e6;
    h(4) = plot(ax, 1e6, T.("heldout_" + target(t) + "_nmse_mean")(k), 'o', 'MarkerEdgeColor', [1 0 0], ...
        'MarkerFaceColor', [1 0 0], 'MarkerSize', 6);
    set(ax, 'XScale', 'log', 'YScale', 'log', 'XLim', [350 1.4e6], 'YLim', [0.015 40], 'XTick', 10.^(3:6), 'YTick', 10.^(-1:1), ...
        'FontName', font, 'FontSize', 12, 'TickDir', 'out');
    title(ax, tlab(t), 'FontName', font, 'FontSize', 16, 'FontWeight', 'normal');
    xlabel(ax, 'Number of training samples', 'FontName', font, 'FontSize', 12);
    if t == 1, ylabel(ax, 'Held-out NMSE', 'FontName', font, 'FontSize', 13); else, ax.YTickLabel = []; end
    if t == 2
        legend(ax, h, [lab, "Network model"], 'Location', 'north', 'NumColumns', 2, 'FontName', font, 'FontSize', 9);
    end
end
end

function parity_panels(tl, y, p, font, tlab)
for t = 1:2
    ax = nexttile(tl);  hold(ax, 'on');  box(ax, 'on');
    sc = 10^floor(log10(max(abs(y(:, t)))));        % axis unit
    a = y(:, t) / sc;  b = p(:, t) / sc;
    lim = reshape(prctile([a; b], [0.1 99.9]), 1, 2);  lim = lim + [-1 1] * 0.03 * diff(lim);
    e = linspace(lim(1), lim(2), 121);
    C = histcounts2(b, a, e, e);                    % rows = predicted
    C(C == 0) = NaN;
    imagesc(ax, e, e, log10(C), 'AlphaData', ~isnan(C));
    colormap(ax, flipud(gray(256)) * 0.85);
    plot(ax, lim, lim, '--', 'Color', [0.85 0.1 0.1], 'LineWidth', 1);
    axis(ax, 'square');  xlim(ax, lim);  ylim(ax, lim);  set(ax, 'YDir', 'normal');
    set(ax, 'FontName', font, 'FontSize', 11, 'TickDir', 'out', 'Layer', 'top');
    title(ax, tlab(t), 'FontName', font, 'FontSize', 16, 'FontWeight', 'normal');
    subtitle(ax, sprintf('NMSE = %.3g', nmse_eq7(p(:, t), y(:, t))), 'FontName', font, 'FontSize', 10);
    xlabel(ax, sprintf('Numeric (\\times10^{%d})', log10(sc)), 'FontName', font, 'FontSize', 12);
    ylabel(ax, sprintf('Predicted (\\times10^{%d})', log10(sc)), 'FontName', font, 'FontSize', 12);
    cb = colorbar(ax);  cb.Label.String = 'log_{10} count';  cb.FontName = font;
end
end

function decile_panels(tl, X, y, p, font, tlab)
names = ["\gamma", "\lambda", "\delta", "\epsilon", "N_1 (\muM)", "N_2 (\muM)"];
logx = [true true false false false false];
col = [0 0.447 0.741; 0.85 0.33 0.10];
v = var(y);
for i = 1:6
    ax = nexttile(tl);  hold(ax, 'on');  box(ax, 'on');
    ed = prctile(X(:, i), 0:10:100);  ed(end) = ed(end) + eps(ed(end));
    [~, ~, bin] = histcounts(X(:, i), ed);
    xm = accumarray(bin, X(:, i), [10 1], @median);
    h = gobjects(1, 2);
    for t = 1:2
        m = accumarray(bin, (p(:, t) - y(:, t)).^2, [10 1], @mean) / v(t);
        h(t) = plot(ax, xm, m, '-o', 'Color', col(t, :), 'MarkerFaceColor', col(t, :), 'MarkerSize', 4, 'LineWidth', 1.3);
    end
    set(ax, 'YScale', 'log', 'FontName', font, 'FontSize', 10, 'TickDir', 'out');
    if logx(i), set(ax, 'XScale', 'log'); end
    xlabel(ax, names(i), 'FontName', font, 'FontSize', 13);
    if mod(i, 3) == 1, ylabel(ax, 'MSE / Var (held-out)', 'FontName', font, 'FontSize', 11); end
    if i == 1, legend(ax, h, tlab, 'Location', 'best', 'FontName', font, 'FontSize', 11); end
end
end

function save_fig(fig, out)
% all text 25% larger than set above (requested 11/13 for the paper's Figs. 5 and 6)
% (sizes recorded first and set from them: axes before their labels, which
% would otherwise follow the axes and be enlarged twice)
h = findall(fig, '-property', 'FontSize');
f0 = arrayfun(@(x) x.FontSize, h);
isax = arrayfun(@(x) isa(x, 'matlab.graphics.axis.Axes'), h);
for k = [find(isax); find(~isax)]'
    h(k).FontSize = 1.25 * f0(k);
end
exportgraphics(fig, out + ".png", 'Resolution', 300, 'BackgroundColor', 'w');
exportgraphics(fig, out + ".pdf", 'ContentType', 'vector', 'BackgroundColor', 'w');
close(fig);
fprintf('wrote %s\n', out);
end

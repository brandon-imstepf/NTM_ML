function plot_sr_lr_results(which)
% Results figures for symbolic regression (Fig. 4) and linear regression
% (Fig. 3), remade (11/19, Suzanne) in the style of the 1-NN and GP figures
% (plot_results_variants.m, v5): the selected model as red dots instead of a
% dashed line, and a panel (b) of held-out predictions against the numeric values
% (density), on the same held-out set as the GP and 1-NN figures
% (data_bias_e5_11.csv, 10^5 samples, never used for fitting or selection).
%
%   sr_results: (a) training and validation NMSE along the Pareto front of the
%               PySR run that gave each selected expression (model_equations.csv
%               and metrics_summary.csv of that run; "validation" is PySR's own
%               validation split); (b) held-out predictions of the selected
%               expressions (make_surrogate 'sr'), with the expression; (c) how
%               many of the 11 runs' "best" expressions (best.txt) contain
%               each input
%   lr_results: (a) forward feature selection: training and 10-fold
%               cross-validation NMSE, mean and SD over the folds, of the
%               feature added at each step (results/lr_selection_curve.csv,
%               lr_selection_curve.py); (b) held-out predictions of the selected
%               10-feature models (make_surrogate 'lr')
%
%   plot_sr_lr_results()  or  plot_sr_lr_results("sr")
% Output: figures/results_variants/sr_results.png/.pdf, lr_results.png/.pdf

if nargin < 1, which = ["sr", "lr"]; end
here = fileparts(mfilename('fullpath'));
nv = fileparts(here);
addpath(nv);
fig_dir = fullfile(here, 'figures', 'results_variants');
if ~isfolder(fig_dir), mkdir(fig_dir); end
P = heldout_predictions(here);
font = 'Times New Roman';
tlab = ["\phi_\tau", "\Delta_\tau"];

for w = string(which)
    if w == "sr"
        % (11/19, user) four panels in a 2 x 2 arrangement, laid out by hand
        % in cm: (a) Pareto front (upper left), (b) expressions (upper right),
        % (c) appearance frequency (lower left), (d) held-out predictions
        % (lower right); each of a, c and d holds a phi_tau and a Delta_tau axes
        % (11/19, user) no title; expressions left-aligned; (d) one y axis
        % (phi_tau) and one colour bar (Delta_tau), shared colour scale
        W = 29.6;  H = 18.4;  xL = [2.0 8.6];  xR = [16.6 22.5];  w = 5.6;  wd = 5.4;
        fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 W H], 'Color', 'w');
        cm = @(x, y, w, h) axes(fig, 'Units', 'centimeters', 'Position', [x y w h]);
        yt = 11.8;  ht = 5.4;  yb = 2.4;       % bottom and height of the top row; bottom of the lower row
        sr_pareto_panels(@(t) cm(xL(t), yt, w, ht), font, tlab);
        sr_equation_panel(cm(xR(1), yt, xR(2) + wd + 1.5 - xR(1), ht));   % over (d): centred on its axes
        sr_freq_panels(@(t) cm(xL(t), yb, w, wd), font, tlab);
        % the colour bar takes its width from the Delta_tau axes: 1.5 cm extra there
        % keeps both square plot boxes wd wide
        parity_panels(@(t) cm(xR(t), yb, wd + 1.5 * (t == 2), wd), P.y, P.sr, font, tlab, strings(1, 2), true);
        ov = cm(0, 0, W, H);  set(ov, 'Visible', 'off', 'XLim', [0 W], 'YLim', [0 H]);
        for L = [struct('s', 'a', 'x', 0.1, 'y', yt + ht + 1.2), struct('s', 'b', 'x', 14.8, 'y', yt + ht + 1.2), ...
                 struct('s', 'c', 'x', 0.1, 'y', yb + wd + 2.0), struct('s', 'd', 'x', 14.8, 'y', yb + wd + 2.0)]
            text(ov, L.x, L.y, L.s, 'VerticalAlignment', 'top', 'FontName', font, 'FontSize', 20, 'FontWeight', 'bold');
        end
        save_fig(fig, fullfile(fig_dir, "sr_results"));
    else
        % (11/19, user) as (d) of the SR figure: (b) with one y axis (phi_tau)
        % and one colour bar (Delta_tau) on a shared colour scale; laid out by
        % hand in cm, no title (as the SR figure)
        W = 17.6;  H = 19.2;  x = [2.0 8.7];  w = 6.2;  wd = 6.2;
        fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 W H], 'Color', 'w');
        cm = @(x, y, w, h) axes(fig, 'Units', 'centimeters', 'Position', [x y w h]);
        yt = 12.2;  ht = 5.6;  yb = 2.4;       % bottom and height of (a); bottom of (b)
        lr_selection_panels(@(t) cm(x(t), yt, w, ht), nv, font, tlab);
        % the colour bar takes its width from the Delta_tau axes: 1.5 cm extra there
        parity_panels(@(t) cm(x(t), yb, wd + 1.5 * (t == 2), wd), P.y, P.lr, font, tlab, strings(1, 2), true);
        ov = cm(0, 0, W, H);  set(ov, 'Visible', 'off', 'XLim', [0 W], 'YLim', [0 H]);
        for L = [struct('s', 'a', 'y', yt + ht + 1.2), struct('s', 'b', 'y', yb + wd + 2.0)]
            text(ov, 0.1, L.y, L.s, 'VerticalAlignment', 'top', 'FontName', font, 'FontSize', 20, 'FontWeight', 'bold');
        end
        save_fig(fig, fullfile(fig_dir, "lr_results"));
    end
end
end

function ax1 = sr_pareto_panels(mk, font, tlab)
% the Pareto front of the run that gave each selected expression; the selected
% expression in red (complexity 41 for phi_tau, from best.txt; 25 for Delta_tau,
% chosen by hand in final_figure_v4.py, W1_COMPLEXITY_OVERRIDE)
root = 'C:\Users\USER\Documents\Code\UCM_AD_ML\datasets\sim_csv_v11\train';
run = {fullfile(root, 'flux_only_bias_train_e3', '2025-10-21_09-36_259200_Info'), ...
       fullfile(root, 'w1_only_bias_train_e3', '2025-10-21_09-14_259200_Info')};
sel = [41 25];
c_train = [0 0.447 0.741];  c_val = [0.466 0.674 0.188];  c_sel = [1 0 0];
for t = 1:2
    E = readtable(fullfile(run{t}, 'model_equations.csv'));
    M = readtable(fullfile(run{t}, 'metrics_summary.csv'));
    [cx, o] = sort(E.complexity);  tr = M.NMSE_Train(o);  va = M.NMSE_Validation(o);
    ax = mk(t);  hold(ax, 'on');  box(ax, 'on');
    if t == 1, ax1 = ax; end
    h1 = plot(ax, cx, tr, '-o', 'Color', c_train, 'MarkerFaceColor', c_train, 'MarkerSize', 4, 'LineWidth', 1.4);
    h2 = plot(ax, cx, va, '-o', 'Color', c_val, 'MarkerFaceColor', c_val, 'MarkerSize', 4, 'LineWidth', 1.4);
    k = cx == sel(t);
    assert(nnz(k) == 1, 'selected complexity %d not on the Pareto front', sel(t));
    h3 = plot(ax, [sel(t) sel(t)], [tr(k) va(k)], 'o', 'MarkerEdgeColor', c_sel, 'MarkerFaceColor', c_sel, 'MarkerSize', 6);
    fprintf('SR %s: complexity %d, training NMSE %.3f, validation NMSE %.3f\n', tlab(t), sel(t), tr(k), va(k));
    set(ax, 'YScale', 'log', 'XLim', [0 65], 'YLim', [0.1 2], 'YTick', [0.1 0.2 0.5 1 2], ...
        'FontName', font, 'FontSize', 12, 'TickDir', 'out');
    title(ax, tlab(t), 'FontName', font, 'FontSize', 16, 'FontWeight', 'normal');
    xlabel(ax, 'Complexity', 'FontName', font, 'FontSize', 12);
    if t == 1, ylabel(ax, 'NMSE', 'FontName', font, 'FontSize', 13); else, ax.YTickLabel = []; end
    if t == 1
        legend(ax, [h1 h2 h3(1)], {'Training', 'Validation', 'Selected model'}, 'Location', 'northeast', ...
            'FontName', font, 'FontSize', 10);
    end
end
end

function ax1 = sr_freq_panels(mk, font, tlab)
% number of the 11 runs' "best" expressions (best.txt) that contain each input
root = 'C:\Users\USER\Documents\Code\UCM_AD_ML\datasets\sim_csv_v11\train';
dirs = ["flux_only_bias_train_e3", "w1_only_bias_train_e3"];
names = ["\gamma", "\lambda", "\delta", "\epsilon", "N_1", "N_2"];
for t = 1:2
    d = dir(fullfile(root, dirs(t), '*_Info', 'best.txt'));
    cnt = zeros(1, 6);
    for f = 1:numel(d)
        s = fileread(fullfile(d(f).folder, d(f).name));
        s = regexp(s, 'Selected Best Function:\s*\d+,\s*([^\r\n]*)', 'tokens', 'once');
        % a term x_k/x_k is 1: not counted as containing x_k (as final_figure_v4.py,
        % which counts after simplification)
        s{1} = regexprep(s{1}, 'x(\d)/x\1(?!\d)', '1');
        for v = 0:5, cnt(v + 1) = cnt(v + 1) + ~isempty(regexp(s{1}, sprintf('x%d(?!\\d)', v), 'once')); end
    end
    fprintf('SR %s: %d runs; inputs present %s\n', tlab(t), numel(d), mat2str(cnt));
    ax = mk(t);  hold(ax, 'on');  box(ax, 'on');
    if t == 1, ax1 = ax; end
    bar(ax, 1:6, cnt, 0.6, 'FaceColor', [0.55 0.55 0.55], 'EdgeColor', 'none');
    for v = 1:6
        text(ax, v, cnt(v) + 0.25, sprintf('%d', cnt(v)), 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', 'FontName', font, 'FontSize', 11);
    end
    set(ax, 'XLim', [0.4 6.6], 'YLim', [0 numel(d) + 2.5], 'XTick', 1:6, 'XTickLabel', names, 'XTickLabelRotation', 0, ...
        'YTick', 0:4:numel(d) + 1, 'FontName', font, 'FontSize', 12, 'TickDir', 'out');
    title(ax, tlab(t), 'FontName', font, 'FontSize', 16, 'FontWeight', 'normal');
    if t == 1, ylabel(ax, 'Appearance frequency', 'FontName', font, 'FontSize', 13);
    else, ax.YTickLabel = []; end
end
end

function lr_selection_panels(mk, nv, font, tlab)
T = readtable(fullfile(nv, 'results', 'lr_selection_curve.csv'), 'TextType', 'string');
target = ["flux", "emc"];
c_train = [0 0.447 0.741];  c_cv = [0.466 0.674 0.188];  c_sel = [1 0 0];
nsel = 10;
for t = 1:2
    r = sortrows(T(T.target == target(t), :), 'n_features');
    ax = mk(t);  hold(ax, 'on');  box(ax, 'on');
    h1 = errorbar(ax, r.n_features, r.train_nmse_mean, r.train_nmse_sd, '-o', 'Color', c_train, ...
        'MarkerFaceColor', c_train, 'MarkerSize', 4, 'LineWidth', 1.4, 'CapSize', 3);
    h2 = errorbar(ax, r.n_features, r.cv_nmse_mean, r.cv_nmse_sd, '-o', 'Color', c_cv, ...
        'MarkerFaceColor', c_cv, 'MarkerSize', 4, 'LineWidth', 1.4, 'CapSize', 3);
    k = r.n_features == nsel;
    h3 = plot(ax, [nsel nsel], [r.train_nmse_mean(k) r.cv_nmse_mean(k)], 'o', 'MarkerEdgeColor', c_sel, ...
        'MarkerFaceColor', c_sel, 'MarkerSize', 6);
    fprintf('LR %s: %d features, training NMSE %.3f, CV NMSE %.3f\n', tlab(t), nsel, r.train_nmse_mean(k), r.cv_nmse_mean(k));
    set(ax, 'XLim', [0 28], 'YLim', [0 1], 'XTick', 0:5:25, 'YTick', 0:0.2:1, ...
        'FontName', font, 'FontSize', 12, 'TickDir', 'out');
    title(ax, tlab(t), 'FontName', font, 'FontSize', 16, 'FontWeight', 'normal');
    xlabel(ax, 'Number of features', 'FontName', font, 'FontSize', 12);
    if t == 1, ylabel(ax, 'NMSE', 'FontName', font, 'FontSize', 13); else, ax.YTickLabel = []; end
    if t == 2
        legend(ax, [h1 h2 h3(1)], {'Training', 'Cross-validation', 'Selected model'}, 'Location', 'northeast', ...
            'FontName', font, 'FontSize', 10);
    end
end
end

function P = heldout_predictions(here)
% the held-out set and the GP/1-NN predictions come from plot_results_variants'
% cache; the SR and LR predictions are cached here
cache = fullfile(here, 'figures', 'cache', 'heldout_pred_lrsr.mat');
if isfile(cache), P = load(cache); return; end
Q = load(fullfile(here, 'figures', 'cache', 'heldout_pred.mat'), 'X', 'y');
cfg = sntm_config();
P.X = Q.X;  P.y = Q.y;
for m = ["sr", "lr"]
    S = make_surrogate(char(m), cfg);
    P.(m) = [S.predict_flux(P.X), S.predict_w1(P.X)];
end
save(cache, '-struct', 'P');
end

function ax1 = parity_panels(mk, y, p, font, tlab, eq, compact)
% as plot_results_variants; eq (LaTeX) is written under the title when given;
% compact (11/19, SR figure): y label and tick labels on the phi_tau panel only
% (each panel is square with equal x and y ranges, so the x ticks give the
% scale) and one colour bar, on the Delta_tau panel, with a colour scale shared
% by both panels
if nargin < 7, compact = false; end
axs = gobjects(1, 2);  cmax = 0;
for t = 1:2
    ax = mk(t);  hold(ax, 'on');  box(ax, 'on');
    if t == 1, ax1 = ax; end
    sc = 10^floor(log10(max(abs(y(:, t)))));        % axis unit
    a = y(:, t) / sc;  b = p(:, t) / sc;
    lim = reshape(prctile([a; b], [0.1 99.9]), 1, 2);  lim = lim + [-1 1] * 0.03 * diff(lim);
    e = linspace(lim(1), lim(2), 121);
    C = histcounts2(b, a, e, e);                    % rows = predicted
    C(C == 0) = NaN;
    imagesc(ax, e, e, log10(C), 'AlphaData', ~isnan(C));
    axs(t) = ax;  cmax = max(cmax, max(log10(C(:)), [], 'omitnan'));
    colormap(ax, flipud(gray(256)) * 0.85);
    plot(ax, lim, lim, '--', 'Color', [0.85 0.1 0.1], 'LineWidth', 1);
    axis(ax, 'square');  xlim(ax, lim);  ylim(ax, lim);  set(ax, 'YDir', 'normal');
    set(ax, 'FontName', font, 'FontSize', 11, 'TickDir', 'out', 'Layer', 'top');
    title(ax, tlab(t), 'FontName', font, 'FontSize', 16, 'FontWeight', 'normal');
    n = nmse_eq7(p(:, t), y(:, t));
    if eq(t) == ""
        subtitle(ax, sprintf('NMSE = %.3g', n), 'FontName', font, 'FontSize', 10);
    else
        subtitle(ax, {char(eq(t)), sprintf('NMSE = %.3g', n)}, 'Interpreter', 'latex', 'FontSize', 10);
    end
    fprintf('held-out %s: NMSE %.4g\n', tlab(t), n);
    xlabel(ax, sprintf('Numeric (\\times10^{%d})', log10(sc)), 'FontName', font, 'FontSize', 12);
    if compact && t == 2
        ax.YTickLabel = [];
    else
        ylabel(ax, sprintf('Predicted (\\times10^{%d})', log10(sc)), 'FontName', font, 'FontSize', 12);
    end
    if ~compact || t == 2
        cb = colorbar(ax);  cb.Label.String = 'log_{10} count';  cb.FontName = font;
    end
end
if compact, set(axs, 'CLim', [0 cmax]); end
end

function ax = sr_equation_panel(ax)
% the two selected expressions, each over its column (inputs z-scored as in fitting)
eq = ["$\phi_\tau = -0.487N_2 - 0.812 + \frac{(N_1 + 3.38)^{-2.64\gamma}}{4.96 - \delta}$", ...
      "$\Delta_\tau = \frac{0.0453\delta - \gamma + 0.0453\lambda - 0.509}{\gamma + 0.691}$"];
set(ax, 'XLim', [0 1], 'YLim', [0 1], 'Visible', 'off');
for t = 1:2      % one under the other, centred over (d) (11/19)
    text(ax, 0.5, 0.78 - 0.56 * (t - 1), eq(t), 'Interpreter', 'latex', 'FontSize', 15, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end
end

function panel_letters(fig, font, s, axs)
% letters at the left edge of the figure, level with the top of each panel's
% outer box (titles included)
drawnow;
for k = 1:numel(s)
    p = axs{k}.OuterPosition;
    annotation(fig, 'textbox', [0.005, p(2) + p(4) - 0.035, 0.05, 0.035], 'String', s(k), 'FontName', font, ...
        'FontSize', 20, 'FontWeight', 'bold', 'EdgeColor', 'none', 'VerticalAlignment', 'top');
end
end

function letters(fig, font, s, y)
for k = 1:numel(s)
    annotation(fig, 'textbox', [0.005 y(k) 0.05 0.05], 'String', s(k), 'FontName', font, 'FontSize', 20, ...
        'FontWeight', 'bold', 'EdgeColor', 'none');
end
end

function save_fig(fig, out)
% all text 25% larger than set above, as plot_results_variants (11/13)
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

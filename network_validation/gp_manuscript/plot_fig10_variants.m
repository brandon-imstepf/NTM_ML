function plot_fig10_variants(variants)
% Layout experiments for Fig. 10 (11/19, user): the grid of plot_fig8_gp
% flipped, months in rows and the models in columns (LR, SR, 1-NN, GP, MLP
% from left to right), at the Primary Parameter Set. Same data, statistics and
% font sizes as plot_fig8_gp (fit through the origin y = b x, R^2 = squared
% Pearson correlation, NMSE of Eq. 7, concentrations in nM). Variants:
%   A  text width (13.1 cm), square panels, b / R^2 / NMSE on three labelled
%      lines above each panel
%   B  landscape page (19.4 x 13.1 cm, for a sideways figure), key
%      (b, R^2, NMSE) under the title and one line of numbers above each panel
%   C  text width, square panels without numbers; the numbers in a table under
%      the grid, one column per model
%   D  text width, one row: the three months overlaid in one panel per model
%      (log-log axes, one colour per month), numbers in a table as in C
%
%   plot_fig10_variants()            all four
%   plot_fig10_variants(["A" "C"])
% Output: figures/fig10_variants/fig10_<variant>.png/.pdf

here = fileparts(mfilename('fullpath'));
nv = fileparts(here);
addpath(nv);
if nargin < 1, variants = ["A", "B", "C", "D"]; end
out_dir = fullfile(here, 'figures', 'fig10_variants');
if ~isfolder(out_dir), mkdir(out_dir); end
src_fmt = 'outputs_consistent/rec_l2f__%s__CA1__g0.005_l0.02_d20_e80.mat';
file = @(m) fullfile(nv, sprintf(src_fmt, m));

% data and statistics, shared by every variant
D.methods = ["lr", "sr", "interp", "gp", "mlp"];      % left to right (model order)
D.labels = ["LR", "SR", "1-NN", "GP", "MLP"];
D.months = [1 3 5];
D.title = 'Network wide NTM simulation regional tau accuracy';
num = load(file("numeric"), 'N', 'trange');
day = num.trange * 180;
hm = arrayfun(@(m) find(abs(day - 30 * m) == min(abs(day - 30 * m)), 1), D.months);
nm = numel(D.methods);  nmo = numel(D.months);
nmax = max(num.N(:, hm), [], 'all');  vmax = nmax;
[D.b, D.r2, D.nmse] = deal(nan(nmo, nm));
D.x = cell(1, nmo);  D.y = cell(nmo, nm);
for c = 1:nmo, D.x{c} = num.N(:, hm(c)) * 1e3; end
for k = 1:nm
    s = load(file(D.methods(k)), 'N');
    for c = 1:nmo
        x = D.x{c};  y = s.N(:, hm(c)) * 1e3;  D.y{c, k} = y;
        v = y(isfinite(y) & y <= 3 * nmax * 1e3);
        if ~isempty(v), vmax = max(vmax, max(v) / 1e3); end
        if ~all(isfinite(y)), continue, end
        D.b(c, k) = sum(x .* y) / sum(x .^ 2);
        Rc = corrcoef(x, y);  D.r2(c, k) = Rc(1, 2)^2;
        D.nmse(c, k) = nmse_eq7(y, x);
    end
end
% linear range and ticks as in plot_fig8_gp
D.lim = [0, 1.05 * vmax * 1e3];
st = 10 .^ floor(log10(D.lim(2))) * [0.1 0.2 0.25 0.4 0.5 1 2 2.5 4 5];
D.tk = 0:st(find(D.lim(2) ./ st <= 3, 1)):D.lim(2);
% log range (variant D): from the decade below the smallest positive numeric
% concentration to the decade above the linear range
xp = cell2mat(D.x(:));  xp = xp(xp > 0);
D.llim = 10 .^ [floor(log10(min(xp))), ceil(log10(D.lim(2)))];
D.llim = [min(D.llim(1), 0.01), max(D.llim(2), 1)];

for v = string(variants)
    draw_variant(v, D, fullfile(out_dir, "fig10_" + v));
end
end

function draw_variant(v, D, out_file)
% font sizes in pt as printed (as plot_fig8_gp)
f_tick = 9.5;  f_box = 10.5;  f_key = 12;  f_row = 13;  f_col = 13.5;  f_lab = 13.5;  f_leg = 12;  f_title = 13;
nc = numel(D.methods);  nmo = numel(D.months);
L0 = 2.0;  R = 0.15;  gx = 0.3;  gy = 0.25;  key = false;  tbl = false;  overlay = false;
switch v
    case "A", W = 13.1;  nr = nmo;  Sb = 1.6;  mode = "stack";  R = 0.45;  gy = 0.35;
    case "B", W = 19.4;  nr = nmo;  Sb = 0.5;  mode = "line";  key = true;  gx = 0.35;  gy = 0.3;  Hmax = 13.1;
    case "C", W = 13.1;  nr = nmo;  Sb = 0.1;  mode = "none";  tbl = true;
    case "D", W = 13.1;  nr = 1;  Sb = 0.1;  mode = "none";  tbl = true;  overlay = true;  gx = 0.55;
end
w = (W - L0 - R - (nc - 1) * gx) / nc;
Tt = 0.75 + 0.65 * key + 0.65;                    % title, key, model names
lineH = 0.45;                                      % table line height
tblH = tbl * ((nmo * 3 + 1) * lineH + (nmo - 1) * 0.2 + 0.4);   % with a header row of model names
legH = 0.75 + 0.55 * overlay;                      % legend (two lines in D)
Bb = 1.25 + legH + tblH;                           % tick labels, x label, legend (11/19: under the panels), table
if v == "B"
    a = (Hmax - Bb - Tt - nr * Sb - (nr - 1) * gy) / nr;   % stretched in x
else
    a = w;                                         % square
end
H = Bb + nr * (a + Sb) + (nr - 1) * gy + Tt;
col = month_colours();

fig = figure('Units', 'centimeters', 'Position', [1 1 W H], 'Color', 'w', ...
    'PaperUnits', 'centimeters', 'PaperSize', [W H], 'PaperPosition', [0 0 W H]);
for r = 1:nr
    for k = 1:nc
        pos = [L0 + (k - 1) * (w + gx), H - Tt - r * (a + Sb) - (r - 1) * gy, w, a];
        ax = axes(fig, 'Units', 'centimeters', 'Position', pos);  hold(ax, 'on');
        set(ax, 'FontSize', f_tick, 'FontName', 'Times', 'Box', 'on', 'TickDir', 'out', 'TickLength', [0.025 0.025]);
        if overlay
            draw_overlay(ax, D, k, col, f_box);
            if k > 1, ax.YTickLabel = []; end
            continue
        end
        set(ax, 'XLim', D.lim, 'YLim', D.lim, 'XTick', D.tk, 'YTick', D.tk);
        if k > 1, ax.YTickLabel = []; end
        if r < nr, ax.XTickLabel = []; end
        if k == 1, ylabel(ax, sprintf('Month %d', D.months(r)), 'FontSize', f_row, 'FontWeight', 'bold', 'FontName', 'Times'); end
        x = D.x{r};  y = D.y{r, k};
        if ~all(isfinite(y))
            text(ax, 0.5, 0.5, 'diverged', 'Units', 'normalized', 'HorizontalAlignment', 'center', ...
                'FontSize', f_box, 'FontName', 'Times', 'Color', [0.6 0 0]);
            continue
        end
        plot(ax, D.lim, D.lim, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.8);
        in = y >= D.lim(1) & y <= D.lim(2);
        scatter(ax, x(in), y(in), 8, [0.12 0.47 0.71], 'filled');
        if any(~in)
            up = y > D.lim(2);
            scatter(ax, x(up), repmat(D.lim(2), nnz(up), 1), 12, [0.12 0.47 0.71], '^', 'filled', 'Clipping', 'off');
            text(ax, 0.05, 0.95, sprintf('%d off\nscale', nnz(~in)), 'Units', 'normalized', ...
                'VerticalAlignment', 'top', 'FontName', 'Times', 'FontSize', f_tick, 'Color', [0.6 0 0]);
        end
        b = D.b(r, k);  xl = [0, max(x)];
        plot(ax, xl, b * xl, 'r-', 'LineWidth', 1);
        switch mode
            case "line"
                s = sprintf('(%.2f, %s, %s)', b, fmt_r2(D.r2(r, k)), fmt2(D.nmse(r, k)));
            case "stack"
                s = {sprintf('b  %.2f', b), sprintf('R^2  %s', fmt_r2(D.r2(r, k))), sprintf('NMSE  %s', fmt2(D.nmse(r, k)))};
            otherwise
                s = '';
        end
        if ~isempty(s)
            text(ax, 0.5, 1.02, s, 'Units', 'normalized', 'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'bottom', 'FontName', 'Times', 'FontSize', f_box);
        end
    end
end

% labels, title, key, table and legend on an invisible full-figure axes
ov = axes(fig, 'Units', 'centimeters', 'Position', [0 0 W H], 'XLim', [0 W], 'YLim', [0 H]);  hold(ov, 'on');
xc = L0 + (nc * w + (nc - 1) * gx) / 2;
yb = H - Tt - nr * (a + Sb) - (nr - 1) * gy;       % bottom of the lowest row
xk = @(k) L0 + (k - 1) * (w + gx) + w / 2;          % centre of column k
for k = 1:nc
    text(ov, xk(k), H - Tt + 0.1, D.labels(k), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'FontName', 'Times', 'FontSize', f_col, 'FontWeight', 'bold');
end
text(ov, xc, yb - 0.6, '[Tau] (nM) - Numeric NTM', 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'top', 'FontName', 'Times', 'FontSize', f_lab);
text(ov, L0 - 1.5, yb + (nr * (a + Sb) + (nr - 1) * gy) / 2, '[Tau] (nM) - sNTM', 'Rotation', 90, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontName', 'Times', 'FontSize', f_lab);
text(ov, W / 2, H - 0.05, D.title, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
    'FontName', 'Times', 'FontSize', f_title, 'FontWeight', 'bold');
if key
    text(ov, xc, H - 0.8, '(b, R^2, NMSE)', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
        'FontName', 'Times', 'FontSize', f_key);
end
if tbl
    % one block of three lines (b, R^2, NMSE) per month, values under each model
    % header: the model names over their columns (11/19, user)
    y = yb - 1.25 - legH - 0.1;
    for k = 1:nc
        text(ov, xk(k), y, D.labels(k), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
            'FontName', 'Times', 'FontSize', f_box, 'FontWeight', 'bold');
    end
    y = y - lineH - 0.05;
    for c = 1:nmo
        names = ["b", "R^2", "NMSE"];
        for q = 1:3
            yy = y - (q - 1) * lineH;
            if q == 1
                tc = 'k';  if overlay, tc = col(c, :); end
                text(ov, 0.1, yy, sprintf('Month %d', D.months(c)), 'VerticalAlignment', 'top', ...
                    'FontName', 'Times', 'FontSize', f_box, 'FontWeight', 'bold', 'Color', tc);
            end
            text(ov, L0 - 0.1, yy, names(q), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
                'FontName', 'Times', 'FontSize', f_box);
            for k = 1:nc
                switch q
                    case 1, s = sprintf('%.2f', D.b(c, k));
                    case 2, s = fmt_r2(D.r2(c, k));
                    case 3, s = fmt2(D.nmse(c, k));
                end
                text(ov, xk(k), yy, s, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
                    'FontName', 'Times', 'FontSize', f_box);
            end
        end
        y = y - 3 * lineH - 0.2;
        if c < nmo
            plot(ov, [0.1, W - R], (y + 0.12) * [1 1], '-', 'Color', [0.8 0.8 0.8], 'LineWidth', 0.5);
        end
    end
end
% legend, drawn by hand so its position is fixed
ylg = yb - 1.25 - 0.35 - 0.55 * overlay;   % (11/19, user) legend right under the x label, above the table
if overlay
    x0 = W / 2 - 5.2;
    for c = 1:nmo
        xx = x0 + (c - 1) * 2.9;
        scatter(ov, xx + 0.2, ylg + 0.55, 14, col(c, :), 'filled');
        plot(ov, xx + [0.45 1.05], (ylg + 0.55) * [1 1], '-', 'Color', col(c, :), 'LineWidth', 1);
        text(ov, xx + 1.15, ylg + 0.55, sprintf('Month %d', D.months(c)), 'VerticalAlignment', 'middle', ...
            'FontName', 'Times', 'FontSize', f_leg);
    end
    plot(ov, x0 + 8.8 + [0 0.7], (ylg + 0.55) * [1 1], '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.8);
    text(ov, x0 + 9.65, ylg + 0.55, 'y = x', 'VerticalAlignment', 'middle', 'FontName', 'Times', 'FontSize', f_leg);
    text(ov, W / 2, ylg - 0.05, 'Points: regions; solid lines: least-squares fit through the origin, y = bx', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontName', 'Times', 'FontSize', f_leg - 1);
else
    x0 = W / 2 - 6.0;
    plot(ov, x0 + [0 0.8], [ylg ylg], '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.8);
    text(ov, x0 + 0.95, ylg, 'y = x', 'VerticalAlignment', 'middle', 'FontName', 'Times', 'FontSize', f_leg);
    x1 = x0 + 2.75;
    plot(ov, x1 + [0 0.8], [ylg ylg], 'r-', 'LineWidth', 1);
    text(ov, x1 + 0.95, ylg, 'Least-squares fit through the origin, y = bx', 'VerticalAlignment', 'middle', ...
        'FontName', 'Times', 'FontSize', f_leg);
end
set(ov, 'XLim', [0 W], 'YLim', [0 H], 'Visible', 'off');
uistack(ov, 'bottom');
print(fig, out_file + ".png", '-dpng', '-r300');
print(fig, out_file + ".pdf", '-dpdf', '-vector');
fprintf('wrote %s.png/.pdf (%.1f x %.1f cm)\n', out_file, W, H);
close(fig);
end

function draw_overlay(ax, D, k, col, f_box)
% variant D: all months in one log-log panel, one colour per month
lim = D.llim;
set(ax, 'XScale', 'log', 'YScale', 'log', 'XLim', lim, 'YLim', lim, ...
    'XTick', [0.01 0.1 1], 'YTick', [0.01 0.1 1], 'XTickLabel', {'0.01', '', '1'}, 'XTickLabelRotation', 0, ...
    'YTickLabel', {'0.01', '0.1', '1'}, 'XMinorTick', 'off', 'YMinorTick', 'off');
plot(ax, lim, lim, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.8);
nout = 0;
for c = numel(D.months):-1:1                       % month 1 drawn last (on top)
    x = D.x{c};  y = D.y{c, k};
    if ~all(isfinite(y)), continue, end
    ok = x > 0 & y >= lim(1) & y <= lim(2);
    nout = nout + nnz(x > 0 & ~ok);
    scatter(ax, x(ok), y(ok), 7, col(c, :), 'filled', 'MarkerFaceAlpha', 0.85);
    xl = [min(x(x > 0)), max(x)];
    plot(ax, xl, D.b(c, k) * xl, '-', 'Color', col(c, :), 'LineWidth', 1);
end
if nout > 0
    text(ax, 0.05, 0.95, sprintf('%d off\nscale', nout), 'Units', 'normalized', 'VerticalAlignment', 'top', ...
        'FontName', 'Times', 'FontSize', f_box - 1, 'Color', [0.6 0 0]);
end
end

function c = month_colours()
% one hue, light to dark with time (months 1, 3, 5)
c = [0.42 0.68 0.84; 0.13 0.44 0.71; 0.03 0.19 0.42];
end

function t = fmt2(v)
t = sprintf('%#.2g', v);
if t(end) == '.', t = t(1:end-1); end
end

function t = fmt_r2(r2)
if round(r2, 2) >= 1, t = sprintf('%.3f', r2); else, t = sprintf('%.2f', r2); end
end

function T = plot_fig8_gp(out_file, src_fmt, fig_title, num_file, interp_file)
% Fig. 10 of the manuscript (network figure; "Fig. 8" in older names): regional
% tau of each sNTM against the numeric NTM at months 1, 3 and 5.
% One point per region. Red line (11/17, advisor): least-squares fit through
% the origin, y = b x with slope b = sum(x y) / sum(x^2), so that b - 1 is the
% proportional bias of the sNTM; dashed line: y = x. All panels share one range
% (square panels). For each model and month the table under the panels gives b,
% R^2 = squared Pearson correlation (11/17) and the NMSE of Eq. 7 over the
% regions (nmse_eq7.m). Concentrations in nM (N is stored in uM).
% Layout (11/19, user: variant C of plot_fig10_variants.m): months in rows,
% models in columns (LR, SR, 1-NN, GP, MLP from left to right), drawn at the
% printed size (31 pc = 13.1 cm text width of sn-jnl) so that the font sizes
% are the sizes on the page; the legend under the panels; the numbers in a
% table under the legend, with the model names over its columns. The earlier
% layout (models in rows, numbers above each panel) is kept as
% plot_fig8_gp_oneline_backup.m in the session scratchpad only.
%
%   plot_fig8_gp();                 % Fig. 8 parameters (outputs_fig8/)
%   plot_fig8_gp(out, 'outputs_consistent/rec_l2f__%s__CA1__g0.005_l0.02_d20_e80.mat')
%       the Primary Parameter Set (main text Fig. 10)
%   plot_fig8_gp(out, src_fmt, title)       title '' for none
%   plot_fig8_gp(out, src_fmt, title, num_file)
%       numeric run from num_file (relative to network_validation/)
%   plot_fig8_gp(out, src_fmt, title, num_file, interp_file)
%       1-NN run from interp_file (edge/ood sets, fig8_cases.m)
%   T: table of method, month, R2 (squared correlation), NMSE, slope b

here = fileparts(mfilename('fullpath'));
nv = fileparts(here);
addpath(nv);
if nargin < 1, out_file = fullfile(here, 'figures', 'fig8_gp'); end
if ~isfolder(fileparts(out_file)), mkdir(fileparts(out_file)); end
if nargin < 2 || isempty(src_fmt), src_fmt = 'outputs_fig8/fig8__%s__ECl__g5e-05_l0.01_d10_e50.mat'; end
if nargin < 3, fig_title = 'Network wide NTM simulation regional tau accuracy'; end
fig_title = string(fig_title);
file = @(m) fullfile(nv, sprintf(src_fmt, m));
if nargin >= 5 && ~isempty(interp_file)
    file0 = file;  file = @(m) pick(m == "interp", fullfile(nv, interp_file), file0(m));
end
if nargin >= 4 && ~isempty(num_file), num_path = fullfile(nv, num_file); else, num_path = file("numeric"); end
% columns in the paper's model order (11/17), left to right
methods = ["lr", "sr", "interp", "gp", "mlp"];
labels = ["LR", "SR", "1-NN", "GP", "MLP"];
has = arrayfun(@(m) isfile(file(m)), methods);
methods = methods(has);  labels = labels(has);
months = [1 3 5];   % (11/18, Suzanne) months 1, 3 and 5 only (before: 1:5)
nc = numel(methods);  nr = numel(months);

% data and statistics
num = load(num_path, 'N', 'trange');
day = num.trange * 180;
hm = arrayfun(@(m) find(abs(day - 30 * m) == min(abs(day - 30 * m)), 1), months);
T = table('Size', [0 5], 'VariableTypes', ["string", "double", "double", "double", "double"], ...
    'VariableNames', ["method", "month", "R2", "NMSE", "slope"]);
X = cell(1, nr);  Y = cell(nr, nc);  [bb, r2, nm] = deal(nan(nr, nc));  dfail = nan(1, nc);
for r = 1:nr, X{r} = num.N(:, hm(r)) * 1e3; end
% one axis range for every panel (11/13): 0 to the largest finite concentration
% of any model at any of the months shown, (11/18) but at most 3 times the
% largest numeric concentration, so that one runaway sNTM does not squeeze the
% others into a corner; points beyond the range are drawn as triangles on the
% edge of the panel (b, R^2 and NMSE use all points)
nmax = max(num.N(:, hm), [], 'all');  vmax = nmax;
for k = 1:nc
    s = load(file(methods(k)), 'N');
    if any(~isfinite(s.N(:))), dfail(k) = day(find(any(~isfinite(s.N), 1), 1)); end
    for r = 1:nr
        x = X{r};  y = s.N(:, hm(r)) * 1e3;  Y{r, k} = y;
        v = y(isfinite(y) & y <= 3 * nmax * 1e3);
        if ~isempty(v), vmax = max(vmax, max(v) / 1e3); end
        if ~all(isfinite(y))
            % the sNTM simulation diverged (NaN/Inf) before this month
            fprintf('%-6s month %d: diverged at day %.1f\n', methods(k), months(r), dfail(k));
            T(end+1, :) = {methods(k), months(r), NaN, NaN, NaN}; %#ok<AGROW>
            continue
        end
        bb(r, k) = sum(x .* y) / sum(x .^ 2);                % least squares through the origin
        Rc = corrcoef(x, y);  r2(r, k) = Rc(1, 2)^2;
        nm(r, k) = nmse_eq7(y, x);
        fprintf('%-6s month %d: b %.3f, R^2 %.4f, NMSE %.3g\n', methods(k), months(r), bb(r, k), r2(r, k), nm(r, k));
        T(end+1, :) = {methods(k), months(r), r2(r, k), nm(r, k), bb(r, k)}; %#ok<AGROW>
    end
end
lim = [0, 1.05 * vmax * 1e3];
% (11/18) tick spacing from the range: the smallest of 1, 2, 2.5, 4, 5 x 10^k
% giving at most 3 intervals
st = 10 .^ floor(log10(lim(2))) * [0.1 0.2 0.25 0.4 0.5 1 2 2.5 4 5];
tk = 0:st(find(lim(2) ./ st <= 3, 1)):lim(2);

% font sizes in pt as printed
f_tick = 9.5;  f_tab = 10.5;  f_row = 13;  f_col = 13.5;  f_lab = 13.5;  f_leg = 12;  f_title = 13;
% layout in cm
W = 13.1;  L0 = 2.0;  R = 0.15;  gx = 0.3;  gy = 0.25;  Sb = 0.1;  lineH = 0.45;
w = (W - L0 - R - (nc - 1) * gx) / nc;  a = w;          % square panels
Tt = 0.75 * (fig_title ~= "") + 0.65;                   % title, model names
legH = 0.75;                                            % legend
tblH = (nr * 3 + 1) * lineH + (nr - 1) * 0.2 + 0.4;     % table with a header row
Bb = 1.25 + legH + tblH;                                % tick labels, x label, legend, table
H = Bb + nr * (a + Sb) + (nr - 1) * gy + Tt;

fig = figure('Units', 'centimeters', 'Position', [1 1 W H], 'Color', 'w', ...
    'PaperUnits', 'centimeters', 'PaperSize', [W H], 'PaperPosition', [0 0 W H]);
for r = 1:nr
    for k = 1:nc
        pos = [L0 + (k - 1) * (w + gx), H - Tt - r * (a + Sb) - (r - 1) * gy, w, a];
        ax = axes(fig, 'Units', 'centimeters', 'Position', pos);  hold(ax, 'on');
        set(ax, 'FontSize', f_tick, 'FontName', 'Times', 'Box', 'on', 'TickDir', 'out', ...
            'XLim', lim, 'YLim', lim, 'XTick', tk, 'YTick', tk, 'TickLength', [0.025 0.025]);
        if k > 1, ax.YTickLabel = []; end
        if r < nr, ax.XTickLabel = []; end
        if k == 1, ylabel(ax, sprintf('Month %d', months(r)), 'FontSize', f_row, 'FontWeight', 'bold', 'FontName', 'Times'); end
        x = X{r};  y = Y{r, k};
        if ~all(isfinite(y))
            text(ax, 0.5, 0.5, sprintf('diverged\n(day %.0f)', dfail(k)), 'Units', 'normalized', ...
                'HorizontalAlignment', 'center', 'FontSize', f_tick, 'FontName', 'Times', 'Color', [0.6 0 0]);
            continue
        end
        plot(ax, lim, lim, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.8);    % y = x
        in = y >= lim(1) & y <= lim(2);
        scatter(ax, x(in), y(in), 8, [0.12 0.47 0.71], 'filled');
        if any(~in)
            % off scale: triangles on the top (above) or bottom (below) edge
            up = y > lim(2);  dn = y < lim(1);
            scatter(ax, x(up), repmat(lim(2), nnz(up), 1), 12, [0.12 0.47 0.71], '^', 'filled', 'Clipping', 'off');
            scatter(ax, x(dn), repmat(lim(1), nnz(dn), 1), 12, [0.12 0.47 0.71], 'v', 'filled', 'Clipping', 'off');
            text(ax, 0.05, 0.95, sprintf('%d off\nscale', nnz(~in)), 'Units', 'normalized', ...
                'VerticalAlignment', 'top', 'FontName', 'Times', 'FontSize', f_tick, 'Color', [0.6 0 0]);
        end
        xl = [0, max(x)];
        plot(ax, xl, bb(r, k) * xl, 'r-', 'LineWidth', 1);
    end
end

% labels, title, legend and table on an invisible full-figure axes
ov = axes(fig, 'Units', 'centimeters', 'Position', [0 0 W H], 'XLim', [0 W], 'YLim', [0 H]);  hold(ov, 'on');
xc = L0 + (nc * w + (nc - 1) * gx) / 2;
yb = H - Tt - nr * (a + Sb) - (nr - 1) * gy;            % bottom of the lowest row
xk = @(k) L0 + (k - 1) * (w + gx) + w / 2;               % centre of column k
for k = 1:nc
    text(ov, xk(k), H - Tt + 0.1, labels(k), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'FontName', 'Times', 'FontSize', f_col, 'FontWeight', 'bold');
end
text(ov, xc, yb - 0.6, '[Tau] (nM) - Numeric NTM', 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'top', 'FontName', 'Times', 'FontSize', f_lab);
text(ov, L0 - 1.5, yb + (nr * (a + Sb) + (nr - 1) * gy) / 2, '[Tau] (nM) - sNTM', 'Rotation', 90, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontName', 'Times', 'FontSize', f_lab);
if fig_title ~= ""
    text(ov, W / 2, H - 0.05, fig_title, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
        'FontName', 'Times', 'FontSize', f_title, 'FontWeight', 'bold');
end
% legend, right under the x label (11/19, user); drawn by hand so its position is fixed
ylg = yb - 1.25 - 0.35;  x0 = W / 2 - 6.0;
plot(ov, x0 + [0 0.8], [ylg ylg], '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.8);
text(ov, x0 + 0.95, ylg, 'y = x', 'VerticalAlignment', 'middle', 'FontName', 'Times', 'FontSize', f_leg);
x1 = x0 + 2.75;
plot(ov, x1 + [0 0.8], [ylg ylg], 'r-', 'LineWidth', 1);
text(ov, x1 + 0.95, ylg, 'Least-squares fit through the origin, y = bx', 'VerticalAlignment', 'middle', ...
    'FontName', 'Times', 'FontSize', f_leg);
% table: a header row of model names (11/19, user), then one block of three
% lines (b, R^2, NMSE) per month, the values under each model's column
y = yb - 1.25 - legH - 0.1;
for k = 1:nc
    text(ov, xk(k), y, labels(k), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
        'FontName', 'Times', 'FontSize', f_tab, 'FontWeight', 'bold');
end
y = y - lineH - 0.05;
names = ["b", "R^2", "NMSE"];
for r = 1:nr
    for q = 1:3
        yy = y - (q - 1) * lineH;
        if q == 1
            text(ov, 0.1, yy, sprintf('Month %d', months(r)), 'VerticalAlignment', 'top', ...
                'FontName', 'Times', 'FontSize', f_tab, 'FontWeight', 'bold');
        end
        text(ov, L0 - 0.1, yy, names(q), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
            'FontName', 'Times', 'FontSize', f_tab);
        for k = 1:nc
            if isnan(bb(r, k))
                s = '-';                                 % diverged
            else
                switch q
                    case 1, s = sprintf('%.2f', bb(r, k));
                    case 2, s = fmt_r2(r2(r, k));
                    case 3, s = fmt2(nm(r, k));
                end
            end
            text(ov, xk(k), yy, s, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
                'FontName', 'Times', 'FontSize', f_tab);
        end
    end
    y = y - 3 * lineH - 0.2;
    if r < nr
        plot(ov, [0.1, W - R], (y + 0.12) * [1 1], '-', 'Color', [0.8 0.8 0.8], 'LineWidth', 0.5);
    end
end
set(ov, 'XLim', [0 W], 'YLim', [0 H], 'Visible', 'off');
uistack(ov, 'bottom');

print(fig, out_file + ".png", '-dpng', '-r300');
print(fig, out_file + ".pdf", '-dpdf', '-vector');
fprintf('wrote %s.png/.pdf (%.1f x %.1f cm)\n', out_file, W, H);
end

function t = fmt2(v)
% two significant figures, keeping a trailing zero (1.0, not 1; 16, not 16.)
t = sprintf('%#.2g', v);
if t(end) == '.', t = t(1:end-1); end
end

function t = fmt_r2(r2)
% R^2 to two decimals, or three where two would round to 1.00
if round(r2, 2) >= 1, t = sprintf('%.3f', r2); else, t = sprintf('%.2f', r2); end
end

function v = pick(c, a, b)
if c, v = a; else, v = b; end
end

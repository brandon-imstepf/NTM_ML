function make_glass_brain_figure(src_fmt, tag, out_file)
% Network simulations figure of the manuscript (label fig:mdlcomp_net; the
% published images/Figure8.png), with the GP:
%   (a) glass brains (N + M, dorsal view, glass_brain_images.m) for the numeric
%       NTM and each sNTM at the seed (day 0) and days 7, 15, 30, 60;
%   (b) regional N (nM) over days 0-60, one panel per model.
% src_fmt: the runs, relative to network_validation/, %s = method (default: the
% Fig. 8 runs in outputs_fig8/). tag names the panel folder
% figures/glass_brain_<tag>/; panels already there are reused (delete the folder
% to redraw). Writes <out_file>.png and .pdf (default figures/glass_brain_<tag>_figure).
%
%   make_glass_brain_figure();
%   make_glass_brain_figure('outputs_consistent/rec_l2f__%s__ECl__g0.005_l0.05_d20_e80.mat', ...
%       'g0.005_l0.05_d20_e80');

here = fileparts(mfilename('fullpath'));
nv = fileparts(here);
if nargin < 1 || isempty(src_fmt), src_fmt = 'outputs_fig8/fig8__%s__ECl__g5e-05_l0.01_d10_e50.mat'; end
if nargin < 2 || isempty(tag), tag = 'fig8'; end
if nargin < 3 || isempty(out_file), out_file = fullfile(here, 'figures', "glass_brain_" + tag + "_figure"); end
panel_dir = fullfile(here, 'figures', "glass_brain_" + tag);

rows = ["numeric", "mlp", "gp", "sr", "lr"];
row_labels = ["NTM" + newline + "(Ground Truth)", "MLP" + newline + "sNTM", "GP" + newline + "sNTM", ...
    "SR" + newline + "sNTM", "LR" + newline + "sNTM"];
days = [0 7 15 30 60];
col_labels = ["Seed", "Day 7", "Day 15", "Day 30", "Day 60"];
b_order = ["lr", "sr", "gp", "mlp", "numeric"];
b_labels = ["LR sNTM", "SR sNTM", "GP sNTM", "MLP sNTM", "NTM" + newline + "(Ground Truth)"];

% Panels: draw the missing ones
pname = @(m, d) fullfile(panel_dir, sprintf('bf_%s_d%d_sagittal.png', m, d));
for m = rows
    need = days(arrayfun(@(d) ~isfile(pname(m, d)), days));
    if m ~= "numeric", need = need(need ~= 0); end
    if ~isempty(need), glass_brain_images(need, m, panel_dir, src_fmt); end
end

% One crop box for all panels (union of non-white pixels), so they align
files = {};
for m = rows
    for d = days
        if isfile(pname(m, d)), files{end+1} = pname(m, d); end %#ok<AGROW>
    end
end
box = [inf -inf inf -inf];
for k = 1:numel(files)
    im = imread(files{k});
    ink = any(im < 245, 3);
    r = find(any(ink, 2));  c = find(any(ink, 1));
    box = [min(box(1), r(1)) max(box(2), r(end)) min(box(3), c(1)) max(box(4), c(end))];
end
prep = @(f) rot90(subsref(imread(f), struct('type', '()', 'subs', {{box(1):box(2), box(3):box(4), ':'}})), 1);
bg = [0.973 0.976 0.984];
% white render background -> box colour, so the panels blend in
blend = @(im) uint8(double(im) .* ~repmat(all(im >= 250, 3), 1, 1, 3) + ...
    repmat(all(im >= 250, 3), 1, 1, 3) .* reshape(255 * bg, 1, 1, 3));

% Layout (normalized figure units)
fw = 18;  fh = 25;   % cm
fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'centimeters', 'Position', [1 1 fw fh]);
font = 'Times New Roman';
edge = [0.08 0.13 0.18];
a_box = [0.01 0.335 0.98 0.625];  b_box = [0.01 0.01 0.98 0.315];
% background boxes on an axes behind everything (filled annotations would cover the axes)
bax = axes(fig, 'Position', [0 0 1 1], 'XLim', [0 1], 'YLim', [0 1], 'Visible', 'off');
rectangle(bax, 'Position', a_box, 'FaceColor', bg, 'EdgeColor', edge, 'LineWidth', 1);
rectangle(bax, 'Position', b_box, 'FaceColor', bg, 'EdgeColor', edge, 'LineWidth', 1);
annotation(fig, 'textbox', [0 0.962 1 0.035], 'String', 'Hippocampal subnetwork NTM simulations', ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontName', font, 'FontSize', 16, 'FontWeight', 'bold');
annotation(fig, 'textbox', [0.02 0.92 0.04 0.03], 'String', 'a', 'EdgeColor', 'none', ...
    'FontName', font, 'FontSize', 16, 'FontWeight', 'bold');
annotation(fig, 'textbox', [0.02 0.29 0.04 0.03], 'String', 'b', 'EdgeColor', 'none', ...
    'FontName', font, 'FontSize', 16, 'FontWeight', 'bold');

% (a) glass brains
x0 = 0.11;  cw = (0.985 - x0) / numel(days);
y_top = 0.915;  rh = (y_top - 0.345) / numel(rows);
for c = 1:numel(days)
    annotation(fig, 'textbox', [x0 + (c - 1) * cw, y_top, cw, 0.025], 'String', col_labels(c), ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontName', font, 'FontSize', 12, 'FontWeight', 'bold');
end
for r = 1:numel(rows)
    yb = y_top - r * rh;
    ax = axes(fig, 'Position', [0.025 yb 0.07 rh], 'Visible', 'off');
    text(ax, 0.5, 0.5, row_labels(r), 'Rotation', 90, 'HorizontalAlignment', 'center', ...
        'FontName', font, 'FontSize', 12, 'FontWeight', 'bold');
    for c = 1:numel(days)
        f = pname(rows(r), days(c));
        if ~isfile(f), continue; end
        ax = axes(fig, 'Position', [x0 + (c - 1) * cw + 0.004, yb + 0.004, cw - 0.008, rh - 0.008]);
        image(ax, blend(prep(f)));  axis(ax, 'image', 'off');
    end
end

% (b) regional N over days 0-60
num = load(fullfile(nv, sprintf(src_fmt, 'numeric')), 'N', 'trange');
day = num.trange * 180;  keep = day <= 60 + 1e-9;
ymax = 0;
D = struct();
for m = b_order
    s = load(fullfile(nv, sprintf(src_fmt, m)), 'N');
    D.(m) = s.N(:, keep) * 1e3;      % uM -> nM
    ymax = max(ymax, max(D.(m)(:)));
end
tab10 = [31 119 180; 255 127 14; 44 160 44; 214 39 40; 148 103 189; 140 86 75; ...
    227 119 194; 127 127 127; 188 189 34; 23 190 207] / 255;
bw = (0.97 - 0.09) / numel(b_order);
for k = 1:numel(b_order)
    ax = axes(fig, 'Position', [0.09 + (k - 1) * bw + 0.012, 0.06, bw - 0.03, 0.205]);
    hold(ax, 'on');
    set(ax, 'ColorOrder', tab10, 'FontName', font, 'FontSize', 9, 'Box', 'on', 'TickDir', 'out');
    plot(ax, day(keep), D.(b_order(k))', 'LineWidth', 0.9);
    xlim(ax, [-2 62]);  ylim(ax, [-0.03 1.05] * ymax);  xticks(ax, 0:15:60);
    title(ax, b_labels(k), 'FontName', font, 'FontSize', 11, 'FontWeight', 'bold');
    if k == 1, ylabel(ax, '[Tau] (nM)', 'FontName', font, 'FontSize', 12, 'FontWeight', 'bold');
    else, ax.YTickLabel = []; end
end
annotation(fig, 'textbox', [0 0.018 1 0.025], 'String', 'Days', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', 'FontName', font, 'FontSize', 12, 'FontWeight', 'bold');

exportgraphics(fig, out_file + ".png", 'Resolution', 300, 'BackgroundColor', 'w');
exportgraphics(fig, out_file + ".pdf", 'ContentType', 'image', 'Resolution', 300, 'BackgroundColor', 'w');
close(fig);
fprintf('wrote %s.png/.pdf\n', out_file);
end

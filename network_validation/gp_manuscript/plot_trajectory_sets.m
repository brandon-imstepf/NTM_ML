function plot_trajectory_sets(groups, t_end)
% Regional tau trajectories (as in panel b of the network-simulations figure,
% make_glass_brain_figure.m) for every network setting of fig8_cases.m: soluble
% tau N (nM) of each of the 30 regions over days 0..t_end, one panel per model
% (LR, SR, GP, MLP sNTM, numeric NTM), on a common y scale.
% The y scale is set by the numeric run (up to 1.15 x its largest value); an
% sNTM that goes above it is not drawn there and the panel says "off scale (max x)". A
% sNTM simulation that diverged (NaN/Inf) is drawn up to that day and the panel
% says "diverged at day d".
% Figures: figures/trajectory_sets/<group>/<setting>__<seed>.png/.pdf and one
% PDF per group (figures/trajectory_sets/<group>.pdf).
%
%   plot_trajectory_sets();                    % all groups, days 0-60
%   plot_trajectory_sets("interior", 180);     % full 180 days

if nargin < 1 || isempty(groups), groups = ["fig8", "interior", "edge", "ood"]; end
if nargin < 2 || isempty(t_end), t_end = 60; end
here = fileparts(mfilename('fullpath'));
nv = fileparts(here);
fig_root = fullfile(here, 'figures', 'trajectory_sets');
if t_end ~= 60, fig_root = fig_root + "_" + t_end + "d"; end

methods = ["lr", "sr", "gp", "mlp", "numeric"];
labels = ["LR sNTM", "SR sNTM", "GP sNTM", "MLP sNTM", "NTM (Ground Truth)"];
tab10 = [31 119 180; 255 127 14; 44 160 44; 214 39 40; 148 103 189; 140 86 75; ...
    227 119 194; 127 127 127; 188 189 34; 23 190 207] / 255;
font = 'Times New Roman';

cases = fig8_cases(groups);
for k = 1:numel(cases)
    c = cases(k);
    num_path = c.num_file;
    if isempty(num_path), num_path = sprintf(c.src_fmt, 'numeric'); end
    if ~isfile(fullfile(nv, num_path)), fprintf('%s %s %s: no numeric run, skipped\n', c.group, c.setting, c.seed); continue; end
    num = load(fullfile(nv, num_path), 'N', 'trange');
    day = num.trange * 180;  keep = day <= t_end + 1e-9;
    ytop = 1.15 * max(num.N(:, keep), [], 'all') * 1e3;

    v = sscanf(char(c.setting), 'g%f_l%f_d%f_e%f')';
    ttl = sprintf('%s, %s seed: \\gamma = %g, \\lambda = %g, \\delta = %g, \\epsilon = %g', ...
        c.group, c.seed, v(1), v(2), v(3), v(4));
    fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'centimeters', 'Position', [1 1 28 8.5]);
    tl = tiledlayout(fig, 1, numel(methods), 'TileSpacing', 'compact', 'Padding', 'compact');
    for m = 1:numel(methods)
        if methods(m) == "numeric", Y = num.N;
        else, s = load(fullfile(nv, sprintf(c.src_fmt, methods(m))), 'N');  Y = s.N; end
        Y = Y(:, keep) * 1e3;
        ax = nexttile(tl);  hold(ax, 'on');
        set(ax, 'ColorOrder', tab10, 'FontName', font, 'FontSize', 10, 'Box', 'on', 'TickDir', 'out');
        Yp = Y;  Yp(~isfinite(Yp) | Yp > ytop) = NaN;     % hide off-scale and diverged values
        plot(ax, day(keep), Yp', 'LineWidth', 0.9);
        xlim(ax, [-0.03 1.03] * t_end);  ylim(ax, [-0.03 1] * ytop);
        xticks(ax, 0:t_end / 4:t_end);
        title(ax, labels(m), 'FontName', font, 'FontSize', 12, 'FontWeight', 'bold');
        if m > 1, ax.YTickLabel = []; end
        bad = find(any(~isfinite(Y), 1), 1);
        note = "";
        if ~isempty(bad), note = sprintf('diverged at day %.1f', day(bad));
        elseif max(Y, [], 'all') > ytop, note = sprintf('off scale (max %.3g nM)', max(Y, [], 'all')); end
        if note ~= ""
            text(ax, 0.5, 0.95, note, 'Units', 'normalized', 'HorizontalAlignment', 'center', ...
                'FontName', font, 'FontSize', 10, 'Color', [0.6 0 0]);
        end
    end
    ylabel(tl, '[Tau] (nM)', 'FontName', font, 'FontSize', 12, 'FontWeight', 'bold');
    xlabel(tl, 'Days', 'FontName', font, 'FontSize', 12, 'FontWeight', 'bold');
    title(tl, ttl, 'FontName', font, 'FontSize', 13, 'FontWeight', 'bold');
    out_dir = fullfile(fig_root, c.group);
    if ~isfolder(out_dir), mkdir(out_dir); end
    out = fullfile(out_dir, c.setting + "__" + c.seed);
    exportgraphics(fig, out + ".png", 'Resolution', 200, 'BackgroundColor', 'w');
    exportgraphics(fig, out + ".pdf", 'ContentType', 'vector', 'BackgroundColor', 'w');
    close(fig);
    fprintf('wrote %s\n', out);
end

% One PDF per group
for g = string(groups)
    f = dir(fullfile(fig_root, g, '*.pdf'));
    if isempty(f), continue; end
    out_pdf = fullfile(fig_root, g + ".pdf");
    if isfile(out_pdf), delete(out_pdf); end
    for k = 1:numel(f)
        im = imread(fullfile(f(k).folder, replace(f(k).name, ".pdf", ".png")));
        fig = figure('Visible', 'off', 'Units', 'pixels', 'Position', [0 0 size(im, 2) size(im, 1)] / 2);
        ax = axes(fig, 'Position', [0 0 1 1]);
        image(ax, im);  axis(ax, 'image', 'off');
        exportgraphics(fig, out_pdf, 'Append', true, 'Resolution', 200);
        close(fig);
    end
    fprintf('wrote %s (%d pages)\n', out_pdf, numel(f));
end
end

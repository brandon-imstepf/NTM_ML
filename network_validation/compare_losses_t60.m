function C = compare_losses_t60(grid_name, n_points, losses)
% Where the minimum of each loss (nm_loss.m) lies against the numeric truth, from
% the MLP simulations saved by scan_loss_t60.m: one parameter at a time, the
% others at their true values, so the minimum of the scanned curve is what a
% search with unlimited restarts would find. The minimum is refined between scan
% points by a parabola through the lowest point and its two neighbours (search
% coordinates: log10 for gamma1 and lambda1). A loss is better when its minimum
% lies closer to the true value.
%
%   C = compare_losses_t60();                        % dt0.225, 49 points, all losses
%   C = compare_losses_t60('dt0.225', 49, ["nmse", "region"]);
%
% C has one row per loss, parameter and truth: min_value, ratio (min / true),
% at_edge (lowest scan point at an end of the range: not refined).
% Writes results/loss_compare_t60.csv and outputs_nm/t60_scan/loss_compare.png.

if nargin < 1 || isempty(grid_name), grid_name = "dt0.225"; end
if nargin < 2 || isempty(n_points), n_points = 49; end
if nargin < 3 || isempty(losses), losses = ["nmse", "region", "log", "late"]; end
grid_name = string(grid_name);  losses = string(losses);
here = fileparts(mfilename('fullpath'));
params = ["gamma1", "lambda1", "delta", "epsilon"];
param_tex = ["\gamma_1", "\lambda_1", "\delta", "\epsilon"];
T = readtable(fullfile(here, 'manifests', 'recovery_truths.csv'), 'TextType', 'string');
T = T(T.convention == "l2f", :);
n_label = "";
if n_points ~= 25, n_label = "__n" + n_points; end

rows = {};
for p = 1:numel(params)
    S = load(fullfile(here, 'outputs_nm', 't60_scan', "scan_" + grid_name + "__" + params(p) + n_label + ".mat"));
    is_log = p <= 2;
    for s = 1:numel(S.truth_id)
        truth = load(fullfile(here, "outputs_t60_" + grid_name, "rec_l2f__numeric__ECl__" + S.truth_id(s) + ".mat"), 'N', 'trange');
        day = truth.trange * 180;
        true_value = T{T.truth_id == S.truth_id(s), params(p)};
        u = S.value(s, :);
        if is_log, u = log10(u); end
        for l = losses
            f = cellfun(@(N) nm_loss(N, truth.N, l, day), S.N(s, :));
            [u_min, at_edge] = parabola_min(u, f);
            v = u_min;
            if is_log, v = 10^u_min; end
            rows(end+1, :) = {l, params(p), S.truth_id(s), true_value, v, v / true_value, at_edge}; %#ok<AGROW>
        end
    end
end
C = cell2table(rows, 'VariableNames', {'loss', 'param', 'truth_id', 'true_value', 'min_value', 'ratio', 'at_edge'});
writetable(C, fullfile(here, 'results', 'loss_compare_t60.csv'));

%% Summary: |min / true - 1| per loss and parameter
fprintf('== Minimum of each loss / true value (%s, %d scan points), MLP vs numeric truth\n', grid_name, n_points);
for p = params
    fprintf('\n-- %s\n%-24s', p, 'truth_id');
    fprintf('%10s', losses);  fprintf('\n');
    ids = unique(C.truth_id, 'stable');
    for id = ids'
        fprintf('%-24s', id);
        for l = losses
            r = C(C.loss == l & C.param == p & C.truth_id == id, :);
            fprintf('%9.3f%s', r.ratio, ternary(r.at_edge, '*', ' '));
        end
        fprintf('\n');
    end
    fprintf('%-24s', 'median |ratio - 1|');
    for l = losses, fprintf('%10.3f', median(abs(C.ratio(C.loss == l & C.param == p) - 1))); end
    fprintf('\n%-24s', 'max |ratio - 1|');
    for l = losses, fprintf('%10.3f', max(abs(C.ratio(C.loss == l & C.param == p) - 1))); end
    fprintf('\n%-24s', 'within 10%');
    for l = losses, fprintf('%10d', sum(abs(C.ratio(C.loss == l & C.param == p) - 1) <= 0.1)); end
    fprintf('\n');
end
fprintf('* lowest scan point at an end of the range (not refined)\n');

%% Figure: min / true per truth, one panel per parameter, one column of dots per loss
ink = [0.23 0.23 0.22];
col = double([0x2a 0x78 0xd6; 0xeb 0x68 0x34; 0x1b 0xaf 0x7a; 0x8a 0x5c 0xd6]) / 255;
mk = 'osd^';
f = figure('Visible', 'off', 'Color', 'w', 'Position', [0 0 1300 420]);
tl = tiledlayout(f, 1, numel(params), 'TileSpacing', 'compact', 'Padding', 'compact');
for p = 1:numel(params)
    ax = nexttile(tl);  hold(ax, 'on');
    set(ax, 'Box', 'off', 'TickDir', 'out', 'FontSize', 9, 'YScale', 'log', 'XColor', [0.32 0.32 0.31], ...
        'YColor', [0.32 0.32 0.31], 'GridColor', [0.85 0.85 0.84], 'GridAlpha', 1);
    grid(ax, 'on');
    yline(ax, 1, '-', 'Color', ink, 'HandleVisibility', 'off');
    yline(ax, [0.9 1.1], ':', 'Color', ink, 'HandleVisibility', 'off');
    h = gobjects(0);
    for l = 1:numel(losses)
        r = C.ratio(C.loss == losses(l) & C.param == params(p));
        x = l + linspace(-0.25, 0.25, numel(r));
        h(l) = plot(ax, x, r, mk(l), 'MarkerFaceColor', col(l, :), 'MarkerEdgeColor', 'w', 'MarkerSize', 7); %#ok<AGROW>
    end
    set(ax, 'XTick', 1:numel(losses), 'XTickLabel', losses, 'XLim', [0.5 numel(losses) + 0.5], 'YLim', [0.3 3]);
    title(ax, param_tex(p), 'FontWeight', 'normal');
    if p == 1, ylabel(ax, 'loss minimum / true value (dotted: \pm10%)'); end
end
title(tl, sprintf('Where each loss is lowest, MLP vs numeric truth, %s grid (8 truths per column)', grid_name), 'FontSize', 12);
exportgraphics(f, fullfile(here, 'outputs_nm', 't60_scan', 'loss_compare.png'), 'Resolution', 150);
close(f);
end

function [u_min, at_edge] = parabola_min(u, f)
% minimum of f over u, refined by a parabola through the lowest point and its
% neighbours (evenly spaced u)
[~, k] = min(f);
at_edge = k == 1 || k == numel(f);
u_min = u(k);
if at_edge, return, end
d = f(k-1) - 2 * f(k) + f(k+1);
if d > 0, u_min = u(k) - (u(k+1) - u(k)) / 2 * (f(k+1) - f(k-1)) / d; end
end

function out = ternary(c, a, b)
if c, out = a; else, out = b; end
end

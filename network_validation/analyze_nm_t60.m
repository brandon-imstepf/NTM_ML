function R = analyze_nm_t60(out_root, start)
% Tables and figures for the 60-day Nelder-Mead recoveries of run_nm_t60.slurm:
% 8 truths x 4 parameters x 3 truth grids (base, dt0.45, dt0.225) x 2 truths
% (numeric, and the MLP's own run as a control).
%   1. sanity: truth file and grid, exitflag, simulations, recovered value
%      within 1% of the range of a bound
%   2. control (MLP truth): the true value should be recovered, NMSE about 0
%   3. numeric truth: recovered / true for the three grids, and the recovered
%      NMSE against the NMSE of the MLP at the true values (nmse_true): below it,
%      the MLP fits the numeric run better away from the truth (surrogate bias);
%      above it, the search missed a better point. Below it does not show the
%      search found the lowest point; 4 does
%   4. against the loss scan of scan_loss_t60.m (default grid, if run): searches
%      that stopped in a local minimum
%   5. figures (out_root/t60_figures): recovered / true per parameter, recovered
%      NMSE vs nmse_true, and per parameter the scanned loss with the search stops
% R is a table with one row per search; also saved as out_root/t60_figures/nm_t60.csv
% and, tracked in git, results/nm_t60_<start>.csv. sim_time_s is the summed
% simulation time of the search (MATLAB tic/toc per simulation, not SLURM wall time).
% start "mid" (default) or "near": the searches of run_nm_t60.slurm ... near, read
% from out_root/t60_<grid>_near; figures then go to out_root/t60_figures_near.
%
%   R = analyze_nm_t60();
%   R = analyze_nm_t60([], 'near');

here = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(out_root), out_root = fullfile(here, 'outputs_nm'); end
if nargin < 2 || isempty(start), start = "mid"; end
start = string(start);
suffix = "";
if start == "near", suffix = "_near"; end
grids = ["base", "dt0.45", "dt0.225"];
grid_points = [84, 133, 265];
truths = ["numeric", "mlp"];
params = ["gamma1", "lambda1", "delta", "epsilon"];
param_tex = ["\gamma_1", "\lambda_1", "\delta", "\epsilon"];
param_range = [5e-5 5e-2; 0.01 0.1; 10 100; 10 100];
T = readtable(fullfile(here, 'manifests', 'recovery_truths.csv'), 'TextType', 'string');
T = T(T.convention == "l2f", :);
fig_dir = fullfile(out_root, "t60_figures" + suffix);
if ~isfolder(fig_dir), mkdir(fig_dir); end

%% Collect one row per search
rows = {};
for g = 1:numel(grids)
    for tr = truths
        for s = 1:height(T)
            % NMSE of the MLP at the true values against this grid's numeric truth
            num = load(fullfile(here, "outputs_t60_" + grids(g), "rec_l2f__numeric__ECl__" + T.truth_id(s) + ".mat"), 'N');
            mlp = load(fullfile(here, "outputs_t60_" + grids(g), "rec_l2f__mlp__ECl__" + T.truth_id(s) + ".mat"), 'N');
            nmse_true = 0;
            if tr == "numeric", nmse_true = sum((mlp.N(:) - num.N(:)).^2) / sum(num.N(:).^2); end
            for p = 1:numel(params)
                f = fullfile(out_root, "t60_" + grids(g) + suffix, tr + "_truth", ...
                    "nm_l2f__mlp__" + T.truth_id(s) + "__" + params(p) + ".mat");
                S = load(f, 'result', 'trace');
                r = S.result;
                % within 1% of the range of a bound, in the search coordinates
                u = @(x) x;
                if p <= 2, u = @log10; end
                at_bound = min(abs(u(r.recovered_value) - u(param_range(p, :)))) <= ...
                    0.01 * abs(diff(u(param_range(p, :))));
                rows(end+1, :) = {grids(g), tr, T.truth_id(s), params(p), r.true_value, ...
                    r.recovered_value, r.recovered_value / r.true_value, r.nmse, nmse_true, ...
                    r.exitflag, r.n_simulations, at_bound, numel(S.trace(1).N(1, :)) == grid_points(g), ...
                    contains(r.truth_file, "outputs_t60_" + grids(g) + "/rec_l2f__" + tr + "__"), ...
                    sum([S.trace.runtime_s]), mean([S.trace.runtime_s]), ...
                    {[S.trace.value]}, {[S.trace.nmse]}}; %#ok<AGROW>
            end
        end
    end
end
R = cell2table(rows, 'VariableNames', {'grid', 'truth', 'truth_id', 'param', 'true_value', ...
    'recovered', 'ratio', 'nmse', 'nmse_true', 'exitflag', 'n_sims', 'at_bound', 'grid_ok', ...
    'truth_file_ok', 'sim_time_s', 'sim_time_mean_s', 'trace_value', 'trace_nmse'});
writetable(removevars(R, {'trace_value', 'trace_nmse'}), fullfile(fig_dir, 'nm_t60.csv'));
% tracked copy for the paper (outputs_nm is not committed)
if ~isfolder(fullfile(here, 'results')), mkdir(fullfile(here, 'results')); end
writetable(removevars(R, {'trace_value', 'trace_nmse'}), ...
    fullfile(here, 'results', "nm_t60_" + start + ".csv"));

%% 1. Sanity
fprintf('== 1. Sanity (%d searches)\n', height(R));
fprintf('truth file and grid as expected: %d of %d\n', sum(R.grid_ok & R.truth_file_ok), height(R));
R.converged = R.exitflag == 1;
G = groupsummary(R, {'grid', 'truth'}, {'mean', 'max'}, 'n_sims');
G.converged = groupsummary(R, {'grid', 'truth'}, 'sum', 'converged').sum_converged;
G.at_bound = groupsummary(R, {'grid', 'truth'}, 'sum', 'at_bound').sum_at_bound;
disp(G(:, {'grid', 'truth', 'GroupCount', 'converged', 'at_bound', 'mean_n_sims', 'max_n_sims'}))
bad = R(R.exitflag ~= 1 | R.at_bound, :);
if ~isempty(bad)
    fprintf('not converged or at a bound:\n');
    disp(bad(:, {'grid', 'truth', 'truth_id', 'param', 'true_value', 'recovered', 'exitflag', 'n_sims', 'at_bound'}))
end

%% 2. Control: MLP truth
fprintf('\n== 2. Control (MLP truth): recovered / true, worst over the 3 grids\n');
C = R(R.truth == "mlp", :);
print_wide(C, T.truth_id, params, @(x) max(abs(x.ratio - 1)), 'max |ratio - 1|', '%10.2g');
fprintf('largest final NMSE: %.2g\n', max(C.nmse));

%% 3. Numeric truth
fprintf('\n== 3. Numeric truth: recovered / true per grid (base | 0.45 | 0.225)\n');
Nm = R(R.truth == "numeric", :);
fprintf('%-24s', '');  fprintf('%-24s', params);  fprintf('\n');
for s = T.truth_id'
    fprintf('%-24s', s);
    for p = params
        x = Nm(Nm.truth_id == s & Nm.param == p, :);
        fprintf('%-24s', sprintf('%.3f|%.3f|%.3f', x.ratio(x.grid == "base"), x.ratio(x.grid == "dt0.45"), x.ratio(x.grid == "dt0.225")));
    end
    fprintf('\n');
end
fprintf('\nrecovered NMSE / NMSE at the true values (dt0.225; < 1 = surrogate bias, > 1 = search missed)\n');
print_wide(Nm(Nm.grid == "dt0.225", :), T.truth_id, params, @(x) x.nmse / x.nmse_true, '', '%10.3f');
fprintf('spread of recovered / true over the 3 grids, max over settings:\n');
for p = params
    sp = arrayfun(@(s) range(Nm.ratio(Nm.truth_id == s & Nm.param == p)), T.truth_id);
    fprintf('  %-8s %.3g\n', p, max(sp));
end

%% 4. Against the loss scan (scan_loss_t60.m, default grid), if it has been run:
% the scan's lowest NMSE is an upper bound on the lowest NMSE over the range,
% so a search that stopped above it missed a lower region (a local minimum)
scan_dir = fullfile(out_root, 't60_scan');
if isfile(fullfile(scan_dir, "scan_base__" + params(1) + ".mat"))
    fprintf('\n== 4. Default-grid searches vs the loss scan (25 values over the range)\n');
    fprintf('   per search: recovered / true, then the scan''s lowest value / true; LOCAL = the search\n');
    fprintf('   stopped more than max(1e-4, 5%%) above the lowest scanned NMSE, i.e. in a local minimum\n');
    fprintf('%-24s', '');  fprintf('%-26s', params);  fprintf('\n');
    for tr = truths
        fprintf('-- %s truth\n', tr);
        for s = T.truth_id'
            fprintf('%-24s', s);
            for p = params
                Sc = load(fullfile(scan_dir, "scan_base__" + p + ".mat"));
                k = find(Sc.truth_id == s);
                e = Sc.("nmse_" + tr)(k, :);
                [e_min, j] = min(e);
                x = R(R.grid == "base" & R.truth == tr & R.truth_id == s & R.param == p, :);
                is_local = x.nmse - e_min > max(1e-4, 0.05 * e_min);
                flag = '';
                if is_local, flag = 'LOCAL'; end
                fprintf('%-26s', sprintf('%.3f  scan %.3f %s', x.ratio, Sc.value(k, j) / x.true_value, flag));
            end
            fprintf('\n');
        end
    end
end

%% 5. Figures
plot_recovery(R, T, params, param_tex, grids, fig_dir);
plot_bias(Nm, params, param_tex, grids, fig_dir);
for p = 1:numel(params)
    plot_loss(R, T, params(p), param_tex(p), param_range(p, :), p <= 2, grids, fig_dir);
end
fprintf('\nfigures in %s\n', fig_dir);
end

%% ---------- helpers ----------
function print_wide(X, ids, params, fun, label, fmt)
% one row per truth, one column per parameter
fprintf('%-24s', label);  fprintf('%10s', params);  fprintf('\n');
for s = ids'
    fprintf('%-24s', s);
    for p = params
        fprintf(fmt, fun(X(X.truth_id == s & X.param == p, :)));
    end
    fprintf('\n');
end
end

function [col, mk, gray] = style()
% reference categorical palette, first three slots (validated all-pairs), one per grid;
% marker shape is the second encoding; gray for the MLP control
col = double([0x2a 0x78 0xd6; 0xeb 0x68 0x34; 0x1b 0xaf 0x7a]) / 255;
mk = ["o", "s", "^"];
gray = [0.55 0.55 0.53];
end

function ax = clean_axes(ax)
set(ax, 'Box', 'off', 'TickDir', 'out', 'FontSize', 9, 'XColor', [0.32 0.32 0.31], ...
    'YColor', [0.32 0.32 0.31], 'GridColor', [0.85 0.85 0.84], 'GridAlpha', 1);
grid(ax, 'on');
end

function label = short_label(id)
% g0.0005_l0.02_d80_e20 -> g5e-4 l.02 d80 e20
t = split(id, "_");
label = sprintf('\\gamma%s \\lambda%s \\delta%s \\epsilon%s', ...
    regexprep(sprintf('%.0e', str2double(extractAfter(t(1), 1))), 'e-0', 'e-'), ...
    extractAfter(t(2), 1), extractAfter(t(3), 1), extractAfter(t(4), 1));
end

function plot_recovery(R, T, params, param_tex, grids, fig_dir)
[col, mk, gray] = style();
f = figure('Visible', 'off', 'Color', 'w', 'Position', [0 0 1400 520]);
tl = tiledlayout(f, 1, numel(params), 'TileSpacing', 'compact', 'Padding', 'compact');
labels = arrayfun(@short_label, T.truth_id, 'UniformOutput', false);
off = [-0.2 0 0.2];
for p = 1:numel(params)
    ax = clean_axes(nexttile(tl));  hold(ax, 'on');
    yline(ax, 1, '-', 'Color', [0.2 0.2 0.2], 'LineWidth', 1, 'HandleVisibility', 'off');
    h = gobjects(0);
    for g = 1:numel(grids)
        for tr = ["mlp", "numeric"]
            x = R(R.grid == grids(g) & R.truth == tr & R.param == params(p), :);
            [~, k] = ismember(x.truth_id, T.truth_id);
            if tr == "mlp"
                hm = plot(ax, k + off(g), x.ratio, mk(g), 'Color', gray, 'MarkerSize', 6, 'LineWidth', 1);
            else
                h(g) = plot(ax, k + off(g), x.ratio, mk(g), 'MarkerFaceColor', col(g, :), ...
                    'MarkerEdgeColor', 'w', 'MarkerSize', 8, 'LineWidth', 1); %#ok<AGROW>
            end
        end
    end
    set(ax, 'YScale', 'log', 'YLim', [0.45 2.4], 'YTick', [0.5 0.75 1 1.5 2], ...
        'XLim', [0.5 height(T) + 0.5], 'XTick', 1:height(T), 'XTickLabel', labels, ...
        'XTickLabelRotation', 60, 'TickLabelInterpreter', 'tex');
    title(ax, param_tex(p), 'FontWeight', 'normal', 'FontSize', 12);
    if p == 1, ylabel(ax, 'recovered / true'); end
end
lg = legend([h, hm], ["numeric truth, default grid", "numeric truth, 0.45 d", "numeric truth, 0.225 d", ...
    "MLP truth (control), all grids"], 'Orientation', 'horizontal', 'Box', 'off');
lg.Layout.Tile = 'south';
title(tl, 'Recovered / true value, 60 days, one parameter searched at a time', 'FontSize', 13);
exportgraphics(f, fullfile(fig_dir, 'recovery_ratio.png'), 'Resolution', 150);
close(f);
end

function plot_bias(Nm, params, param_tex, grids, fig_dir)
[col, mk] = style();
f = figure('Visible', 'off', 'Color', 'w', 'Position', [0 0 1400 400]);
tl = tiledlayout(f, 1, numel(params), 'TileSpacing', 'compact', 'Padding', 'compact');
lims = [min([Nm.nmse; Nm.nmse_true]) / 2, max([Nm.nmse; Nm.nmse_true]) * 2];
for p = 1:numel(params)
    ax = clean_axes(nexttile(tl));  hold(ax, 'on');
    plot(ax, lims, lims, '-', 'Color', [0.2 0.2 0.2], 'HandleVisibility', 'off');
    h = gobjects(0);
    for g = 1:numel(grids)
        x = Nm(Nm.grid == grids(g) & Nm.param == params(p), :);
        h(g) = plot(ax, x.nmse_true, x.nmse, mk(g), 'MarkerFaceColor', col(g, :), ...
            'MarkerEdgeColor', 'w', 'MarkerSize', 8); %#ok<AGROW>
    end
    set(ax, 'XScale', 'log', 'YScale', 'log', 'XLim', lims, 'YLim', lims);
    axis(ax, 'square');
    title(ax, param_tex(p), 'FontWeight', 'normal', 'FontSize', 12);
    xlabel(ax, 'NMSE of the MLP at the true values');
    if p == 1, ylabel(ax, 'NMSE at the recovered value'); end
end
lg = legend(h, ["default grid", "0.45 d", "0.225 d"], 'Orientation', 'horizontal', 'Box', 'off');
lg.Layout.Tile = 'south';
title(tl, 'Numeric truth: below the diagonal the MLP fits better away from the truth; above it the search missed', 'FontSize', 12);
exportgraphics(f, fullfile(fig_dir, 'nmse_recovered_vs_true.png'), 'Resolution', 150);
close(f);
end

function plot_loss(R, T, param, ptex, prange, is_log, grids, fig_dir)
% The loss over the whole range (scan_loss_t60.m, default grid, if it has been
% run) and where each search stopped. NMSE is floored at 1e-6 for the log axis.
[col, mk, gray] = style();
ink = [0.23 0.23 0.22];
floor_nmse = 1e-6;
scan_file = fullfile(fileparts(fig_dir), 't60_scan', "scan_base__" + param + ".mat");
has_scan = isfile(scan_file);
if has_scan, Sc = load(scan_file); end
f = figure('Visible', 'off', 'Color', 'w', 'Position', [0 0 1400 640]);
tl = tiledlayout(f, 2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
for s = 1:height(T)
    ax = clean_axes(nexttile(tl));  hold(ax, 'on');
    hs = gobjects(0);
    if has_scan
        k = find(Sc.truth_id == T.truth_id(s));
        hs(1) = plot(ax, Sc.value(k, :), max(Sc.nmse_numeric(k, :), floor_nmse), '-', 'Color', ink, 'LineWidth', 2);
        hs(2) = plot(ax, Sc.value(k, :), max(Sc.nmse_mlp(k, :), floor_nmse), '--', 'Color', gray, 'LineWidth', 1.5);
    end
    h = gobjects(0);
    for g = 1:numel(grids)
        c = R(R.grid == grids(g) & R.truth == "mlp" & R.param == param & R.truth_id == T.truth_id(s), :);
        hc = plot(ax, c.recovered, max(c.nmse, floor_nmse), mk(g), 'Color', gray, 'MarkerSize', 7, 'LineWidth', 1);
        x = R(R.grid == grids(g) & R.truth == "numeric" & R.param == param & R.truth_id == T.truth_id(s), :);
        h(g) = plot(ax, x.recovered, x.nmse, mk(g), 'MarkerFaceColor', col(g, :), ...
            'MarkerEdgeColor', 'w', 'MarkerSize', 8); %#ok<AGROW>
    end
    ht = xline(ax, x.true_value, ':', 'Color', ink, 'LineWidth', 1.5);
    set(ax, 'YScale', 'log', 'XLim', prange, 'YLim', [floor_nmse / 2, 3]);
    if is_log, set(ax, 'XScale', 'log'); end
    title(ax, short_label(T.truth_id(s)), 'FontWeight', 'normal', 'FontSize', 9, 'Interpreter', 'tex');
    if mod(s, 4) == 1, ylabel(ax, 'NMSE'); end
    if s > 4, xlabel(ax, ptex); end
end
names = ["stop, numeric truth, default grid", "stop, numeric truth, 0.45 d", "stop, numeric truth, 0.225 d", "stop, MLP truth (control)", "true value"];
handles = [h, hc, ht];
if has_scan
    handles = [hs, handles];
    names = ["loss vs numeric truth (scan)", "loss vs MLP truth (scan)", names];
end
lg = legend(handles, names, 'Box', 'off', 'NumColumns', 4);
lg.Layout.Tile = 'south';
title(tl, sprintf('Loss over the range of %s (others at their true values, 60 days) and where each search stopped', ptex), 'FontSize', 13);
exportgraphics(f, fullfile(fig_dir, "loss_" + param + ".png"), 'Resolution', 150);
close(f);
end

function R = analyze_nm_consistent()
% Tables and figures for the network recoveries of nm_consistent.m (MLP sNTM
% fitted to the consistent numeric runs, 180 days), for the parameter-recovery
% report sent with the reviewer response.
%
% Tables (results/, tracked in git):
%   nm_consistent_single.csv   one row per one-parameter search: seed, truth
%                              (numeric or mlp control), recovered value and
%                              recovered / true, NMSE at the recovered value and
%                              of the MLP at the true values, simulations, exit flag
%   nm_consistent_scan.csv     per truth and parameter: the lowest point of the
%                              25-point loss scan (refined by a parabola through it
%                              and its neighbours in the search coordinate), for
%                              the numeric and the MLP truth
%   nm_consistent_joint.csv    one row per four-parameter search
% Figures (outputs_nm/consistent/figures/, not in git):
%   scan_<param>.png           loss along one parameter, one panel per truth
%   ratio_single.png           recovered / true, all one-parameter searches
%   search_example.png         the four one-parameter searches of one truth
%   joint_trace.png            the four-parameter searches, value / true per simulation
%   joint_regions.png          regional tau for one four-parameter search
%
%   R = analyze_nm_consistent();

here = fileparts(mfilename('fullpath'));
root = fullfile(here, 'outputs_nm', 'consistent');
fig_dir = fullfile(root, 'figures');
if ~isfolder(fig_dir), mkdir(fig_dir); end
T = readtable(fullfile(here, 'manifests', 'recovery_truths.csv'), 'TextType', 'string');
T = T(T.convention == "l2f", :);
params = ["gamma1", "lambda1", "delta", "epsilon"];
ptex = ["\gamma_1", "\lambda_1", "\delta", "\epsilon"];
is_log = [true, true, false, false];
prange = [5e-5 5e-2; 0.01 0.1; 10 100; 10 100];
seeds = ["ECl", "CA1", "RSPd"];
truth_label = arrayfun(@(k) sprintf('\\gamma_1=%g, \\lambda_1=%g, \\delta=%g, \\epsilon=%g', ...
    T.gamma1(k), T.lambda1(k), T.delta(k), T.epsilon(k)), 1:height(T), 'UniformOutput', false);
col = struct('ECl', [42 120 214]/255, 'CA1', [235 104 52]/255, 'RSPd', [27 175 122]/255);
grey = [0.45 0.45 0.45];

%% NMSE of the MLP at the true values against each numeric truth
nmse_true = NaN(height(T), numel(seeds));
for s = 1:height(T)
    for d = 1:numel(seeds)
        num = load(truth_file("numeric", seeds(d), T.truth_id(s)), 'N');
        mlp = load(truth_file("mlp", seeds(d), T.truth_id(s)), 'N');
        nmse_true(s, d) = nm_loss(mlp.N, num.N, "nmse");
    end
end

%% One-parameter searches
rows = cell(0, 12);
% truth: numeric (start mid-range), mlp (control), twostage (numeric truth, start at the scan's lowest point)
dirs = struct('numeric', "numeric_truth", 'mlp', "mlp_truth", 'twostage', "numeric_truth_twostage");
for tr = ["numeric", "mlp", "twostage"]
    for d = 1:numel(seeds)
        if tr ~= "numeric" && seeds(d) ~= "ECl", continue; end
        for s = 1:height(T)
            for p = 1:4
                f = fullfile(root, seeds(d), dirs.(tr), "nm_l2f__mlp__" + T.truth_id(s) + "__" + params(p) + ".mat");
                if ~isfile(f) || isempty(whos('-file', f, 'result')), continue; end
                S = load(f, 'result');  r = S.result;
                nt = 0;  if tr ~= "mlp", nt = nmse_true(s, d); end
                rows(end+1, :) = {seeds(d), tr, T.truth_id(s), params(p), r.true_value, r.recovered_value, ...
                    r.recovered_value / r.true_value, rel_err(r.recovered_value, r.true_value, p), r.nmse, nt, ...
                    r.n_simulations, r.exitflag}; %#ok<AGROW>
            end
        end
    end
end
R.single = cell2table(rows, 'VariableNames', {'seed', 'truth', 'truth_id', 'param', 'true_value', 'recovered', ...
    'ratio', 'error_frac_range', 'nmse', 'nmse_true', 'n_sims', 'exitflag'});
writetable(R.single, fullfile(here, 'results', 'nm_consistent_single.csv'));

%% Loss scans
rows = cell(0, 10);
scan = struct();
for s = 1:height(T)
    for p = 1:4
        f = fullfile(root, 'scan', "scan__" + T.truth_id(s) + "__" + params(p) + ".mat");
        if ~isfile(f), continue; end
        S = load(f, 'value', 'nmse_numeric', 'nmse_mlp', 'true_value');
        scan(s, p).value = S.value;  scan(s, p).num = S.nmse_numeric;  scan(s, p).mlp = S.nmse_mlp;
        [xn, ln] = scan_min(S.value, S.nmse_numeric, is_log(p));
        [xm, lm] = scan_min(S.value, S.nmse_mlp, is_log(p));
        rows(end+1, :) = {T.truth_id(s), params(p), S.true_value, xn, xn / S.true_value, ...
            rel_err(xn, S.true_value, p), ln, nmse_true(s, 1), xm / S.true_value, lm}; %#ok<AGROW>
    end
end
R.scan = cell2table(rows, 'VariableNames', {'truth_id', 'param', 'true_value', 'scan_min_numeric', ...
    'ratio_numeric', 'error_frac_range_numeric', 'nmse_min_numeric', 'nmse_true', 'ratio_mlp', 'nmse_min_mlp'});
writetable(R.scan, fullfile(here, 'results', 'nm_consistent_scan.csv'));

%% Four-parameter searches
rows = {};
joint = struct('trace', {}, 'result', {});
for s = 1:height(T)
    f = fullfile(root, 'ECl', 'numeric_truth', "nm_l2f__mlp__" + T.truth_id(s) + "__gamma1+lambda1+delta+epsilon.mat");
    if ~isfile(f) || isempty(whos('-file', f, 'result')), continue; end
    S = load(f, 'result', 'trace', 'best');
    joint(s).trace = S.trace;  joint(s).result = S.result;  joint(s).best = S.best;
    r = S.result;
    rows(end+1, :) = {T.truth_id(s), r.recovered_value(1) / T.gamma1(s), r.recovered_value(2) / T.lambda1(s), ...
        r.recovered_value(3) / T.delta(s), r.recovered_value(4) / T.epsilon(s), r.nmse, nmse_true(s, 1), ...
        r.n_simulations, r.exitflag, sum([S.trace.runtime_s])}; %#ok<AGROW>
end
if ~isempty(rows)
    R.joint = cell2table(rows, 'VariableNames', {'truth_id', 'ratio_gamma1', 'ratio_lambda1', 'ratio_delta', ...
        'ratio_epsilon', 'nmse', 'nmse_true', 'n_sims', 'exitflag', 'sim_time_s'});
    writetable(R.joint, fullfile(here, 'results', 'nm_consistent_joint.csv'));
end

%% Figure: loss scans, one figure per parameter
for p = 1:4
    fig = figure('Visible', 'off', 'Position', [0 0 1500 640]);
    tl = tiledlayout(fig, 2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
    for s = 1:height(T)
        ax = nexttile(tl);  hold(ax, 'on');  box(ax, 'on');  grid(ax, 'on');
        if s > size(scan, 1) || isempty(scan(s, p).value), continue; end
        v = scan(s, p).value;
        h1 = plot(ax, v, scan(s, p).num, '-', 'Color', 'k', 'LineWidth', 1.5);
        h2 = plot(ax, v, scan(s, p).mlp, '--', 'Color', grey, 'LineWidth', 1);
        tv = T{s, params(p)};
        h3 = xline(ax, tv, ':', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.2);
        hs = gobjects(0);
        k = false;
        if height(R.single) > 0
            k = R.single.seed == "ECl" & R.single.truth == "numeric" & R.single.truth_id == T.truth_id(s) & R.single.param == params(p);
        end
        if any(k)
            hs = plot(ax, R.single.recovered(k), R.single.nmse(k), 'o', 'MarkerSize', 8, ...
                'MarkerFaceColor', [0.85 0.1 0.1], 'MarkerEdgeColor', 'w', 'LineWidth', 1);
        end
        set(ax, 'YScale', 'log');
        if is_log(p), set(ax, 'XScale', 'log'); end
        xlim(ax, prange(p, :));  ylim(ax, [1e-5 3]);   % the control curve is ~0 at the true value
        title(ax, truth_label{s}, 'FontWeight', 'normal', 'FontSize', 8);
        xlabel(ax, ptex(p));
        if mod(s, 4) == 1, ylabel(ax, 'NMSE'); end
    end
    leg = {'against the numeric NTM', 'against the MLP sNTM at the true values', 'true value'};
    hh = [h1, h2, h3];
    if ~isempty(hs), hh(end+1) = hs; leg{end+1} = 'Nelder-Mead estimate (ECl)'; end
    lg = legend(hh, leg, 'Orientation', 'horizontal');  lg.Layout.Tile = 'south';
    exportgraphics(fig, fullfile(fig_dir, "scan_" + params(p) + ".png"), 'Resolution', 200);
    close(fig);
end

%% Figure: recovered / true for every one-parameter search
if height(R.single) > 0
fig = figure('Visible', 'off', 'Position', [0 0 1500 420]);
tl = tiledlayout(fig, 1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
mk = struct('ECl', 'o', 'CA1', 's', 'RSPd', 'd');
for p = 1:4
    ax = nexttile(tl);  hold(ax, 'on');  box(ax, 'on');  grid(ax, 'on');
    yline(ax, 1, '-', 'Color', [0.6 0.6 0.6]);
    h = gobjects(0);  lab = {};
    for d = 1:numel(seeds)
        k = R.single.seed == seeds(d) & R.single.truth == "numeric" & R.single.param == params(p);
        x = arrayfun(@(id) find(T.truth_id == id), R.single.truth_id(k)) + (d - 2) * 0.18;
        y = R.single.ratio(k);
        if isempty(x), x = NaN; y = NaN; end   % keeps the legend entry
        h(end+1) = plot(ax, x, y, mk.(seeds(d)), 'MarkerSize', 7, 'MarkerFaceColor', col.(seeds(d)), ...
            'MarkerEdgeColor', 'w', 'LineStyle', 'none'); %#ok<AGROW>
        lab{end+1} = "Nelder-Mead, numeric truth, " + seeds(d); %#ok<AGROW>
    end
    k = R.scan.param == params(p);
    x = arrayfun(@(id) find(T.truth_id == id), R.scan.truth_id(k));
    y = R.scan.ratio_numeric(k);
    if isempty(x), x = NaN; y = NaN; end
    h(end+1) = plot(ax, x, y, '_', 'Color', 'k', 'MarkerSize', 14, 'LineWidth', 1.5);
    lab{end+1} = 'lowest point of the loss scan (ECl)';
    k = R.single.truth == "twostage" & R.single.param == params(p);
    if any(k)
        x = arrayfun(@(id) find(T.truth_id == id), R.single.truth_id(k)) + 0.32;
        h(end+1) = plot(ax, x, R.single.ratio(k), '^', 'MarkerSize', 7, 'MarkerFaceColor', 'k', ...
            'MarkerEdgeColor', 'w', 'LineStyle', 'none');
        lab{end+1} = 'two-stage (scan, then Nelder-Mead), ECl';
    end
    k = R.single.seed == "ECl" & R.single.truth == "mlp" & R.single.param == params(p);
    if any(k)
        x = arrayfun(@(id) find(T.truth_id == id), R.single.truth_id(k));
        h(end+1) = plot(ax, x, R.single.ratio(k), 'x', 'Color', grey, 'MarkerSize', 7, 'LineWidth', 1.2);
        lab{end+1} = 'Nelder-Mead, MLP truth (control), ECl';
    end
    set(ax, 'YScale', 'log', 'XTick', 1:height(T), 'XTickLabel', compose('%d', 1:height(T)));
    xlim(ax, [0.4, height(T) + 0.6]);  ylim(ax, [0.25 4]);
    set(ax, 'YTick', [0.25 0.5 1 2 4]);
    title(ax, ptex(p), 'FontWeight', 'normal');
    xlabel(ax, 'parameter setting');
    if p == 1, ylabel(ax, 'recovered / true'); end
end
lg = legend(h, lab, 'Orientation', 'horizontal', 'NumColumns', 3);  lg.Layout.Tile = 'south';
exportgraphics(fig, fullfile(fig_dir, 'ratio_single.png'), 'Resolution', 200);
close(fig);
end

%% Figure: the four one-parameter searches of one truth (ECl, numeric truth)
ex = 5;   % g0.005_l0.02_d80_e20
fig = figure('Visible', 'off', 'Position', [0 0 1500 560]);
tl = tiledlayout(fig, 2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
for p = 1:4
    f = fullfile(root, 'ECl', 'numeric_truth', "nm_l2f__mlp__" + T.truth_id(ex) + "__" + params(p) + ".mat");
    ax1 = nexttile(tl, p);  ax2 = nexttile(tl, p + 4);
    if ~isfile(f), continue; end
    S = load(f, 'trace', 'result');
    v = [S.trace.value];  l = [S.trace.nmse];
    hold(ax1, 'on');  box(ax1, 'on');  grid(ax1, 'on');
    plot(ax1, 1:numel(v), v, '-o', 'Color', col.ECl, 'MarkerFaceColor', col.ECl, 'MarkerSize', 4);
    yline(ax1, T{ex, params(p)}, ':', 'true', 'LineWidth', 1.2);
    if is_log(p), set(ax1, 'YScale', 'log'); end
    title(ax1, ptex(p), 'FontWeight', 'normal');
    if p == 1, ylabel(ax1, 'value simulated'); end
    hold(ax2, 'on');  box(ax2, 'on');  grid(ax2, 'on');
    plot(ax2, 1:numel(l), l, 'o', 'Color', col.ECl, 'MarkerFaceColor', col.ECl, 'MarkerSize', 4);
    plot(ax2, 1:numel(l), cummin(l), '-', 'Color', 'k');
    yline(ax2, nmse_true(ex, 1), '--', 'MLP at true values', 'Color', grey);
    set(ax2, 'YScale', 'log');
    xlabel(ax2, 'simulation');
    if p == 1, ylabel(ax2, 'NMSE'); end
end
title(tl, "One-parameter Nelder-Mead searches, " + truth_label{ex} + ", ECl seed", 'FontSize', 10);
exportgraphics(fig, fullfile(fig_dir, 'search_example.png'), 'Resolution', 200);
close(fig);

%% Figure: four-parameter searches
if ~isempty(joint)
    fig = figure('Visible', 'off', 'Position', [0 0 1500 640]);
    tl = tiledlayout(fig, 2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
    for s = 1:height(T)
        ax = nexttile(tl);  hold(ax, 'on');  box(ax, 'on');  grid(ax, 'on');
        if s > numel(joint) || isempty(joint(s).trace), continue; end
        V = vertcat(joint(s).trace.value);   % simulations x 4
        tv = [T.gamma1(s), T.lambda1(s), T.delta(s), T.epsilon(s)];
        nm = [joint(s).trace.nmse];
        bestsofar = zeros(size(V));
        ib = 1;
        for k = 1:size(V, 1)
            if nm(k) <= nm(ib), ib = k; end
            bestsofar(k, :) = V(ib, :);
        end
        pc = lines(4);
        for p = 1:4
            plot(ax, 1:size(V, 1), bestsofar(:, p) / tv(p), '-', 'Color', pc(p, :), 'LineWidth', 1.4);
        end
        yline(ax, 1, ':', 'Color', 'k');
        set(ax, 'YScale', 'log');  ylim(ax, [0.1 10]);
        title(ax, sprintf('%s\nfinal NMSE %.2g (MLP at truth %.2g)', truth_label{s}, joint(s).result.nmse, nmse_true(s, 1)), ...
            'FontWeight', 'normal', 'FontSize', 8);
        xlabel(ax, 'simulation');
        if mod(s, 4) == 1, ylabel(ax, 'best value so far / true'); end
    end
    lg = legend(ptex, 'Orientation', 'horizontal');  lg.Layout.Tile = 'south';
    exportgraphics(fig, fullfile(fig_dir, 'joint_trace.png'), 'Resolution', 200);
    close(fig);

    % regional tau for one four-parameter search
    ex = 5;
    if numel(joint) >= ex && ~isempty(joint(ex).trace)
        num = load(truth_file("numeric", "ECl", T.truth_id(ex)), 'N', 'trange', 'region_names');
        mlp = load(truth_file("mlp", "ECl", T.truth_id(ex)), 'N');
        Nb = joint(ex).best.N;
        day = num.trange * 180;
        fig = figure('Visible', 'off', 'Position', [0 0 1500 1100]);
        tl = tiledlayout(fig, 5, 6, 'TileSpacing', 'compact', 'Padding', 'compact');
        for i = 1:size(num.N, 1)
            ax = nexttile(tl);  hold(ax, 'on');  box(ax, 'on');
            h1 = plot(ax, day, num.N(i, :), '-', 'Color', 'k', 'LineWidth', 1.6);
            h2 = plot(ax, day, mlp.N(i, :), '--', 'Color', col.ECl, 'LineWidth', 1.2);
            h3 = plot(ax, day, Nb(i, :), ':', 'Color', [0.85 0.1 0.1], 'LineWidth', 1.6);
            title(ax, strrep(string(num.region_names{i}), '_', ' '), 'FontWeight', 'normal', 'FontSize', 7);
            set(ax, 'FontSize', 7);  xlim(ax, [0 180]);
        end
        lg = legend([h1 h2 h3], {'numeric NTM (truth)', 'MLP sNTM at the true values', ...
            'MLP sNTM at the four-parameter estimate'}, 'Orientation', 'horizontal');
        lg.Layout.Tile = 'south';
        xlabel(tl, 'day');  ylabel(tl, 'soluble tau N');
        exportgraphics(fig, fullfile(fig_dir, 'joint_regions.png'), 'Resolution', 200);
        close(fig);
    end
end
if height(R.single) > 0, disp(R.single(R.single.truth == "numeric" & R.single.seed == "ECl", :)), end
disp(R.scan)
if isfield(R, 'joint'), disp(R.joint), end

    function f = truth_file(method, seed, id)
        f = fullfile(here, 'outputs_consistent', "rec_l2f__" + method + "__" + seed + "__" + id + ".mat");
    end
    function e = rel_err(x, t, p)
        if is_log(p), e = abs(log10(x) - log10(t)) / diff(log10(prange(p, :)));
        else, e = abs(x - t) / diff(prange(p, :)); end
    end
end

function [x_min, l_min] = scan_min(v, l, is_log)
% lowest point of a scan, refined by a parabola through it and its neighbours
% (in log10 for gamma1, lambda1); the parabola's vertex is kept only if it lies
% between the neighbours and is below the lowest point
u = v;  if is_log, u = log10(v); end
[l_min, k] = min(l);  x_min = v(k);
if k == 1 || k == numel(l), return; end
c = polyfit(u(k-1:k+1), l(k-1:k+1), 2);
if c(1) <= 0, return; end
uv = -c(2) / (2 * c(1));
lv = polyval(c, uv);
if uv > u(k-1) && uv < u(k+1) && lv < l_min && lv > 0
    l_min = lv;  x_min = uv;  if is_log, x_min = 10^uv; end
end
end

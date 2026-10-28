function rn = analyze_network_recovery(conv, manifest_file, rec_dir)
% Single-parameter recovery on the network, from the simulations of a
% make_recovery_manifest.py manifest.
%
%   analyze_network_recovery('l2eq')                                       % recovery_single_l2eq.csv
%   analyze_network_recovery('l2eq', 'manifests/recovery_pilot_l2eq.csv')
%
% For each truth (manifests/recovery_truths.csv), model and parameter, the
% candidates are the manifest rows whose other three parameters equal the
% truth's. Loss of a candidate = NMSE of its N against the numeric truth's N over
% all regions and time points, sum (s - n)^2 / sum n^2 (Eq. 7 of the paper).
% Control loss = the same against the surrogate's own run at the true values.
% Recovered value = grid minimum, refined by a parabola through it and its two
% neighbours in search coordinates (log10 for gamma1 and lambda1). A curve that
% does not change (max - min <= 1e-12 max) is flat: no recovered value.
% Missing simulations are NaN and are counted, so this runs on a partial set.
%
% conv          'l2eq' or 'l2f'
% manifest_file default manifests/recovery_single_<conv>.csv
% rec_dir       folder of the candidate (and l2f truth) runs, default outputs_recovery/
%
% Writes to <rec_dir>/analysis_<manifest name>/:
%   recovery_network.mat     rn (curves and recovered values)
%   recovery_summary.csv     one row per truth, model and parameter
%   single_<truth>.png       loss vs each parameter (rows) for each model (columns)
%   true_vs_recovered.png    recovered vs true value, all truths

here = fileparts(mfilename('fullpath'));
if nargin < 2 || isempty(manifest_file)
    manifest_file = fullfile(here, 'manifests', sprintf('recovery_single_%s.csv', conv));
end
if nargin < 3, rec_dir = fullfile(here, 'outputs_recovery'); end
[~, manifest_name] = fileparts(manifest_file);
outdir = fullfile(rec_dir, ['analysis_' manifest_name]);
if ~isfolder(outdir), mkdir(outdir); end

model_names = ["mlp", "sr", "gp", "lr"];
model_labels = ["MLP", "SR", "GP", "LR"];
model_colors = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"];
markers = {'o', 's', 'd', '^'};
param_names = ["gamma1", "lambda1", "delta", "epsilon"];
param_labels = {'\gamma_1', '\lambda_1', '\delta', '\epsilon'};
param_log = [true, true, false, false];
param_range = [5e-5 5e-2; 0.01 0.1; 10 100; 10 100];

T = readtable(fullfile(here, 'manifests', 'recovery_truths.csv'), 'TextType', 'string');
T = T(T.convention == conv, :);
C = readtable(manifest_file, 'TextType', 'string', 'Delimiter', ',');
model_keep = ismember(model_names, unique(C.method));
model_names = model_names(model_keep);  model_labels = model_labels(model_keep);
model_colors = model_colors(model_keep);  markers = markers(model_keep);
n_truths = height(T);  n_models = numel(model_names);  n_params = numel(param_names);
same = @(a, b) abs(a - b) <= 1e-6 * abs(b);

%% truths
truth_N = cell(n_truths, 1);
for t = 1:n_truths
    f = fullfile(here, T.truth_dir(t), T.truth_case_id(t) + ".mat");
    if strcmp(T.truth_dir(t), "outputs_recovery"), f = fullfile(rec_dir, T.truth_case_id(t) + ".mat"); end
    if isfile(f)
        s = load(f, 'N');
        truth_N{t} = s.N;
    else
        fprintf('truth %s missing (%s)\n', T.truth_id(t), f);
    end
end

%% candidate losses
% every manifest row's N is loaded once; NaN when the run does not exist yet
cand_N = cell(height(C), 1);
for r = 1:height(C)
    f = fullfile(rec_dir, C.case_id(r) + ".mat");
    if isfile(f)
        s = load(f, 'N');
        cand_N{r} = s.N;
    end
end
fprintf('%d of %d simulations found\n', sum(~cellfun(@isempty, cand_N)), height(C));
nmse = @(S, N) sum((S(:) - N(:)).^2) / sum(N(:).^2);

rn = struct('conv', conv, 'manifest', string(manifest_file), 'truth_ids', T.truth_id, ...
    'true_params', [T.gamma1, T.lambda1, T.delta, T.epsilon], ...
    'model_labels', model_labels, 'model_colors', model_colors, ...
    'param_names', param_names, 'param_labels', {param_labels}, 'param_log', param_log);
rn.grid = cell(n_params, n_models, n_truths);      % candidate values
rn.loss = cell(n_params, n_models, n_truths);      % NMSE vs numeric truth
rn.control = cell(n_params, n_models, n_truths);   % NMSE vs surrogate at truth
[rn.hat, rn.min_loss, rn.flat, rn.hat_control] = deal(nan(n_params, n_models, n_truths));
summary = table();
for t = 1:n_truths
    theta = rn.true_params(t, :);
    for m = 1:n_models
        for i = 1:n_params
            others = setdiff(1:n_params, i);
            rows = C.method == model_names(m);
            for j = others
                rows = rows & same(C.(param_names(j)), theta(j));
            end
            rows = find(rows);
            [x, order] = sort(C.(param_names(i))(rows));
            rows = rows(order);
            at_truth = rows(same(x, theta(i)));
            [loss, control] = deal(nan(size(x)));
            for k = 1:numel(rows)
                S = cand_N{rows(k)};
                if isempty(S), continue, end
                if ~isempty(truth_N{t}), loss(k) = nmse(S, truth_N{t}); end
                if ~isempty(at_truth) && ~isempty(cand_N{at_truth(1)})
                    control(k) = nmse(S, cand_N{at_truth(1)});
                end
            end
            rn.grid{i,m,t} = x;  rn.loss{i,m,t} = loss;  rn.control{i,m,t} = control;
            [rn.hat(i,m,t), rn.min_loss(i,m,t), rn.flat(i,m,t)] = gridMinimum(x, loss, param_log(i));
            rn.hat_control(i,m,t) = gridMinimum(x, control, param_log(i));
            if param_log(i), err = log10(rn.hat(i,m,t) / theta(i)); else, err = rn.hat(i,m,t) - theta(i); end
            summary = [summary; table(model_names(m), T.truth_id(t), param_names(i), theta(i), ...
                rn.hat(i,m,t), rn.min_loss(i,m,t), rn.flat(i,m,t), err, rn.hat_control(i,m,t), ...
                sum(isfinite(loss)), numel(x), 'VariableNames', {'model', 'truth_id', 'param', ...
                'true_value', 'recovered_value', 'min_loss', 'flat', 'error', 'recovered_control', ...
                'n_found', 'n_grid'})]; %#ok<AGROW>
        end
    end
end
writetable(summary, fullfile(outdir, 'recovery_summary.csv'));
save(fullfile(outdir, 'recovery_network.mat'), 'rn');

%% loss curves, one figure per truth
for t = 1:n_truths
    if all(cellfun(@(l) all(isnan(l)), rn.loss(:,:,t)) & cellfun(@(l) all(isnan(l)), rn.control(:,:,t)), 'all')
        continue % nothing run for this truth yet
    end
    theta = rn.true_params(t, :);
    fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 max(10*n_models + 2, 34) 8.5*n_params]);
    tl = tiledlayout(fig, n_params, n_models, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, sprintf('%s (%s):  \\gamma_1 = %.3g,  \\lambda_1 = %.3g,  \\delta = %.3g,  \\epsilon = %.3g', ...
        strrep(T.truth_id(t), '_', '\_'), conv, theta), 'FontSize', 14);
    for i = 1:n_params
        for m = 1:n_models
            ax = nexttile(tl, (i-1)*n_models + m); hold(ax, 'on');
            x = rn.grid{i,m,t};
            loss = rn.loss{i,m,t};  control = rn.control{i,m,t};
            % the control is 0 at the truth (the run against itself), which a log axis
            % cannot show: draw it at the bottom of the axis instead
            positive = [loss(loss > 0); control(control > 0)];
            y_floor = 10^floor(log10(min(positive)) - 0.5);
            control(control == 0) = y_floor;
            xline(ax, theta(i), '--k', 'LineWidth', 1);
            plot(ax, x, control, '-', 'Color', [0.55 0.55 0.55], 'LineWidth', 1);
            plot(ax, x, loss, '-o', 'Color', hex2rgb(model_colors(m)), 'LineWidth', 2, ...
                'MarkerSize', 4, 'MarkerFaceColor', hex2rgb(model_colors(m)));
            control_txt = '';
            if ~isnan(rn.hat_control(i,m,t)), control_txt = sprintf('  (control %.3g)', rn.hat_control(i,m,t)); end
            if rn.flat(i,m,t)
                title(ax, sprintf('%s:  %s flat (not identifiable)', model_labels(m), param_labels{i}), 'FontWeight', 'normal');
            elseif ~isnan(rn.hat(i,m,t))
                plot(ax, rn.hat(i,m,t), rn.min_loss(i,m,t), 'o', 'MarkerFaceColor', 'r', ...
                    'MarkerEdgeColor', 'k', 'MarkerSize', 7);
                title(ax, sprintf('%s:  %s = %.3g%s', model_labels(m), param_labels{i}, rn.hat(i,m,t), control_txt), ...
                    'FontWeight', 'normal');
            else
                title(ax, model_labels(m), 'FontWeight', 'normal');
            end
            if param_log(i), set(ax, 'XScale', 'log'); end
            set(ax, 'YScale', 'log', 'Box', 'off', 'TickDir', 'out');
            xlim(ax, param_range(i,:));
            ylim(ax, [y_floor, Inf]);
            grid(ax, 'on');
            xlabel(ax, param_labels{i});
            if m == 1, ylabel(ax, 'NMSE'); end
            hold(ax, 'off');
        end
    end
    ax1 = nexttile(tl, 1); hold(ax1, 'on');
    h = [plot(ax1, NaN, NaN, '-o', 'Color', [0.3 0.3 0.3], 'LineWidth', 2), ...
         plot(ax1, NaN, NaN, '-', 'Color', [0.55 0.55 0.55], 'LineWidth', 1), ...
         plot(ax1, NaN, NaN, '--k', 'LineWidth', 1), ...
         plot(ax1, NaN, NaN, 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k')];
    lg = legend(h, {'NMSE vs numeric truth', 'NMSE vs surrogate at truth (control; 0 at the truth drawn at the axis bottom)', ...
        'true value', 'recovered (minimum NMSE)'}, ...
        'Orientation', 'horizontal', 'Box', 'off', 'FontSize', 11);
    lg.Layout.Tile = 'south';
    hold(ax1, 'off');
    exportgraphics(fig, fullfile(outdir, sprintf('single_%s.png', T.truth_id(t))), 'Resolution', 130);
    close(fig);
end

%% true vs recovered, all truths
fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 10*n_params + 2 11]);
tl = tiledlayout(fig, 1, n_params, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('Network recovery (%s): recovered vs true', conv), 'FontSize', 14);
for i = 1:n_params
    ax = nexttile(tl); hold(ax, 'on');
    plot(ax, param_range(i,:), param_range(i,:), '-k', 'LineWidth', 1, 'DisplayName', 'recovered = true');
    for m = 1:n_models
        plot(ax, rn.true_params(:,i), squeeze(rn.hat(i,m,:)), markers{m}, 'MarkerSize', 8, ...
            'MarkerFaceColor', hex2rgb(model_colors(m)), 'MarkerEdgeColor', 'w', 'DisplayName', model_labels(m));
    end
    if param_log(i), set(ax, 'XScale', 'log', 'YScale', 'log'); end
    set(ax, 'Box', 'off', 'TickDir', 'out');
    xlim(ax, param_range(i,:));  ylim(ax, param_range(i,:));
    axis(ax, 'square');  grid(ax, 'on');
    xlabel(ax, ['true ' param_labels{i}]);  ylabel(ax, ['recovered ' param_labels{i}]);
    if i == 1, legend(ax, 'Location', 'northwest', 'Box', 'off'); end
    hold(ax, 'off');
end
exportgraphics(fig, fullfile(outdir, 'true_vs_recovered.png'), 'Resolution', 130);
close(fig);
fprintf('Wrote %s\n', outdir);
end

function [hat, min_loss, flat] = gridMinimum(x, loss, is_log)
% Grid minimum of loss(x), refined by a parabola through it and its two
% neighbours (in log10 x for log-searched parameters), kept between the neighbours.
[hat, min_loss, flat] = deal(NaN, NaN, false);
ok = isfinite(loss);
if sum(ok) < 3, return, end
x = x(ok);  loss = loss(ok);
if max(loss) - min(loss) <= 1e-12 * max(abs(loss))
    flat = true;
    return
end
u = x;
if is_log, u = log10(x); end
[min_loss, k] = min(loss);
hat_u = u(k);
if k > 1 && k < numel(u)
    p = polyfit(u(k-1:k+1), loss(k-1:k+1), 2);
    if p(1) > 0
        hat_u = min(max(-p(2) / (2*p(1)), u(k-1)), u(k+1));
        min_loss = max(polyval(p, hat_u), 0); % NMSE >= 0
    end
end
hat = hat_u;
if is_log, hat = 10^hat_u; end
end

function rgb = hex2rgb(hex)
hex = char(hex);
rgb = hex2dec({hex(2:3), hex(4:5), hex(6:7)})' / 255;
end

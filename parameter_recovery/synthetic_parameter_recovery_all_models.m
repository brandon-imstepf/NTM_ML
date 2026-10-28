% Synthetic Parameter Recovery - All Models %
% 10/16/26, outputs added 10/20/26 %
% Recovers every parameter alone (single-parameter) and every unique pair of
% parameters (joint) from the locked single-edge sets, with every surrogate,
% using one shared loss.
%
% One "recovery" = one locked set: gamma1, lambda1, delta, epsilon are fixed
% and the set has many (N1, N2, flux, W1) points. The loss for a candidate
% parameter value is the mean, over all points of the set, of the squared
% error of a scaled output. It is computed for each output separately and
% for both together:
%     flux:  loss = mean( (flux_pred_scaled - flux_scaled).^2 )      (phi_tau)
%     w1:    loss = mean( (w1_pred_scaled - w1_scaled).^2 )          (Delta_tau)
%     both:  loss = (flux loss + w1 loss) / 2
% Scaled = z-scored with the SR target scaling (flux_scaling_params.mat,
% emc_scaling_params.mat), as in the GP and flux-only scripts, so losses are
% on the same scale for every model.
%   Single p:     the other three parameters held at their true values.
%   Joint (p, q): the other two parameters held at their true values.
% A parameter the model's flux does not depend on (e.g. lambda1 and epsilon
% for SR) leaves the loss flat along it. It is flagged "flat" and its
% recovered value is NaN (not identifiable).
%
% The models come from network_validation/make_surrogate.m, the same code the
% network validation uses, so a change there changes both pipelines.
%
% Outputs (in outdir/<output>/, output = flux, w1, both):
%   recovery_summary.csv                one row per model x set x recovery x recovered parameter
%   true_vs_recovered_single.png        true vs recovered, single-parameter, all models
%   true_vs_recovered_joint.png         true vs recovered, each joint pair, all models
%   by_set/single_<set>.png             per set: loss curve for each parameter (rows) x model (columns)
%   by_set/joint_<set>.png              per set: loss surface for each pair (rows) x model (columns)
%   by_group/single_mean_<group>.png    per group of sets (interior, edge, ood): the same,
%   by_group/joint_mean_<group>.png     averaged relative to the truth
%   recovery_curves.mat                 the curves and surfaces behind by_set/ and by_group/,
%                                       so the figures can be redrawn without rerunning:
%                                       load it, then recoveryOverviewPlot(rc, s, outdir)
%                                       or recoveryGroupMeanPlot(rc, group, outdir)
clc;
clear;
close all;

%% 1. Settings (change here; affects every model)
model_names  = ["mlp", "sr", "gp", "lr"];   % make_surrogate method names
model_labels = ["MLP", "SR", "GP", "LR"];   % names used in plots and the table
model_colors = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"];
model_markers = {'o', 's', 'd', '^'};

locked_dir = fullfile('locked_edge', 'outputs'); % locked sets, one CSV per set
outdir = 'recovery_all_models';

% outputs the loss is computed on (one full recovery each)
outputs = struct( ...
    'name',  {"flux",       "w1",          "both"}, ...
    'label', {'\phi_\tau', '\Delta_\tau', '\phi_\tau and \Delta_\tau'});
% outputs to run: flux and w1, each on its own ("both" was run on 10/20/26 and is
% kept in recovery_all_models/both/, but recovery is reported per output)
outputs = outputs(ismember([outputs.name], ["flux", "w1"]));

ncores = 8;            % parfor workers (one locked set each); 0 = one set at a time
numiterations = 200;   % points on each single-parameter loss curve (plot only)
gridN = 30;            % joint grid points per axis (grid search, then fminsearch polish)

% Parameters to recover, with their training ranges.
%   col: column in the set = argument position in flux_fun
%   log: true = search in log10, false = search on a linear scale
params = struct( ...
    'name',  {"gamma1",       "lambda1",     "delta",     "epsilon"}, ...
    'label', {'\gamma_1',     '\lambda_1',   '\delta',    '\epsilon'}, ...
    'col',   {1,              2,             3,           4}, ...
    'range', {[5e-5, 5e-2],   [0.01, 0.1],   [10, 100],   [10, 100]}, ...
    'log',   {true,           true,          false,       false});
% joint recoveries: every unique pair of params (row k = [a, b], indices into params)
pairs = nchoosek(1:numel(params), 2);

%% 2. Load the models
addpath(fullfile('..', 'network_validation'));
cfg = sntm_config();
surrogates = cell(numel(model_names), 1);
for m = 1:numel(model_names)
    surrogates{m} = make_surrogate(char(model_names(m)), cfg);
end

% load matching z-score mean used to fit (SR target scaling, shared by all models)
scaling.flux = load(fullfile("flux_scaling_params.mat"));
scaling.w1 = load(fullfile("emc_scaling_params.mat"));

%% 3. Load the locked sets
% columns: [gamma1, lambda1, delta, epsilon, N1, N2, FValue, W1]
set_files = dir(fullfile(locked_dir, '*.csv'));
if isempty(set_files)
    error('No locked sets found in %s', locked_dir);
end
datasets = cell(numel(set_files), 1);
set_ids = strings(numel(set_files), 1);
for s = 1:numel(set_files)
    datasets{s} = readmatrix(fullfile(set_files(s).folder, set_files(s).name));
    [~, set_ids(s)] = fileparts(set_files(s).name);
    % drop points the generator could not solve (written as NaN)
    datasets{s} = datasets{s}(~isnan(datasets{s}(:,7)) & ~isnan(datasets{s}(:,8)), :);
end
fprintf('%d locked sets: %s\n', numel(set_ids), strjoin(set_ids, ', '));

%% 4. Search grids (identical every run)
for i = 1:numel(params)
    params(i).curve_grid = searchSpace(params(i), numiterations);
    params(i).joint_grid = searchSpace(params(i), gridN);
end

%% 5. Recovery: every model x every locked set (recoverOneSet, below)
n_models = numel(model_names);
n_sets = numel(datasets);
n_params = numel(params);
n_pairs = size(pairs, 1);
if ncores > 0 && isempty(gcp('nocreate'))
    parpool(ncores);
end

set_groups = unique(string(regexp(set_ids, '^[^_]+', 'match', 'once')), 'stable');

for o = 1:numel(outputs)
    output = outputs(o);
    out_o = fullfile(outdir, output.name);
    if ~isfolder(out_o), mkdir(out_o); end

    % everything the figures need, saved in step 8 (model x set)
    rc = struct();
    rc.output = output.name;
    rc.output_label = output.label;
    rc.model_labels = model_labels;
    rc.model_colors = model_colors;
    rc.set_ids = set_ids;
    rc.params = params;
    rc.pairs = pairs;
    rc.true_params = zeros(n_sets, n_params);                       % [gamma1, lambda1, delta, epsilon] per set
    rc.single_curves = zeros(numiterations, n_params, n_models, n_sets);
    rc.single_hat = zeros(n_params, n_models, n_sets);               % recovered value (NaN if flat)
    rc.single_loss = zeros(n_params, n_models, n_sets);              % loss at the recovered value
    rc.single_flat = false(n_params, n_models, n_sets);
    rc.joint_surfaces = zeros(gridN, gridN, n_pairs, n_models, n_sets);
    rc.joint_hat = zeros(n_pairs, 2, n_models, n_sets);              % recovered [a, b] (NaN if flat)
    rc.joint_loss = zeros(n_pairs, n_models, n_sets);
    rc.joint_flat = false(n_pairs, 2, n_models, n_sets);             % [a flat, b flat]

    for m = 1:n_models
        sur = surrogates{m};
        fprintf('%s, %s: recovering %d sets\n', output.name, model_labels(m), n_sets);
        out = cell(n_sets, 1);
        parfor (s = 1:n_sets, ncores)
            out{s} = recoverOneSet(sur, scaling, datasets{s}, params, pairs, output.name);
            fprintf('[%s | %s] done\n', model_labels(m), set_ids(s));
        end
        for s = 1:n_sets
            rc.true_params(s,:) = out{s}.theta_true;
            rc.single_curves(:,:,m,s) = out{s}.single_curves;
            rc.single_hat(:,m,s) = out{s}.single_hat;
            rc.single_loss(:,m,s) = out{s}.single_loss;
            rc.single_flat(:,m,s) = out{s}.single_flat;
            rc.joint_surfaces(:,:,:,m,s) = out{s}.joint_surfaces;
            rc.joint_hat(:,:,m,s) = out{s}.joint_hat;
            rc.joint_loss(:,m,s) = out{s}.joint_loss;
            rc.joint_flat(:,:,m,s) = out{s}.joint_flat;
        end
    end

    %% 6. Summary table (one row per recovered parameter)
    results = {};
    for m = 1:n_models
        for s = 1:n_sets
            for i = 1:n_params
                results(end+1,:) = {model_labels(m), set_ids(s), "single", params(i).name, params(i).name, ...
                    rc.true_params(s,i), rc.single_hat(i,m,s), rc.single_loss(i,m,s), rc.single_flat(i,m,s)}; %#ok<SAGROW>
            end
            for k = 1:n_pairs
                recovery = params(pairs(k,1)).name + "_" + params(pairs(k,2)).name;
                for ab = 1:2
                    i = pairs(k,ab);
                    results(end+1,:) = {model_labels(m), set_ids(s), "joint", recovery, params(i).name, ...
                        rc.true_params(s,i), rc.joint_hat(k,ab,m,s), rc.joint_loss(k,m,s), rc.joint_flat(k,ab,m,s)}; %#ok<SAGROW>
                end
            end
        end
    end
    summary = cell2table(results, 'VariableNames', {'model', 'set_id', 'type', 'recovery', 'param', ...
        'true_value', 'recovered_value', 'min_loss', 'flat'});
    summary.log_error = abs(log10(summary.recovered_value) - log10(summary.true_value));
    % at_bound: the recovered value sits on the edge of the search range, so the
    % minimum may lie outside it (or the true value is on the edge itself)
    summary.at_bound = false(height(summary), 1);
    for i = 1:n_params
        rows = summary.param == params(i).name;
        u = toSearch(params(i), summary.recovered_value(rows));
        lo = toSearch(params(i), params(i).range(1));
        hi = toSearch(params(i), params(i).range(2));
        summary.at_bound(rows) = abs(u - lo) < 1e-3 * (hi - lo) | abs(u - hi) < 1e-3 * (hi - lo);
    end
    writetable(summary, fullfile(out_o, 'recovery_summary.csv'));

    %% 7. True vs recovered, all models
    trueVsRecoveredPlot(summary, model_labels, model_colors, model_markers, params, pairs, out_o);

    %% 8. Two figures per set, all models (recoveryOverviewPlot.m)
    save(fullfile(out_o, 'recovery_curves.mat'), 'rc');
    for s = 1:n_sets
        recoveryOverviewPlot(rc, s, out_o);
    end

    %% 9. Mean figures per group of sets (recoveryGroupMeanPlot.m)
    % group = start of the set_id up to the first underscore: interior, edge, ood, ...
    % groups with one set (e.g. paper_fig8) are skipped
    for g = 1:numel(set_groups)
        recoveryGroupMeanPlot(rc, set_groups(g), out_o);
    end
end % outputs


%% Local functions

function out = recoverOneSet(sur, scaling, dataset, params, pairs, output)
% Every single-parameter and joint recovery for one model, one locked set and
% one output ("flux", "w1" or "both")

    theta_true = dataset(1,1:4);   % [gamma1, lambda1, delta, epsilon]
    n1 = dataset(:,5);
    n2 = dataset(:,6);
    fs = scaling.flux;
    ws = scaling.w1;
    flux_scaled = ( dataset(:,7) - fs.target_mean ) / fs.target_std;
    w1_scaled = ( dataset(:,8) - ws.target_mean ) / ws.target_std;

    % surrogate predictions, scaled like the data
    flux_pred_scaled = @(theta) ...
        (sur.flux_fun(theta(1),theta(2),theta(3),theta(4),n1,n2) - fs.target_mean) / fs.target_std;
    w1_pred_scaled = @(theta) ...
        (sur.w1_fun(theta(1),theta(2),theta(3),theta(4),n1,n2,[]) - ws.target_mean) / ws.target_std;
    % loss for a candidate theta = [gamma1, lambda1, delta, epsilon], over all points of the set
    flux_loss = @(theta) mean( (flux_pred_scaled(theta) - flux_scaled).^2 );
    w1_loss = @(theta) mean( (w1_pred_scaled(theta) - w1_scaled).^2 );
    switch output
        case "flux", loss_fun = flux_loss;
        case "w1",   loss_fun = w1_loss;
        case "both", loss_fun = @(theta) (flux_loss(theta) + w1_loss(theta)) / 2;
    end

    n_params = numel(params);
    n_pairs = size(pairs, 1);
    out.theta_true = theta_true;

    % single-parameter: p alone, the rest at their true values
    out.single_curves = zeros(numel(params(1).curve_grid), n_params);
    out.single_hat = zeros(n_params, 1);
    out.single_loss = zeros(n_params, 1);
    out.single_flat = false(n_params, 1);
    for i = 1:n_params
        p = params(i);
        single_loss = @(x) loss_fun(setParams(theta_true, p.col, x));
        curve = arrayfun(single_loss, p.curve_grid);
        flat = isFlat(curve, 1);
        u_hat = fminbnd(@(u) single_loss(fromSearch(p, u)), ...
            toSearch(p, p.range(1)), toSearch(p, p.range(2)));
        out.single_curves(:,i) = curve;
        out.single_loss(i) = single_loss(fromSearch(p, u_hat));
        out.single_flat(i) = flat;
        out.single_hat(i) = fromSearch(p, u_hat);
        if flat, out.single_hat(i) = NaN; end
    end

    % joint: every pair (a, b), the other two at their true values
    gridN = numel(params(1).joint_grid);
    out.joint_surfaces = zeros(gridN, gridN, n_pairs);
    out.joint_hat = zeros(n_pairs, 2);
    out.joint_loss = zeros(n_pairs, 1);
    out.joint_flat = false(n_pairs, 2);
    for k = 1:n_pairs
        pa = params(pairs(k,1));
        pb = params(pairs(k,2));
        joint_loss = @(x, y) loss_fun(setParams(theta_true, [pa.col, pb.col], [x, y]));

        % grid search: a down the rows, b across the columns
        [A, B] = ndgrid(pa.joint_grid, pb.joint_grid);
        loss_surface = zeros(size(A));
        for j = 1:numel(A)
            loss_surface(j) = joint_loss(A(j), B(j));
        end
        [~, lin_idx] = min(loss_surface(:));
        [ai, bi] = ind2sub(size(loss_surface), lin_idx);

        % polish the best grid point with fminsearch in search coordinates, clamped to the range
        polish_fun = @(u) joint_loss(fromSearch(pa, u(1)), fromSearch(pb, u(2)));
        u_hat = fminsearch(polish_fun, [toSearch(pa, pa.joint_grid(ai)), toSearch(pb, pb.joint_grid(bi))]);
        u_hat(1) = min(max(u_hat(1), toSearch(pa, pa.range(1))), toSearch(pa, pa.range(2)));
        u_hat(2) = min(max(u_hat(2), toSearch(pb, pb.range(1))), toSearch(pb, pb.range(2)));
        hat = [fromSearch(pa, u_hat(1)), fromSearch(pb, u_hat(2))];

        % flat along a = the loss does not change down any column, and likewise for b
        flat = [isFlat(loss_surface, 1), isFlat(loss_surface, 2)];
        out.joint_surfaces(:,:,k) = loss_surface;
        out.joint_loss(k) = joint_loss(hat(1), hat(2)); % loss at the clamped point
        hat(flat) = NaN;
        out.joint_hat(k,:) = hat;
        out.joint_flat(k,:) = flat;
    end
end

function flat = isFlat(values, dim)
% true if values do not change at all along dimension dim (to 1e-12 relative)
    if dim == 1
        ref = values(1,:);
    else
        ref = values(:,1);
    end
    flat = all(abs(values - ref) <= 1e-12 * abs(ref), 'all');
end

function theta = setParams(theta, cols, values)
% theta with the entries at cols replaced by values
    theta(cols) = values;
end

function grid = searchSpace(p, n)
% n points evenly spaced in search coordinates over p.range
    grid = fromSearch(p, linspace(toSearch(p, p.range(1)), toSearch(p, p.range(2)), n)');
end

function u = toSearch(p, x)
% value -> search coordinate (log10 for log-searched parameters)
    if p.log, u = log10(x); else, u = x; end
end

function x = fromSearch(p, u)
% search coordinate -> value
    if p.log, x = 10.^u; else, x = u; end
end

function trueVsRecoveredPlot(summary, model_labels, model_colors, model_markers, params, pairs, outdir)
% y = x is perfect recovery. Flat (not identifiable) recoveries are NaN, so they do not appear.
%   true_vs_recovered_single.png: one panel per parameter
%   true_vs_recovered_joint.png:  one row per pair, one panel per parameter of the pair

    % single
    fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 10*numel(params) 10]);
    tl = tiledlayout(fig, 1, numel(params), 'TileSpacing', 'compact', 'Padding', 'compact');
    for i = 1:numel(params)
        ax = nexttile(tl);
        rows = summary.type == "single" & summary.param == params(i).name;
        scatterPanel(ax, summary(rows,:), params(i), model_labels, model_colors, model_markers);
        title(ax, [params(i).label ', single-parameter'], 'FontWeight', 'normal');
    end
    lg = legend(nexttile(tl, 1), 'Orientation', 'horizontal', 'Box', 'off'); lg.Layout.Tile = 'north';
    exportgraphics(fig, fullfile(outdir, 'true_vs_recovered_single.png'), 'Resolution', 200);
    close(fig);

    % joint
    fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 20 9*size(pairs,1)]);
    tl = tiledlayout(fig, size(pairs,1), 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    for k = 1:size(pairs,1)
        recovery = params(pairs(k,1)).name + "_" + params(pairs(k,2)).name;
        for ab = 1:2
            p = params(pairs(k,ab));
            other = params(pairs(k,3-ab));
            ax = nexttile(tl);
            rows = summary.type == "joint" & summary.recovery == recovery & summary.param == p.name;
            scatterPanel(ax, summary(rows,:), p, model_labels, model_colors, model_markers);
            title(ax, [p.label ', joint with ' other.label], 'FontWeight', 'normal');
        end
    end
    lg = legend(nexttile(tl, 1), 'Orientation', 'horizontal', 'Box', 'off'); lg.Layout.Tile = 'north';
    exportgraphics(fig, fullfile(outdir, 'true_vs_recovered_joint.png'), 'Resolution', 200);
    close(fig);
end

function scatterPanel(ax, rows, p, model_labels, model_colors, model_markers)
% true vs recovered for one parameter, one marker style per model
    hold(ax, 'on');
    plot(ax, p.range, p.range, ':', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
    for m = 1:numel(model_labels)
        r = rows(rows.model == model_labels(m), :);
        scatter(ax, r.true_value, r.recovered_value, 60, hex2rgb(model_colors(m)), 'filled', ...
            model_markers{m}, 'MarkerEdgeColor', 'w', 'DisplayName', model_labels(m));
    end
    if p.log
        set(ax, 'XScale', 'log', 'YScale', 'log');
    end
    set(ax, 'XLim', p.range, 'YLim', p.range, 'Box', 'off', 'TickDir', 'out');
    axis(ax, 'square'); grid(ax, 'on');
    xlabel(ax, 'true'); ylabel(ax, 'recovered');
    hold(ax, 'off');
end

function rgb = hex2rgb(hex)
    h = char(hex);
    rgb = hex2dec({h(2:3), h(4:5), h(6:7)})' / 255;
end

function D = fig8_offset_diagnostics()
% Why the MLP and GP sNTMs of the new Fig. 8 runs (outputs_fig8/, gamma1 5e-5,
% lambda 0.01, delta 10, epsilon 50, ECl seed, init_rescale 0.01) hold a
% constant offset in regional tau (MLP about 0.89 x, GP about 1.23 x numeric).
%
% 1. Error map along the numeric run: at every time point, the MLP and GP
%    phi_tau and Delta_tau (W1) at the edge inputs of the numeric run (active
%    edges: C > 0 and N_i > 0 or N_j > 0), against the numeric F and W1 at the
%    same inputs. NMSE (Eq. 7) over all active edges, over the edges of the seed
%    region and over the others; the inputs' position in the training data (N
%    above the 99.9th percentile or the maximum of the MLP training N).
%    W1 is not computed by the simulator at the last time point.
% 2. Seed-region close-up: N of the seed and its 6 most strongly connected
%    neighbours over days 0-10, and the mass budget: node mass sum_i V_i (N_i + M_i)
%    of each run, and the numeric run's stored edge mass (sum of EdgeMass) and its
%    release, edge mass(0) - edge mass(t), against the numeric node-mass gain.
% Writes results/fig8_offset_errors_by_time.csv, results/fig8_offset_seed_edges_t0.csv,
% results/fig8_offset_mass.csv and outputs_fig8/figures/offset_*.png.
%
%   D = fig8_offset_diagnostics();

here = fileparts(mfilename('fullpath'));
cfg = sntm_config();
fdir = fullfile(here, 'outputs_fig8', 'figures');
if ~isfolder(fdir), mkdir(fdir); end
file = @(m) fullfile(here, 'outputs_fig8', "fig8__" + m + "__ECl__g5e-05_l0.01_d10_e50.mat");
num = load(file("numeric"));
run = struct('mlp', load(file("mlp")), 'gp', load(file("gp")));
p = num.meta;
day = num.trange(:)' * 180;
nt = numel(day);
[~, seed] = max(num.N(:, 1));
fprintf('seed region %d (%s), N0 = %.3g uM\n', seed, num.region_names{seed}, num.N(seed, 1));

% training N (MLP training data)
Ntrain = [];
for f = dir(fullfile(cfg.ntm_ml_root, 'data', 'data_neural_network_model', 'training_data_e6', 'data_bias_e5_*.csv'))'
    T = readmatrix(fullfile(f.folder, f.name));
    Ntrain = [Ntrain; T(:, 5); T(:, 6)]; %#ok<AGROW>
end
pct = prctile(Ntrain, [99 99.9]);  Nmax = max(Ntrain);
fprintf('training N: 99th %.3g, 99.9th %.3g, max %.3g uM\n', pct, Nmax);

%% 1. Error map
sur = struct('mlp', make_surrogate('mlp', cfg), 'gp', make_surrogate('gp', cfg));
[ei, ej] = find(num.C > 0);
seed_edge = ei == seed | ej == seed;
rows = {};
for k = 1:nt
    a = num.N(ei, k);  b = num.N(ej, k);
    on = a > 0 | b > 0;
    Fn = arrayfun(@(q) num.F(ei(q), ej(q), k), (1:numel(ei))');
    Wn = arrayfun(@(q) num.W1(ei(q), ej(q), k), (1:numel(ei))');
    x = [a(on); b(on)];
    row = {day(k), nnz(on), nnz(on & seed_edge), max(x) / Nmax, mean(x > pct(2)), mean(x > Nmax)};
    for m = ["mlp", "gp"]
        Fs = sur.(m).flux_fun(p.gamma1, p.lambda1, p.delta, p.epsilon, a, b);
        Ws = sur.(m).w1_fun(p.gamma1, p.lambda1, p.delta, p.epsilon, a, b, []);
        sel = {on, on & seed_edge, on & ~seed_edge};
        for c = 1:3
            row{end+1} = nmse_eq7(Fs(sel{c}), Fn(sel{c})); %#ok<AGROW>
            if k < nt, row{end+1} = nmse_eq7(Ws(sel{c}), Wn(sel{c})); else, row{end+1} = NaN; end %#ok<AGROW>
        end
        % summed W1 over active edges, surrogate / numeric (the edge-mass coefficient level)
        if k < nt, row{end+1} = sum(Ws(on)) / sum(Wn(on)); else, row{end+1} = NaN; end %#ok<AGROW>
    end
    rows(end+1, :) = row; %#ok<AGROW>
end
names = {'day', 'n_active', 'n_active_seed', 'maxN_over_train_max', 'frac_N_above_p999', 'frac_N_above_max'};
for m = ["mlp", "gp"]
    for c = ["all", "seed", "other"]
        names{end+1} = char(m + "_flux_nmse_" + c); %#ok<AGROW>
        names{end+1} = char(m + "_w1_nmse_" + c); %#ok<AGROW>
    end
    names{end+1} = char(m + "_w1_sum_ratio"); %#ok<AGROW>
end
D.time = cell2table(rows, 'VariableNames', names);
writetable(D.time, fullfile(here, 'results', 'fig8_offset_errors_by_time.csv'));

% seed edges at t0, one row per edge
k = 1;  q = find(seed_edge & (num.N(ei, k) > 0 | num.N(ej, k) > 0));
a = num.N(ei(q), k);  b = num.N(ej(q), k);
E = table(ei(q), ej(q), a, b, 'VariableNames', {'i', 'j', 'N_i', 'N_j'});
E.F_numeric = arrayfun(@(r) num.F(ei(r), ej(r), k), q);
E.W1_numeric = arrayfun(@(r) num.W1(ei(r), ej(r), k), q);
for m = ["mlp", "gp"]
    E.("F_" + m) = sur.(m).flux_fun(p.gamma1, p.lambda1, p.delta, p.epsilon, a, b);
    E.("W1_" + m) = sur.(m).w1_fun(p.gamma1, p.lambda1, p.delta, p.epsilon, a, b, []);
end
D.seed_edges_t0 = E;
writetable(E, fullfile(here, 'results', 'fig8_offset_seed_edges_t0.csv'));

%% 2. Mass budget and seed close-up
load(fullfile(cfg.matdir, 'DefaultAtlas.mat'), 'DefaultAtlas');
load(fullfile(cfg.matdir, 'CCF_labels.mat'), 'CCF_labels');
inds = ismember(CCF_labels(:,3), 'Hippocampus') | ismember(CCF_labels(:,1), 'Piriform area') | ...
    ismember(CCF_labels(:,3), 'Retrosplenial Area');
V = DefaultAtlas.volumes(inds);  V = V(:);
node = @(r) sum(V .* (r.N + r.M), 1);
em_num = squeeze(sum(sum(num.EdgeMass, 1, 'omitnan'), 2))';
M = table(day', node(num)', node(run.mlp)', node(run.gp)', em_num', 'VariableNames', ...
    {'day', 'node_mass_numeric', 'node_mass_mlp', 'node_mass_gp', 'edge_mass_numeric'});
D.mass = M;
writetable(M, fullfile(here, 'results', 'fig8_offset_mass.csv'));
fprintf('node mass at day 0: %.4g (all runs: %.4g, %.4g)\n', M.node_mass_numeric(1), M.node_mass_mlp(1), M.node_mass_gp(1));
fprintf('numeric edge mass: day 0 %.4g, day 150 %.4g; edge release %.4g vs node gain %.4g\n', em_num(1), ...
    em_num(find(day >= 150, 1)), em_num(1) - em_num(find(day >= 150, 1)), ...
    M.node_mass_numeric(find(day >= 150, 1)) - M.node_mass_numeric(1));

%% Figures
col = struct('num', [0 0 0], 'mlp', [0.17 0.40 0.69], 'gp', [0.11 0.62 0.47]);
fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [1 1 26 16], 'Color', 'w');
tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
dd = max(day, 0.1);   % log axis
for c = 1:2
    tgt = ["flux", "w1"];  tl_ = ["\phi_\tau", "\Delta_\tau (W_1)"];
    ax = nexttile(tl);  hold(ax, 'on');  box(ax, 'on');  grid(ax, 'on');
    for m = ["mlp", "gp"]
        plot(ax, dd, D.time.(m + "_" + tgt(c) + "_nmse_seed"), '-', 'Color', col.(m), 'LineWidth', 1.5);
        plot(ax, dd, D.time.(m + "_" + tgt(c) + "_nmse_other"), '--', 'Color', col.(m), 'LineWidth', 1.2);
    end
    set(ax, 'XScale', 'log', 'YScale', 'log');  xlim(ax, [0.1 180]);
    title(ax, "NMSE of " + tl_(c) + " at the numeric run's edge inputs", 'FontWeight', 'normal');
    xlabel(ax, 'day (day 0 drawn at 0.1)');  ylabel(ax, 'NMSE');
    legend(ax, {'MLP, seed edges', 'MLP, other edges', 'GP, seed edges', 'GP, other edges'}, 'Location', 'best', 'FontSize', 7);
end
ax = nexttile(tl);  hold(ax, 'on');  box(ax, 'on');  grid(ax, 'on');
plot(ax, dd, D.time.maxN_over_train_max, 'k-', 'LineWidth', 1.5);
yline(ax, pct(2) / Nmax, ':', '99.9th percentile');  yline(ax, 1, '--', 'training max');
set(ax, 'XScale', 'log');  xlim(ax, [0.1 180]);
title(ax, 'largest N at an active edge / largest training N', 'FontWeight', 'normal');  xlabel(ax, 'day');
ax = nexttile(tl);  hold(ax, 'on');  box(ax, 'on');  grid(ax, 'on');
plot(ax, dd, D.time.mlp_w1_sum_ratio, '-', 'Color', col.mlp, 'LineWidth', 1.5);
plot(ax, dd, D.time.gp_w1_sum_ratio, '-', 'Color', col.gp, 'LineWidth', 1.5);
yline(ax, 1, ':');
set(ax, 'XScale', 'log');  xlim(ax, [0.1 180]);
title(ax, '\Sigma W_1 over active edges, surrogate / numeric', 'FontWeight', 'normal');  xlabel(ax, 'day');
legend(ax, {'MLP', 'GP'}, 'Location', 'best');
exportgraphics(fig, fullfile(fdir, 'offset_error_map.png'), 'Resolution', 200);  close(fig);

% seed close-up and mass budget
[~, order] = sort(num.C(seed, :)' + num.C(:, seed), 'descend');
nb = order(order ~= seed);  nb = nb(1:6);
show = day <= 10;
fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [1 1 28 18], 'Color', 'w');
tl = tiledlayout(fig, 3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
for r = [seed; nb]'
    ax = nexttile(tl);  hold(ax, 'on');  box(ax, 'on');
    plot(ax, day(show), num.N(r, show), 'k-', 'LineWidth', 1.6);
    plot(ax, day(show), run.mlp.N(r, show), '--', 'Color', col.mlp, 'LineWidth', 1.3);
    plot(ax, day(show), run.gp.N(r, show), ':', 'Color', col.gp, 'LineWidth', 1.6);
    t = strrep(string(num.region_names{r}), '_', ' ');  if r == seed, t = "seed: " + t; end
    title(ax, t, 'FontWeight', 'normal', 'FontSize', 8);  set(ax, 'FontSize', 7);
end
xlabel(tl, 'day');  ylabel(tl, 'soluble tau N (\muM)');
ax = nexttile(tl, [1 2]);  hold(ax, 'on');  box(ax, 'on');  grid(ax, 'on');
plot(ax, day, M.node_mass_numeric / M.node_mass_numeric(1), 'k-', 'LineWidth', 1.6);
plot(ax, day, M.node_mass_mlp / M.node_mass_numeric(1), '--', 'Color', col.mlp, 'LineWidth', 1.4);
plot(ax, day, M.node_mass_gp / M.node_mass_numeric(1), ':', 'Color', col.gp, 'LineWidth', 1.6);
plot(ax, day, 1 + (em_num(1) - em_num) / M.node_mass_numeric(1), '-', 'Color', [0.7 0.3 0.3], 'LineWidth', 1);
set(ax, 'XScale', 'log');  xlim(ax, [0.4 180]);
legend(ax, {'node mass, numeric', 'node mass, MLP', 'node mass, GP', '1 + edge mass released (numeric)'}, 'Location', 'best', 'FontSize', 7);
title(ax, 'node mass \Sigma V(N + M) / numeric node mass at day 0', 'FontWeight', 'normal');  xlabel(ax, 'day');
exportgraphics(fig, fullfile(fdir, 'offset_seed_mass.png'), 'Resolution', 200);  close(fig);
disp(D.time(ismember(round(D.time.day, 2), [0 0.45 0.9 2.25 4.5 9.9 29.7 150.3]), ...
    {'day', 'maxN_over_train_max', 'mlp_flux_nmse_seed', 'mlp_w1_nmse_seed', 'mlp_w1_nmse_other', 'mlp_w1_sum_ratio', 'gp_w1_nmse_seed', 'gp_w1_sum_ratio'}))
disp(E)
end

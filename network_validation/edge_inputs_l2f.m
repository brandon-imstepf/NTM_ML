function E = edge_inputs_l2f(n_gp, run_set)
% Single-edge inputs met in the network simulations of the 8 l2f truths
% (numeric NTM, 180 days; seeds as in network_checks_l2f), for the reviewer
% responses R2.M2 (interpolation or extrapolation), R2.M4 (edge errors) and
% R1.M1 (GP predictive variance). The inputs of edge (i,j) at time t are
% (gamma1, lambda1, delta, epsilon, N_i(t), N_j(t)) on the edges of the
% connectome (C > 0) with N_i > 0 or N_j > 0, the edges the simulator solves.
%   Coverage: N_i, N_j against the N columns (NRow, NCol) of the MLP training
%     data (data_bias_e5_*.csv, 1.1e6 rows): fraction of inputs with one node
%     empty (N = 0), below the training minimum and its 1st percentile, above the
%     training maximum, its 99th and its 99.9th percentile.
%   Edge error: each surrogate's phi_tau and Delta_tau at the numeric states,
%     against the numeric F and W1 at the same states:
%     sum (s - n)^2 / sum n^2 over edges and time points; also restricted to
%     the low-N edge of the training domain (edges with a nonzero N below the
%     training 1st percentile) and to the other edges.
%   GP variance: predictive standard deviation of the GP flux model at n_gp
%     random network inputs per simulation and at n_gp rows of held-out data
%     (data_bias_e5_2.csv, not in the GP's training subset), and its Spearman
%     correlation with the GP's absolute flux error at the network inputs.
% Writes results/edge_coverage_l2f.csv, results/edge_error_l2f.csv and
% results/gp_variance_l2f.csv; with run_set = "consistent" (l2f_runs.m) the runs
% in outputs_consistent/, into the same files with suffix _consistent.
%
%   E = edge_inputs_l2f();        % n_gp = 2000, run_set "original"
%   E = edge_inputs_l2f(500);
%   E = edge_inputs_l2f(2000, "consistent");

if nargin < 1, n_gp = 2000; end
if nargin < 2, run_set = "original"; end
here = fileparts(mfilename('fullpath'));
cfg = sntm_config();
T = readtable(fullfile(here, 'manifests', 'recovery_truths.csv'), 'TextType', 'string');
T = T(T.convention == "l2f", :);
methods = ["mlp", "sr", "gp", "lr"];
seeds = ["ECl", "CA1", "RSPd"];
[seed_dirs, tag] = l2f_runs(run_set);
seed_file = @(k, m, id) fullfile(here, seed_dirs(k), "rec_l2f__" + m + "__" + seeds(k) + "__" + id + ".mat");

% N columns of the MLP training data
files = dir(fullfile(cfg.ntm_ml_root, 'data', 'data_neural_network_model', 'training_data_e6', 'data_bias_e5_*.csv'));
Ntrain = [];
for f = files'
    D = readmatrix(fullfile(f.folder, f.name));
    Ntrain = [Ntrain; D(:, 5); D(:, 6)]; %#ok<AGROW>
end
lo = min(Ntrain);  hi = max(Ntrain);  pct = prctile(Ntrain, [1 99 99.9]);  p999 = pct(3);
fprintf('training N (%d values): min %.3g, 1st/99th/99.9th pct %.3g/%.3g/%.3g, max %.3g\n', numel(Ntrain), lo, pct, hi);

sur = struct();
for m = methods, sur.(m) = make_surrogate(char(m), cfg); end
G = load(cfg.gp_model_file);
rng(1);
H = readmatrix(fullfile(cfg.ntm_ml_root, 'data', 'data_neural_network_model', 'training_data_e6', 'data_bias_e5_2.csv'));
H = H(randperm(size(H, 1), n_gp), :);
[~, sd_held] = predict(G.stzd_flux_gprModel, H(:, 1:6));
sd_held = sd_held * G.flux_std;

cov_rows = {};  err_rows = {};  gp_rows = {};
for k = 1:numel(seeds)
for s = 1:height(T)
    id = T.truth_id(s);
    if ~all(arrayfun(@(m) isfile(seed_file(k, m, id)), ["numeric", methods])), continue, end
    num = load(seed_file(k, "numeric", id), 'N', 'F', 'W1', 'C', 'meta');
    p = num.meta;
    [i, j] = find(num.C > 0);
    nt = size(num.N, 2);
    Ni = num.N(i, :);  Nj = num.N(j, :);
    active = Ni > 0 | Nj > 0;
    a = Ni(active);  b = Nj(active);
    one_empty = mean(a == 0 | b == 0);
    x = [a(a > 0); b(b > 0)];
    cov_rows(end+1, :) = {seeds(k), id, nnz(active), one_empty, mean(x < lo), mean(x < pct(1)), mean(x > hi), mean(x > pct(2)), mean(x > p999), ...
        max(x) / hi}; %#ok<AGROW>

    % Surrogates at the numeric states, on the active edge-time points
    lin = sub2ind(size(num.C), repmat(i, 1, nt), repmat(j, 1, nt));
    tt = repmat(1:nt, numel(i), 1);
    Fn = reshape(num.F, [], nt);  Wn = reshape(num.W1, [], nt);
    Fn = Fn(sub2ind(size(Fn), lin(active), tt(active)));
    Wn = Wn(sub2ind(size(Wn), lin(active), tt(active)));
    % W1 at the last time point is not computed by the simulator
    w_ok = tt(active) < nt;
    low = (a > 0 & a < pct(1)) | (b > 0 & b < pct(1));
    for m = methods
        Fs = sur.(m).flux_fun(p.gamma1, p.lambda1, p.delta, p.epsilon, a, b);
        Ws = sur.(m).w1_fun(p.gamma1, p.lambda1, p.delta, p.epsilon, a, b, []);
        err_rows(end+1, :) = {seeds(k), id, m, nmse(Fs, Fn), nmse(Ws(w_ok), Wn(w_ok)), mean(low), ...
            nmse(Fs(low), Fn(low)), nmse(Fs(~low), Fn(~low)), ...
            nmse(Ws(w_ok & low), Wn(w_ok & low)), nmse(Ws(w_ok & ~low), Wn(w_ok & ~low))}; %#ok<AGROW>
    end

    % GP predictive standard deviation at random network inputs
    pick = randperm(numel(a), min(n_gp, numel(a)));
    X = [repmat([p.gamma1, p.lambda1, p.delta, p.epsilon], numel(pick), 1), a(pick), b(pick)];
    [mu, sd] = predict(G.stzd_flux_gprModel, X);
    mu = mu * G.flux_std + G.flux_mean;  sd = sd * G.flux_std;
    err = abs(mu - Fn(pick));
    outside = a(pick) > p999 | b(pick) > p999 | a(pick) == 0 | b(pick) == 0;
    gp_rows(end+1, :) = {seeds(k), id, numel(pick), median(sd), median(sd_held), corr(sd, err, 'Type', 'Spearman'), ...
        median(sd(outside)), median(sd(~outside)), median(err(outside)), median(err(~outside))}; %#ok<AGROW>
end
end

E.coverage = cell2table(cov_rows, 'VariableNames', {'seed', 'truth_id', 'n_inputs', 'frac_one_node_empty', ...
    'frac_N_below_train_min', 'frac_N_below_p1', 'frac_N_above_train_max', 'frac_N_above_p99', 'frac_N_above_p999', 'max_N_over_train_max'});
E.error = cell2table(err_rows, 'VariableNames', {'seed', 'truth_id', 'method', 'flux_nmse', 'w1_nmse', ...
    'frac_low_N', 'flux_nmse_low_N', 'flux_nmse_other', 'w1_nmse_low_N', 'w1_nmse_other'});
E.gp = cell2table(gp_rows, 'VariableNames', {'seed', 'truth_id', 'n', 'sd_median_network', 'sd_median_heldout', ...
    'spearman_sd_abs_err', 'sd_median_outside', 'sd_median_inside', 'abs_err_median_outside', 'abs_err_median_inside'});
writetable(E.coverage, fullfile(here, 'results', "edge_coverage_l2f" + tag + ".csv"));
writetable(E.error, fullfile(here, 'results', "edge_error_l2f" + tag + ".csv"));
writetable(E.gp, fullfile(here, 'results', "gp_variance_l2f" + tag + ".csv"));

fprintf('\n== Coverage (fractions of active edge inputs; N fractions over nonzero N)\n');  disp(E.coverage)
fprintf('== Edge error at the numeric states, median [min, max] over simulations\n');
for m = methods
    e = E.error(E.error.method == m, :);
    fprintf('%-4s flux %s   W1 %s\n', m, mmm(e.flux_nmse), mmm(e.w1_nmse));
end
fprintf('\n== GP predictive s.d. (flux units)\n');  disp(E.gp)
end

function v = nmse(s, n)
v = sum((s(:) - n(:)).^2) / sum(n(:).^2);
end

function s = mmm(x)
s = sprintf('%.3g [%.3g, %.3g]', median(x, 'omitnan'), min(x), max(x));
end

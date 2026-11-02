function scan_loss_t60(param, n_points, grid)
% The Nelder-Mead loss of run_nm_t60.slurm on a fixed grid of values over the
% whole training range, one parameter at a time, the others at their true
% values: shows the whole loss curve, so a search that stopped in a local
% minimum can be told from one that found the lowest point. Each MLP
% simulation is scored against both 60-day truths of the grid (numeric and the
% MLP's own run), with the same NMSE as nelder_mead_network_recovery.m.
% Values are evenly spaced in the search coordinates (log10 for gamma1 and
% lambda1, linear for delta and epsilon).
%
%   scan_loss_t60('lambda1')               % 25 values, default grid
%   scan_loss_t60('epsilon', 25, 'dt0.225')
%
% Writes outputs_nm/t60_scan/scan_<grid>__<param>[__n<n_points>].mat (__n only
% when n_points is not 25): truth_id, value (truths x points), nmse_numeric,
% nmse_mlp, runtime_s, and N (truths x points cell of the simulated N, so other
% losses can be scored later: compare_losses_t60.m). Saved after every truth;
% a rerun skips truths already done.
%
%   scan_loss_t60('gamma1', 49, 'dt0.225')   % for compare_losses_t60

if nargin < 2 || isempty(n_points), n_points = 25; end
if nargin < 3 || isempty(grid), grid = "base"; end
param = string(param);  grid = string(grid);
params = ["gamma1", "lambda1", "delta", "epsilon"];
param_range = [5e-5 5e-2; 0.01 0.1; 10 100; 10 100];
p_idx = find(params == param);
assert(~isempty(p_idx), 'param must be among %s', strjoin(params, ', '));
if p_idx <= 2
    values = logspace(log10(param_range(p_idx, 1)), log10(param_range(p_idx, 2)), n_points);
else
    values = linspace(param_range(p_idx, 1), param_range(p_idx, 2), n_points);
end

setenv('SNTM_T_END', '60');
switch grid
    case "base", setenv('SNTM_DT_MAX', '');
    case "dt0.45", setenv('SNTM_DT_MAX', '0.45');
    case "dt0.225", setenv('SNTM_DT_MAX', '0.225');
    otherwise, error('unknown grid %s', grid);
end
here = fileparts(mfilename('fullpath'));
cfg = sntm_config();
addpath(fullfile(cfg.sim_root, 'Code', 'Transport_Network'));
sur = make_surrogate('mlp', cfg);
out_dir = fullfile(here, 'outputs_nm', 't60_scan');
if ~isfolder(out_dir), mkdir(out_dir); end
n_label = "";
if n_points ~= 25, n_label = "__n" + n_points; end
out_file = fullfile(out_dir, "scan_" + grid + "__" + param + n_label + ".mat");

T = readtable(fullfile(here, 'manifests', 'recovery_truths.csv'), 'TextType', 'string');
T = T(T.convention == "l2f", :);
truth_id = T.truth_id;
value = repmat(values, height(T), 1);
nmse_numeric = nan(height(T), n_points);
nmse_mlp = nan(height(T), n_points);
runtime_s = nan(height(T), n_points);
N_sim = cell(height(T), n_points);
if isfile(out_file)
    S = load(out_file);
    assert(isequal(S.value, value), 'saved scan has other values; delete %s', out_file);
    nmse_numeric = S.nmse_numeric;  nmse_mlp = S.nmse_mlp;  runtime_s = S.runtime_s;
    if isfield(S, 'N'), N_sim = S.N; end
end

for s = 1:height(T)
    if all(isfinite(nmse_numeric(s, :))) && all(~cellfun(@isempty, N_sim(s, :))), continue, end
    num = load(fullfile(here, "outputs_t60_" + grid, "rec_l2f__numeric__ECl__" + truth_id(s) + ".mat"), 'N', 'meta');
    mlp = load(fullfile(here, "outputs_t60_" + grid, "rec_l2f__mlp__ECl__" + truth_id(s) + ".mat"), 'N');
    theta_true = [T.gamma1(s), T.lambda1(s), T.delta(s), T.epsilon(s)];
    for k = 1:n_points
        theta = theta_true;  theta(p_idx) = values(k);
        p = struct('gamma1', theta(1), 'lambda1', theta(2), 'lambda2', num.meta.lambda2, ...
            'delta', theta(3), 'epsilon', theta(4));
        t0 = tic;
        out = sntm_simulate(cfg, sur, true, p, num.meta.seed, num.meta.subset, num.meta.init_rescale, k);
        runtime_s(s, k) = toc(t0);
        N = out.Predicted.N;
        N_sim{s, k} = N;
        assert(isequal(size(N), size(num.N)), 'simulation and truth have different time points');
        nmse_numeric(s, k) = sum((N(:) - num.N(:)).^2) / sum(num.N(:).^2);
        nmse_mlp(s, k) = sum((N(:) - mlp.N(:)).^2) / sum(mlp.N(:).^2);
    end
    N = N_sim;
    save(out_file, 'truth_id', 'value', 'nmse_numeric', 'nmse_mlp', 'runtime_s', 'grid', 'param', 'N');
    fprintf('[%s %s] %s done (%.0f s)\n', grid, param, truth_id(s), sum(runtime_s(s, :)));
end
end

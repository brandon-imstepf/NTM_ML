function W = w1_check(id, seed, n_states)
% Delta_tau (W1) of the numeric network run with per-edge flux solves
% (../outputs_fixedflux) against the single-edge code at the same states.
% n_states edge states (i, j, step) are drawn with probability ~ W1^2. At each:
%   W_net    network W1 (MassBalanceCalculator_explicit_W_faster, lambda2 = 0.01)
%   W_se_l1  SingleEdgeWCalculator, lambda2 = lambda1 (the training-data convention)
%   W_se_01  SingleEdgeWCalculator, lambda2 = 0.01 (the network runs' convention)
%   W_mlp    MLP prediction
% both single-edge calls given the network's (per-edge) flux and the network tolerances.
% Writes ../results/w1_check_<seed>_<id>.csv.
%
%   W = w1_check('g0.0005_l0.05_d80_e20', 'ECl', 40);

if nargin < 2, seed = 'ECl'; end
if nargin < 3, n_states = 40; end
here = fileparts(mfilename('fullpath'));
nv = fileparts(here);
addpath(nv);
cfg = sntm_config();
addpath(fullfile(cfg.sim_root, 'Code', 'Transport_Network'));
r = load(fullfile(nv, 'outputs_fixedflux', ['rec_l2f__numeric__' seed '__' id '.mat']), 'N', 'F', 'W1', 'C', 'meta');
p = r.meta;
[i, j] = find(r.C > 0);
nt = size(r.N, 2);
cand = zeros(0, 3);
for h = 1:nt-1
    on = r.N(i, h) > 0 | r.N(j, h) > 0;
    cand = [cand; i(on), j(on), repmat(h, nnz(on), 1)]; %#ok<AGROW>
end
w = r.W1(sub2ind(size(r.W1), cand(:, 1), cand(:, 2), cand(:, 3))).^2;
rng(3);
pick = arrayfun(@(~) find(cumsum(w) / sum(w) >= rand, 1), 1:n_states)';
fixed = {'beta', cfg.beta, 'gamma2', cfg.gamma2, 'frac', cfg.frac, 'L_int', cfg.L_int, 'L1', cfg.L1, 'L2', cfg.L2, ...
    'L_ais', cfg.L_ais, 'L_syn', cfg.L_syn, 'resmesh', cfg.resmesh, 'reltol', cfg.reltol, 'abstol', cfg.abstol, ...
    'fsolvetol', cfg.fsolvetol, 'gamma1', p.gamma1, 'lambda1', p.lambda1, 'delta', p.delta, 'epsilon', p.epsilon};
M = load(cfg.mlp_weights_file);
rows = zeros(n_states, 7);
for k = 1:n_states
    c = cand(pick(k), :);  n1 = r.N(c(1), c(3));  n2 = r.N(c(2), c(3));  F = r.F(c(1), c(2), c(3));
    Wl1 = SingleEdgeWCalculator(cfg.matdir, F, fixed{:}, 'lambda2', p.lambda1, 'tau_x0', n1, 'tau_xL', n2, 'use_sr_w1', 0);
    W01 = SingleEdgeWCalculator(cfg.matdir, F, fixed{:}, 'lambda2', 0.01, 'tau_x0', n1, 'tau_xL', n2, 'use_sr_w1', 0);
    rows(k, :) = [c, r.W1(c(1), c(2), c(3)), Wl1, W01, mlp_forward(M.w1, [p.gamma1, p.lambda1, p.delta, p.epsilon, n1, n2])];
end
W = array2table(rows, 'VariableNames', {'i', 'j', 'step', 'W_net', 'W_se_l1', 'W_se_01', 'W_mlp'});
writetable(W, fullfile(nv, 'results', sprintf('w1_check_%s_%s.csv', seed, id)));
d = @(s, n) sum((s - n).^2) / sum(n.^2);
fprintf('%s %s, %d states weighted by W1^2 (NMSE):\n', seed, id, n_states);
fprintf('  network vs single-edge (l2 = l1): %.3g\n  network vs single-edge (l2 = 0.01): %.3g\n', ...
    d(W.W_net, W.W_se_l1), d(W.W_net, W.W_se_01));
fprintf('  single-edge l2 = 0.01 vs l2 = l1: %.3g\n  MLP vs single-edge (l2 = l1): %.3g\n  MLP vs network: %.3g\n', ...
    d(W.W_se_01, W.W_se_l1), d(W.W_mlp, W.W_se_l1), d(W.W_mlp, W.W_net));
end

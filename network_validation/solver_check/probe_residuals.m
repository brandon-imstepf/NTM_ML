% Residuals of the two candidate fluxes of edge 4 -> 15 (step 57), solved alone,
% with and without ode45's non-negativity constraint, and fsolve's exit flag in
% the original joint solve. Needs probe_target.m's saved .mat; set NTM_SIM_ROOT
% to the jtorok_dev_orig copy.
here = fileparts(mfilename('fullpath'));
addpath(here);
nv = fileparts(here);
addpath(nv);
cfg = sntm_config();
addpath(fullfile(cfg.sim_root, 'Code', 'Transport_Network'));
id = 'g0.0005_l0.02_d80_e20';  h = 57;  tgt = 15;
S = load(fullfile(here, sprintf('probe_%s_h%d_t%d.mat', id, h, tgt)));
r = load(fullfile(nv, 'outputs_recovery', ['rec_l2f__numeric__ECl__' id '.mat']), 'N', 'C', 'meta');
p = r.meta;  N = r.N(:, h);
% a 1-edge network: only source 4 and target 15 nonzero, so the target's stack holds
% the edges from 4 and from any other source; restrict to edge 4 via a masked Adj
common = {'beta', cfg.beta, 'gamma1', p.gamma1, 'gamma2', cfg.gamma2, 'delta', p.delta, 'epsilon', p.epsilon, ...
    'lambda1', p.lambda1, 'lambda2', p.lambda2, 'frac', cfg.frac, 'L_int', cfg.L_int, 'L1', cfg.L1, 'L2', cfg.L2, ...
    'L_ais', cfg.L_ais, 'L_syn', cfg.L_syn, 'resmesh', cfg.resmesh, 'connectome_subset', 'Hippocampus+PC+RSP', ...
    'reltol', cfg.reltol, 'abstol', cfg.abstol, 'fsolvetol', cfg.fsolvetol, 'only_target', tgt};
Adj = r.C > 0;
X0 = N .* Adj;  Nt = N;     % edge 4 is first in region 15's stack; one value = that edge alone
cands = [S.Aa(1), S.Ab(1)];
fprintf('edge 4 -> 15, N_src %.3g, N_tgt %.3g\n', N(4), N(tgt));
for nn = [0 1]
    for A = cands
        [~, ~, e] = nfc_probe(X0, Nt, cfg.matdir, common{:}, 'nonneg_n', nn, 'probe_A', A);
        fprintf('  A = %.4g, nonneg %d: residual %.3g, min n %.3g\n', A, nn, e.res, e.min_n);
    end
end
% fsolve exit flag and final residual of the original joint solve of region 15
[Fj, ~, ij] = nfc_probe(N .* Adj, N, cfg.matdir, common{:});
fprintf('joint solve (original) flux of edge 4: %.4g; fsolve exitflag %d after %d iterations; final residuals %s\n', ...
    Fj(4, tgt), ij.exitflag, ij.iterations, mat2str(ij.fval', 3));

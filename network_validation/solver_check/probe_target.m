% Solver-disagreement probe on the edges into one region at one stored state.
% Uses nfc_probe.m (patched copy of the ORIGINAL NetworkFluxCalculator.m).
% Make nfc_probe.m first (from the ORIGINAL simulator code):
%   python -I make_probe.py <jtorok_dev_orig>/Code/Transport_Network/NetworkFluxCalculator.m nfc_probe.m
% Set NTM_SIM_ROOT to the jtorok_dev_orig copy so the single-edge code is the original too.
here = fileparts(mfilename('fullpath'));
addpath(here);
nv = fileparts(here);
addpath(nv);
cfg = sntm_config();
addpath(fullfile(cfg.sim_root, 'Code', 'Transport_Network'));
id = getenv_or('PROBE_TRUTH', 'g0.0005_l0.02_d80_e20');
h = str2double(getenv_or('PROBE_H', '57'));
tgt = str2double(getenv_or('PROBE_TARGET', '15'));
r = load(fullfile(nv, 'outputs_recovery', ['rec_l2f__numeric__ECl__' id '.mat']), 'N', 'F', 'C', 'meta');
p = r.meta;  Adj = r.C > 0;  N = r.N(:, h);
common = {'beta', cfg.beta, 'gamma1', p.gamma1, 'gamma2', cfg.gamma2, 'delta', p.delta, 'epsilon', p.epsilon, ...
    'lambda1', p.lambda1, 'lambda2', p.lambda2, 'frac', cfg.frac, 'L_int', cfg.L_int, 'L1', cfg.L1, 'L2', cfg.L2, ...
    'L_ais', cfg.L_ais, 'L_syn', cfg.L_syn, 'resmesh', cfg.resmesh, 'connectome_subset', 'Hippocampus+PC+RSP', ...
    'reltol', cfg.reltol, 'abstol', cfg.abstol, 'fsolvetol', cfg.fsolvetol, 'only_target', tgt};
run = @(varargin) nfc_probe(N .* Adj, N, cfg.matdir, common{:}, varargin{:});

[Fa, ~, ia] = run('probe_mode', 'joint', 'nonneg_n', 1);          % original behaviour
nE = numel(ia.rows);
[Fb, ~] = run('probe_mode', 'per_edge', 'nonneg_n', 1);          % each edge alone
[Fc, ~] = run('probe_mode', 'joint', 'nonneg_n', nE);            % joint, all components non-negative
Aa = Fa(ia.rows, tgt);  Ab = Fb(ia.rows, tgt);  Ac = Fc(ia.rows, tgt);
% residual and smallest n along the edge, without any non-negativity constraint
[~, ~, ea] = run('probe_mode', 'joint', 'nonneg_n', 0, 'probe_A', Aa);
[~, ~, eb] = run('probe_mode', 'joint', 'nonneg_n', 0, 'probe_A', Ab);
% the single-edge code at the same (N_i, N_j) and the same tolerances
Fs = zeros(nE, 1);
fixed = {'beta', cfg.beta, 'gamma2', cfg.gamma2, 'frac', cfg.frac, 'L_int', cfg.L_int, 'L1', cfg.L1, 'L2', cfg.L2, ...
    'L_ais', cfg.L_ais, 'L_syn', cfg.L_syn, 'resmesh', cfg.resmesh, 'reltol', cfg.reltol, 'abstol', cfg.abstol, ...
    'fsolvetol', cfg.fsolvetol, 'gamma1', p.gamma1, 'lambda1', p.lambda1, 'lambda2', p.lambda2, ...
    'delta', p.delta, 'epsilon', p.epsilon};
for e = 1:nE
    o = SingleEdgeFluxCalculator_WrapperFun(cfg.matdir, fixed{:}, 'tau_x0', ia.B(e), 'tau_xL', ia.C(e), ...
        'sim_no', e, 'use_sr_flux', 0, 'use_sr_w1', 0);
    Fs(e) = o.Predicted.F;
end
stored = r.F(ia.rows, tgt, h);

fprintf('\n%s, step %d, edges into region %d (%d active; stack order = row order)\n', id, h, tgt, nE);
fprintf('%4s %4s %10s %10s | %10s %10s %10s %10s | %10s %10s | %10s %10s\n', 'pos', 'src', 'N_src', 'N_tgt', ...
    'stored', 'joint', 'per_edge', 'joint_nn', 'single', 'resid(jt)', 'minN(jt)', 'minN(pe)');
for e = 1:nE
    fprintf('%4d %4d %10.3g %10.3g | %10.3g %10.3g %10.3g %10.3g | %10.3g %10.2g | %10.3g %10.3g\n', e, ia.rows(e), ...
        ia.B(e), ia.C(e), stored(e), Aa(e), Ab(e), Ac(e), Fs(e), ea.res(e), ea.min_n(e), eb.min_n(e));
end
save(fullfile(here, sprintf('probe_%s_h%d_t%d.mat', id, h, tgt)), 'Aa', 'Ab', 'Ac', 'Fs', 'stored', 'ia', 'ea', 'eb');

function v = getenv_or(name, default)
v = getenv(name); if isempty(v), v = default; end
end

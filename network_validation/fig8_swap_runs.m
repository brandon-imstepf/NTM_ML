function fig8_swap_runs(which)
% Swap test for the offset of the new Fig. 8 runs (FIG8_OFFSET.md): at the Fig. 8
% parameters (gamma1 5e-5, lambda 0.01 = lambda2, delta 10, epsilon 50, ECl seed,
% init_rescale 0.01, 180 days, default grid, joint flux solve), take phi_tau and
% Delta_tau (W1) from different sources:
%   "A"  numeric phi_tau, MLP Delta_tau     -> outputs_fig8/fig8__swapA_numflux_mlpw1__...
%   "B"  MLP phi_tau, numeric Delta_tau     -> outputs_fig8/fig8__swapB_mlpflux_numw1__...
% In B the numeric W1 is computed by the simulator from the edge profile of the
% MLP flux (W1 is a function of the flux and the profile).
% Saved in the format of run_sntm_cases.m; a finished run is skipped.
% Compare with fig8__numeric__ and fig8__mlp__ (offset: sum N MLP / numeric = 0.89).
%
%   fig8_swap_runs("A");  fig8_swap_runs("B");   % one MATLAB process each

here = fileparts(mfilename('fullpath'));
cfg = sntm_config();
addpath(fullfile(cfg.sim_root, 'Code', 'Transport_Network'));
which = string(which);
switch which
    case "A", use = [0 1];  tag = "swapA_numflux_mlpw1";
    case "B", use = [1 0];  tag = "swapB_mlpflux_numw1";
    otherwise, error('which is "A" or "B"');
end
out_file = fullfile(here, 'outputs_fig8', "fig8__" + tag + "__ECl__g5e-05_l0.01_d10_e50.mat");
if isfile(out_file), fprintf('%s exists, skipping\n', out_file); return; end
maxNumCompThreads(1);
sur = make_surrogate('mlp', cfg);
p = struct('gamma1', 5e-5, 'lambda1', 0.01, 'lambda2', 0.01, 'delta', 10, 'epsilon', 50);
init_rescale = 0.01;
t0 = tic;
out = sntm_simulate(cfg, sur, use, p, "ECl", "Hippocampus+PC+RSP", init_rescale, 1);
res = struct();
res.meta = struct('case_id', "fig8__" + tag + "__ECl__g5e-05_l0.01_d10_e50", 'method', tag, ...
    'use_flux_surrogate', use(1), 'use_w1_surrogate', use(2), 'surrogate', "mlp", 'seed', "ECl", ...
    'subset', "Hippocampus+PC+RSP", 'gamma1', p.gamma1, 'lambda1', p.lambda1, 'delta', p.delta, ...
    'epsilon', p.epsilon, 'init_rescale', init_rescale, 'lambda2', p.lambda2);
res.N = out.Predicted.N;  res.M = out.Predicted.M;  res.F = out.Predicted.F;  res.W1 = out.Predicted.W1;
if isfield(out.Predicted, 'EdgeMass'), res.EdgeMass = out.Predicted.EdgeMass; end
res.trange = out.Sim.trange;  res.region_names = out.Sim.region_names;  res.C = out.Sim.C;
res.runtime_s = toc(t0);
res.cfg = rmfield(cfg, 'seeds');  res.cfg.init_rescale = init_rescale;  res.cfg.lambda2 = p.lambda2;
save(out_file, '-struct', 'res', '-v7.3');
fprintf('[%s] done in %.0f s\n', tag, res.runtime_s);
end

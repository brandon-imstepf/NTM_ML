function S = flux_solve_survey(every, seeds)
% How often the numeric network flux solve fails to converge, and what it does
% to the flux. For the numeric runs of the 8 l2f truths (180 days), at every
% 'every'-th time step, the stored state N is re-solved with nfc_fast_probe.m
% (instrumented copy of the edited NetworkFluxCalculator.m; same flux as the
% original): once jointly per target region, as the simulator does, and once
% with one scalar fsolve per edge. Per state:
%   n_edges         edges solved
%   n_unconv_joint  edges whose region's joint fsolve ended with exitflag <= 0
%   n_unconv_edge   edges whose own fsolve ended with exitflag <= 0
%   max_res_joint   largest |shooting residual| of the joint solve / max N
%   max_res_edge    the same for the per-edge solves
%   nmse_joint_edge sum (F_joint - F_edge)^2 / sum F_edge^2 over the edges
%   stored_vs_joint NMSE of the stored flux against the joint re-solve (check: ~0)
% Writes ../results/flux_solve_survey.csv. Make nfc_fast_probe.m first:
%   python -I make_fast_probe.py <jtorok_dev>/Code/Transport_Network/NetworkFluxCalculator.m nfc_fast_probe.m
%
%   S = flux_solve_survey();              % every = 5, seed ECl
%   S = flux_solve_survey(1, ["ECl", "CA1", "RSPd"]);

if nargin < 1, every = 5; end
if nargin < 2, seeds = "ECl"; end
here = fileparts(mfilename('fullpath'));
nv = fileparts(here);
addpath(here, nv);
cfg = sntm_config();
addpath(fullfile(cfg.sim_root, 'Code', 'Transport_Network'));
T = readtable(fullfile(nv, 'manifests', 'recovery_truths.csv'), 'TextType', 'string');
T = T(T.convention == "l2f", :);
dir_of = containers.Map({'ECl', 'CA1', 'RSPd'}, {'outputs_recovery', 'outputs_seeds', 'outputs_seeds'});
nmse = @(s, n) sum((s(:) - n(:)).^2) / sum(n(:).^2);

rows = {};
for sd = seeds
for s = 1:height(T)
    id = T.truth_id(s);
    r = load(fullfile(nv, dir_of(char(sd)), "rec_l2f__numeric__" + sd + "__" + id + ".mat"), 'N', 'F', 'C', 'meta', 'trange');
    p = r.meta;  Adj = r.C > 0;
    common = {'beta', cfg.beta, 'gamma1', p.gamma1, 'gamma2', cfg.gamma2, 'delta', p.delta, 'epsilon', p.epsilon, ...
        'lambda1', p.lambda1, 'lambda2', p.lambda2, 'frac', cfg.frac, 'L_int', cfg.L_int, 'L1', cfg.L1, 'L2', cfg.L2, ...
        'L_ais', cfg.L_ais, 'L_syn', cfg.L_syn, 'resmesh', cfg.resmesh, 'connectome_subset', char(p.subset), ...
        'reltol', cfg.reltol, 'abstol', cfg.abstol, 'fsolvetol', cfg.fsolvetol};
    for h = 1:every:size(r.N, 2)
        N = r.N(:, h);
        [Fj, ~, ij] = nfc_fast_probe(N .* Adj, N, cfg.matdir, common{:});
        [Fe, ~, ie] = nfc_fast_probe(N .* Adj, N, cfg.matdir, common{:}, 'per_edge', 1);
        on = ~isnan(ij.flag);
        rows(end+1, :) = {sd, id, h, r.trange(h) * 180, nnz(on), nnz(ij.flag(on) <= 0), nnz(ie.flag(on) <= 0), ...
            max(abs(ij.res(on))) / max(N), max(abs(ie.res(on))) / max(N), nmse(Fj(on), Fe(on)), ...
            nmse(r.F(:, :, h), Fj)}; %#ok<AGROW>
        fprintf('%s %s step %3d: %3d edges, unconverged joint %3d, per-edge %3d, NMSE joint vs per-edge %.3g\n', ...
            sd, id, h, rows{end, 5}, rows{end, 6}, rows{end, 7}, rows{end, 10});
    end
end
end
S = cell2table(rows, 'VariableNames', {'seed', 'truth_id', 'step', 'day', 'n_edges', 'n_unconv_joint', ...
    'n_unconv_edge', 'max_res_joint', 'max_res_edge', 'nmse_joint_edge', 'stored_vs_joint'});
writetable(S, fullfile(nv, 'results', 'flux_solve_survey.csv'));
end

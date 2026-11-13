function D = time_step_consistent()
% Time-step sensitivity of the MLP sNTM in the surrogates' training convention
% (l2f_runs.m, "consistent"), for R2.M4. The 8 l2f settings, seed ECl, are run
% for 60 days with W1 and edge mass at lambda2 = lambda1 (SNTM_LAMBDA2_W1), on
% the default grid and with the time step limited to 0.225 days (SNTM_DT_MAX),
% into outputs_t60c_base/ and outputs_t60c_dt0.225/ (finished runs are skipped).
% Per setting, at the shared time points up to day 59:
%   mlp_base_vs_dt0225  NMSE (Eq. 7, nmse_eq7.m) of the MLP sNTM on the default grid
%                       against the MLP sNTM on the 0.225-day grid
%   mlp_vs_numeric      NMSE of the MLP sNTM against the numeric NTM on the
%                       default grid (the 180-day runs in outputs_consistent/)
% Writes results/network_dt_mlp_l2f_consistent.csv.
%
%   D = time_step_consistent();

here = fileparts(mfilename('fullpath'));
manifest = fullfile(here, 'manifests', 'consistent_l2f_surrogate.csv');
M = readtable(manifest, 'TextType', 'string');
rows = find(M.method == "mlp" & M.seed == "ECl");
grids = ["base", "dt0.225"];
dt_max = ["", "0.225"];

setenv('SNTM_T_END', '60');
setenv('SNTM_LAMBDA2_W1', 'lambda1');
for g = 1:2
    setenv('SNTM_DT_MAX', dt_max(g));
    setenv('SNTM_OUT_DIR', fullfile(here, "outputs_t60c_" + grids(g)));
    for r = rows'
        run_sntm_cases(manifest, r, r);
    end
end
setenv('SNTM_T_END', '');  setenv('SNTM_LAMBDA2_W1', '');  setenv('SNTM_DT_MAX', '');  setenv('SNTM_OUT_DIR', '');

nmse = @nmse_eq7;   % Eq. 7 (sse_rel.m before 2026-11-10)
out = {};
for r = rows'
    id = M.case_id(r);
    b = load(fullfile(here, 'outputs_t60c_base', id + ".mat"), 'N', 'trange');
    f = load(fullfile(here, 'outputs_t60c_dt0.225', id + ".mat"), 'N', 'trange');
    n = load(fullfile(here, 'outputs_consistent', strrep(id, "__mlp__", "__numeric__") + ".mat"), 'N', 'trange');
    [~, ib, jf] = intersect(round(b.trange * 180, 6), round(f.trange * 180, 6));
    keep = b.trange(ib) * 180 <= 59;
    ib = ib(keep);  jf = jf(keep);
    % the 60-day default grid is the first part of the 180-day grid
    assert(all(abs(n.trange(ib) - b.trange(ib)) < 1e-12));
    out(end+1, :) = {M.setting_id(r), nmse(b.N(:, ib), f.N(:, jf)), nmse(b.N(:, ib), n.N(:, ib))}; %#ok<AGROW>
end
D = cell2table(out, 'VariableNames', {'setting_id', 'mlp_base_vs_dt0225', 'mlp_vs_numeric'});
writetable(D, fullfile(here, 'results', 'network_dt_mlp_l2f_consistent.csv'));
disp(D)
end

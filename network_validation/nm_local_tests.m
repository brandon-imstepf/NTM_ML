function nm_local_tests(test_dir)
% Nelder-Mead network recovery tests on one truth (g0.0005_l0.02_d80_e20, l2f),
% run on the PC on 10/29/26, and their figures:
%   test A  MLP fitted to the MLP's own run at the true values: the loss is 0
%           at the truth, so this tests the optimizer alone
%   test B  MLP fitted to the numeric run: tests recovery with surrogate error
%           (no lambda1 search: test A showed lambda1 needs a better start)
%   lambda1 at fixed values, to draw the loss along lambda1 in both figures
%
% test_dir holds truth/ (the MLP and numeric runs at the true values, from
% run_sntm_cases rows recovery_single_l2f.csv 6 and recovery_truth_l2f.csv 1),
% and gets testA/, testB/, l1_fixed/ and figures/. Every part skips what is
% already saved (nelder_mead_network_recovery.m replays finished searches
% without simulating), so on the saved results this only redraws the figures.
% From scratch, the searches take about 1-1.5 h each and the lambda1 runs
% about 45 min, one after the other.
%
%   nm_local_tests()   % uses outputs_nm/local_tests_10-29

here = fileparts(mfilename('fullpath'));
if nargin < 1, test_dir = fullfile(here, 'outputs_nm', 'local_tests_10-29'); end
truth_id = 'g0.0005_l0.02_d80_e20';
truth_mlp = fullfile(test_dir, 'truth', 'rec_l2f__mlp__ECl__g0.0005_l0.02_d80_e20.mat');
truth_numeric = fullfile(test_dir, 'truth', 'rec_l2f__numeric__ECl__g0.0005_l0.02_d80_e20.mat');

%% 1. Searches (nelder_mead_network_recovery.m reads the truth and output folder from these)
runSearches(truth_id, ["gamma1", "lambda1", "delta", "epsilon"], truth_mlp, fullfile(test_dir, 'testA'));
runSearches(truth_id, ["gamma1", "delta", "epsilon"], truth_numeric, fullfile(test_dir, 'testB'));

%% 2. MLP at fixed lambda1 values, the other parameters at the truth
l1_dir = fullfile(test_dir, 'l1_fixed');
if ~isfolder(l1_dir), mkdir(l1_dir); end
l1_values = [0.01 0.0125 0.015 0.0175 0.0225 0.025 0.045 0.06 0.08 0.1];
cfg = sntm_config();
addpath(fullfile(cfg.sim_root, 'Code', 'Transport_Network'));
truth = load(truth_mlp, 'meta');
sur = [];
for l1 = l1_values
    f = fullfile(l1_dir, sprintf('l1_%g.mat', l1));
    if isfile(f), continue, end
    if isempty(sur), sur = make_surrogate('mlp', cfg); end
    p = struct('gamma1', 5e-4, 'lambda1', l1, 'lambda2', truth.meta.lambda2, 'delta', 80, 'epsilon', 20);
    out = sntm_simulate(cfg, sur, true, p, truth.meta.seed, truth.meta.subset, truth.meta.init_rescale, 1);
    N = out.Predicted.N;
    save(f, 'l1', 'N');
    fprintf('lambda1 = %g done\n', l1);
end

%% 3. Figures
extra = struct('param', {}, 'value', {}, 'N', {});
for l1 = l1_values
    S = load(fullfile(l1_dir, sprintf('l1_%g.mat', l1)));
    extra(end+1) = struct('param', 'lambda1', 'value', S.l1, 'N', S.N); %#ok<AGROW>
end
fig_dir = fullfile(test_dir, 'figures');
plotNelderMeadRecovery(fullfile(test_dir, 'testA'), fig_dir, 'testA_mlp_truth', extra, truth_mlp);
% against the numeric truth, every MLP run with only lambda1 changed is a point
% on the curve: the fixed values, test A's lambda1 search, and the truth itself
A = load(fullfile(test_dir, 'testA', sprintf('nm_l2f__mlp__%s__lambda1.mat', truth_id)), 'trace');
for t = A.trace
    extra(end+1) = struct('param', 'lambda1', 'value', t.value, 'N', t.N); %#ok<AGROW>
end
M = load(truth_mlp, 'N');
extra(end+1) = struct('param', 'lambda1', 'value', 0.02, 'N', M.N);
plotNelderMeadRecovery(fullfile(test_dir, 'testB'), fig_dir, 'testB_numeric_truth', extra, truth_numeric);
end

function runSearches(truth_id, params, truth_file, out_dir)
    setenv('NM_TRUTH_FILE', truth_file);
    setenv('NM_OUT_DIR', out_dir);
    for p = params
        nelder_mead_network_recovery(truth_id, p);
    end
    setenv('NM_TRUTH_FILE', '');
    setenv('NM_OUT_DIR', '');
end

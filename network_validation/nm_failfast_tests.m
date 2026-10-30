function nm_failfast_tests(parts, test_dir, t_end, start_frac)
% Quick ("fail fast") Nelder-Mead checks on the truth of nm_local_tests.m
% (g0.0005_l0.02_d80_e20, l2f): searches started close to the truth, with small
% first steps, and shortened simulations. The parameters not searched are fixed
% at the truth.
%   A1  gamma1 alone,              MLP fitted to the MLP's own run (loss 0 at the truth)
%   B1  gamma1 alone,              MLP fitted to the numeric run
%   A2  gamma1 and delta together, MLP fitted to the MLP's own run
%   B2  gamma1 and delta together, MLP fitted to the numeric run
% A tests the search alone; in B the minimum of the loss is where the MLP fits
% the numeric run best, which need not be the truth.
%
% The truths are made from the 180-day runs in outputs_nm/local_tests_10-29/truth:
% the numeric run (11 h) is cut to the time points up to t_end, the MLP run is
% simulated again up to t_end (so its loss at the truth is exactly 0). Every part
% skips what is already saved; on the saved results this only redraws the figures.
%
%   nm_failfast_tests()               % all four, outputs_nm/failfast_50d, 50 days, start 0.98 x truth
%   nm_failfast_tests(["A1", "A2"])   % some of them (e.g. in two MATLAB sessions at once)
%
% parts       any of "A1", "B1", "A2", "B2" (default all)
% test_dir    output folder (gets truth/, testA/, testB/, figures/)
% t_end       days simulated (SNTM_T_END)
% start_frac  start of the searches as a fraction of the true values
%
% First steps: 0.005 decades in gamma1 (1.2%), 1 in delta (1.25%).

here = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(parts), parts = ["A1", "B1", "A2", "B2"]; end
if nargin < 3 || isempty(t_end), t_end = 50; end
if nargin < 2 || isempty(test_dir), test_dir = fullfile(here, 'outputs_nm', sprintf('failfast_%gd', t_end)); end
if nargin < 4 || isempty(start_frac), start_frac = 0.98; end
parts = string(parts);
truth_id = 'g0.0005_l0.02_d80_e20';
case_name = 'ECl__g0.0005_l0.02_d80_e20.mat';
full_dir = fullfile(here, 'outputs_nm', 'local_tests_10-29', 'truth');
truth_dir = fullfile(test_dir, 'truth');
if ~isfolder(truth_dir), mkdir(truth_dir); end
searches = struct( ...
    'part',  {"A1", "B1", "A2", "B2"}, ...
    'truth', {"mlp", "numeric", "mlp", "numeric"}, ...
    'param', {"gamma1", "gamma1", ["gamma1", "delta"], ["gamma1", "delta"]}, ...
    'true',  {5e-4, 5e-4, [5e-4, 80], [5e-4, 80]}, ...
    'step',  {0.005, 0.005, [0.005, 1], [0.005, 1]});
searches = searches(ismember([searches.part], parts));

setenv('SNTM_T_END', num2str(t_end));
cfg = sntm_config();
addpath(fullfile(cfg.sim_root, 'Code', 'Transport_Network'));

%% 1. Truths up to t_end
truth_file = struct('numeric', fullfile(truth_dir, ['rec_l2f__numeric__' case_name]), ...
                    'mlp', fullfile(truth_dir, ['rec_l2f__mlp__' case_name]));
if any([searches.truth] == "numeric") && ~isfile(truth_file.numeric)
    S = load(fullfile(full_dir, ['rec_l2f__numeric__' case_name]), 'N', 'trange', 'meta', 'region_names');
    keep = ismember(S.trange, cfg.trange);
    assert(isequal(S.trange(keep), cfg.trange), 'the full run does not have the time points of SNTM_T_END');
    S.N = S.N(:, keep);  S.trange = S.trange(keep);
    save(truth_file.numeric, '-struct', 'S');
end
if any([searches.truth] == "mlp") && ~isfile(truth_file.mlp)
    S = load(fullfile(full_dir, ['rec_l2f__mlp__' case_name]), 'meta', 'region_names');
    m = S.meta;
    p = struct('gamma1', m.gamma1, 'lambda1', m.lambda1, 'lambda2', m.lambda2, 'delta', m.delta, 'epsilon', m.epsilon);
    t0 = tic;
    out = sntm_simulate(cfg, make_surrogate('mlp', cfg), true, p, m.seed, m.subset, m.init_rescale, 1);
    S.runtime_s = toc(t0);
    S.N = out.Predicted.N;  S.trange = cfg.trange;
    save(truth_file.mlp, '-struct', 'S');
    fprintf('MLP truth up to %g days: %.1f min\n', t_end, S.runtime_s / 60);
end

%% 2. Searches from start_frac x truth
for s = searches
    test = "testA";
    if s.truth == "numeric", test = "testB"; end
    setenv('NM_TRUTH_FILE', truth_file.(s.truth));
    setenv('NM_OUT_DIR', fullfile(test_dir, test));
    setenv('NM_START', strjoin(compose('%.10g', start_frac * s.true), ','));
    setenv('NM_STEP', strjoin(compose('%g', s.step), ','));
    nelder_mead_network_recovery(truth_id, s.param);
end
for v = ["NM_TRUTH_FILE", "NM_OUT_DIR", "NM_START", "NM_STEP", "SNTM_T_END"], setenv(v, ''); end

%% 3. Figures (of every finished search in test_dir, not only parts)
fig_dir = fullfile(test_dir, 'figures');
for test = ["testA", "testB"]
    tf = truth_file.mlp;  name = "mlp_truth";
    if test == "testB", tf = truth_file.numeric;  name = "numeric_truth"; end
    d = fullfile(test_dir, test);
    one = fullfile(d, sprintf('nm_l2f__mlp__%s__gamma1.mat', truth_id));
    if isfile(one) && hasResult(one)
        plotNelderMeadRecovery(d, fig_dir, "failfast_" + test + "_" + name, [], tf);
    end
    two = fullfile(d, sprintf('nm_l2f__mlp__%s__gamma1+delta.mat', truth_id));
    if isfile(two) && hasResult(two)
        plotNelderMeadRecovery2D(two, fig_dir, "failfast_" + test + "_" + name + "_gamma1_delta", tf);
    end
end
end

function tf = hasResult(f)
    tf = ~isempty(whos('-file', f, 'result'));
end

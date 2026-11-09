function nm_consistent(worker, n_workers, kinds)
% Network parameter recovery with the MLP sNTM against the numeric truths of
% the consistent runs (outputs_consistent/: per-edge flux solve, W1 at
% lambda2 = lambda1, the surrogates' training convention; see
% solver_check/SOLVER_CHECK.md, sections 11-12). 180 days, default time grid.
% Replaces, for the recovery document, the searches of run_nm_t60.slurm,
% whose numeric truths used the joint flux solve (unconverged at gamma1 5e-4,
% delta 80, epsilon 20) and W1 at lambda2 = 0.01.
%
% Task kinds (all for the 8 l2f truths of manifests/recovery_truths.csv):
%   "scan"     loss along each parameter: MLP sNTM at 25 values evenly spaced
%              over the training range (log10 for gamma1 and lambda1), the
%              other three at their true values; ECl seed. Scored against the
%              numeric truth and against the MLP's own run at the true values.
%   "single"   Nelder-Mead (nelder_mead_network_recovery.m), one parameter at a
%              time, the others at their true values; start in the middle of
%              the range; numeric truth; seeds ECl, CA1, RSPd.
%   "control"  as "single", ECl, against the MLP's own run (loss 0 at the truth:
%              tests the search alone).
%   "joint"    Nelder-Mead over all four parameters together, start in the
%              middle of the range, numeric truth, ECl, at most 160 simulations.
%   "twostage" as "single", ECl, numeric truth, but started from the lowest
%              point of the 25-point loss scan ("scan", which must be done),
%              with a first step of half the scan spacing; the scan is the
%              first stage of the search, Nelder-Mead the second.
%              Output: ECl/numeric_truth_twostage/.
% Tasks are numbered in the order scan, single, control, twostage, joint (within a kind:
% truth, then parameter, then seed); worker w of n runs tasks w, w + n, ...
% Finished tasks are skipped, so a rerun resumes.
%
%   nm_consistent(1, 1)                       % everything, one process
%   nm_consistent(2, 6, ["single", "joint"])  % worker 2 of 6
%
% Outputs in outputs_nm/consistent/ (not in git):
%   scan/scan__<truth_id>__<param>.mat         value, nmse_numeric, nmse_mlp, N, runtime_s
%   <seed>/numeric_truth/nm_l2f__mlp__<truth_id>__<param>.mat   (single, joint)
%   ECl/mlp_truth/nm_l2f__mlp__<truth_id>__<param>.mat          (control)
% Tables and figures: analyze_nm_consistent.m.

if nargin < 1 || isempty(worker), worker = 1; end
if nargin < 2 || isempty(n_workers), n_workers = 1; end
if nargin < 3 || isempty(kinds), kinds = ["scan", "single", "control", "joint", "twostage"]; end
kinds = string(kinds);
maxNumCompThreads(1);

here = fileparts(mfilename('fullpath'));
setenv('SNTM_LAMBDA2_W1', 'lambda1');
setenv('SNTM_T_END', '');
setenv('SNTM_DT_MAX', '');
setenv('SNTM_FLUX_PER_EDGE', '');   % surrogate runs: the flux comes from the MLP
out_root = fullfile(here, 'outputs_nm', 'consistent');
T = readtable(fullfile(here, 'manifests', 'recovery_truths.csv'), 'TextType', 'string');
T = T(T.convention == "l2f", :);
params = ["gamma1", "lambda1", "delta", "epsilon"];
seeds = ["ECl", "CA1", "RSPd"];

tasks = struct('kind', {}, 'truth_id', {}, 'param', {}, 'seed', {});
for s = 1:height(T)
    for p = params, tasks(end+1) = struct('kind', "scan", 'truth_id', T.truth_id(s), 'param', p, 'seed', "ECl"); end %#ok<AGROW>
end
for s = 1:height(T)
    for p = params
        for sd = seeds, tasks(end+1) = struct('kind', "single", 'truth_id', T.truth_id(s), 'param', p, 'seed', sd); end %#ok<AGROW>
    end
end
for s = 1:height(T)
    for p = params, tasks(end+1) = struct('kind', "control", 'truth_id', T.truth_id(s), 'param', p, 'seed', "ECl"); end %#ok<AGROW>
end
for s = 1:height(T)
    for p = params, tasks(end+1) = struct('kind', "twostage", 'truth_id', T.truth_id(s), 'param', p, 'seed', "ECl"); end %#ok<AGROW>
end
for s = 1:height(T)
    tasks(end+1) = struct('kind', "joint", 'truth_id', T.truth_id(s), 'param', "gamma1+lambda1+delta+epsilon", 'seed', "ECl"); %#ok<AGROW>
end

for k = worker:n_workers:numel(tasks)
    t = tasks(k);
    if ~ismember(t.kind, kinds), continue; end
    truth_file = @(m) fullfile(here, 'outputs_consistent', "rec_l2f__" + m + "__" + t.seed + "__" + t.truth_id + ".mat");
    fprintf('\n=== task %d of %d: %s %s %s %s ===\n', k, numel(tasks), t.kind, t.truth_id, t.param, t.seed);
    switch t.kind
        case "scan"
            scan_one(here, out_root, T(T.truth_id == t.truth_id, :), t.param, truth_file);
        case {"single", "control", "joint", "twostage"}
            truth = "numeric";
            if t.kind == "control", truth = "mlp"; end
            out_dir = fullfile(out_root, t.seed, truth + "_truth");
            setenv('NM_START', '');  setenv('NM_STEP', '');
            if t.kind == "twostage"
                out_dir = out_dir + "_twostage";
                S = load(fullfile(out_root, 'scan', "scan__" + t.truth_id + "__" + t.param + ".mat"), 'value', 'nmse_numeric');
                [~, k0] = min(S.nmse_numeric);
                is_log = ismember(t.param, ["gamma1", "lambda1"]);
                u = S.value;  if is_log, u = log10(u); end
                setenv('NM_START', sprintf('%.12g', S.value(k0)));
                setenv('NM_STEP', sprintf('%.12g', (u(2) - u(1)) / 2));
            end
            f = fullfile(out_dir, "nm_l2f__mlp__" + t.truth_id + "__" + t.param + ".mat");
            if isfile(f) && ~isempty(whos('-file', f, 'result')), fprintf('done already\n'); continue; end
            setenv('NM_TRUTH_FILE', truth_file(truth));
            setenv('NM_OUT_DIR', out_dir);
            setenv('NM_MAX_EVALS', '');
            nelder_mead_network_recovery(t.truth_id, split(t.param, "+")', 'mlp', 'l2f');
    end
end
end

function scan_one(here, out_root, row, param, truth_file)
params = ["gamma1", "lambda1", "delta", "epsilon"];
param_range = [5e-5 5e-2; 0.01 0.1; 10 100; 10 100];
n_points = 25;
out_dir = fullfile(out_root, 'scan');
if ~isfolder(out_dir), mkdir(out_dir); end
f = fullfile(out_dir, "scan__" + row.truth_id + "__" + param + ".mat");
if isfile(f), fprintf('done already\n'); return; end
p_idx = find(params == param);
if p_idx <= 2
    value = logspace(log10(param_range(p_idx, 1)), log10(param_range(p_idx, 2)), n_points);
else
    value = linspace(param_range(p_idx, 1), param_range(p_idx, 2), n_points);
end
num = load(truth_file("numeric"), 'N', 'meta');
own = load(truth_file("mlp"), 'N');
cfg = sntm_config();
addpath(fullfile(cfg.sim_root, 'Code', 'Transport_Network'));
sur = make_surrogate('mlp', cfg);
theta = [row.gamma1, row.lambda1, row.delta, row.epsilon];
N = cell(1, n_points);  nmse_numeric = NaN(1, n_points);  nmse_mlp = NaN(1, n_points);  runtime_s = NaN(1, n_points);
for j = 1:n_points
    th = theta;  th(p_idx) = value(j);
    p = struct('gamma1', th(1), 'lambda1', th(2), 'lambda2', num.meta.lambda2, 'delta', th(3), 'epsilon', th(4));
    t0 = tic;
    out = sntm_simulate(cfg, sur, true, p, num.meta.seed, num.meta.subset, num.meta.init_rescale, j);
    runtime_s(j) = toc(t0);
    N{j} = out.Predicted.N;
    nmse_numeric(j) = nm_loss(N{j}, num.N, "nmse");
    nmse_mlp(j) = nm_loss(N{j}, own.N, "nmse");
    fprintf('  %s = %-10.4g NMSE numeric %.4g, MLP %.4g (%.1f s)\n', param, value(j), nmse_numeric(j), nmse_mlp(j), runtime_s(j));
end
truth_id = row.truth_id;  true_value = theta(p_idx);
save(f, 'truth_id', 'param', 'true_value', 'value', 'nmse_numeric', 'nmse_mlp', 'N', 'runtime_s', '-v7.3');
end

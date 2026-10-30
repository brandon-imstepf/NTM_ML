function nelder_mead_network_recovery(truth_id, param, method, conv)
% Parameter recovery on the network with Nelder-Mead (fminsearch): the
% parameters in param are searched together, the others are fixed at their
% true values. Each loss evaluation is one surrogate network simulation, run
% here, one after the other, as fminsearch asks for them.
%
% Loss = NMSE of the surrogate's N against the numeric truth's N over all
% regions and time points, sum (s - n)^2 / sum n^2 (normalized by the sum of
% squares, not the variance of Eq. 7 of the paper), the same loss as
% analyze_network_recovery.m.
% The search is in log10 for gamma1 and lambda1 and linear for delta and
% epsilon, starting from the middle of the training range. fminsearch has no
% bounds: a point outside the range is simulated at the nearest bound and
% its loss gets + (distance outside / range width), which pushes the
% simplex back inside.
%
%   nelder_mead_network_recovery('g0.0005_l0.05_d20_e80', 'lambda1')              % MLP, l2f
%   nelder_mead_network_recovery('g0.0005_l0.05_d20_e80', 'lambda1', 'mlp', 'l2eq')
%   nelder_mead_network_recovery('g0.0005_l0.05_d20_e80', ["gamma1", "delta"])     % two together
%
% truth_id  a truth of manifests/recovery_truths.csv
% param     gamma1, lambda1, delta or epsilon, or several of them (string array)
% method    make_surrogate method name (default mlp)
% conv      l2f (lambda2 = 0.01) or l2eq (lambda2 = lambda1); default l2f
%
% Writes outputs_nm/nm_<conv>__<method>__<truth_id>__<param>.mat (several
% parameters joined by +) with
%   trace   every simulation: search coordinates u, parameter values, NMSE, N, runtime
%   best    N, F, W1, M of the lowest-NMSE simulation
%   result  recovered values, NMSE, fminsearch exit flag and output, when finished
% NM_TRUTH_FILE (a .mat with N and meta) replaces the truth run, NM_OUT_DIR the
% output folder, NM_MAX_EVALS the simulation limit. NM_START sets the start
% values of the parameters (default: middle of the range) and NM_STEP the first
% steps of the search in the search coordinates (decades for gamma1 and
% lambda1; default 5% of the start coordinate), one comma-separated value per
% parameter; with NM_STEP the stopping width is also at most a tenth of the step.
% SNTM_T_END (sntm_config.m) shortens the simulations; the truth must have the
% same time points.
% The file is saved after every simulation. Rerunning after a timeout replays
% the search from the saved simulations (fminsearch is deterministic), so only
% the remaining ones are run.

if nargin < 3 || isempty(method), method = 'mlp'; end
if nargin < 4 || isempty(conv), conv = 'l2f'; end
truth_id = string(truth_id);  param = string(param);  method = char(method);  conv = string(conv);
param_label = strjoin(param, '+');

%% 1. Settings
param_names = ["gamma1", "lambda1", "delta", "epsilon"];
param_log = [true, true, false, false];
param_range = [5e-5 5e-2; 0.01 0.1; 10 100; 10 100];   % training range
max_evals = 40 * numel(param);   % simulations at most (about 5 min each for MLP); NM_MAX_EVALS overrides (tests)
if ~isempty(getenv('NM_MAX_EVALS')), max_evals = str2double(getenv('NM_MAX_EVALS')); end
tol_x_frac = 1e-3;       % stop when the simplex is narrower than this fraction of the range ...
tol_fun = 1e-6;          % ... and its NMSE values differ by less than this

%% 2. Truth and surrogate
here = fileparts(mfilename('fullpath'));
cfg = sntm_config();
addpath(fullfile(cfg.sim_root, 'Code', 'Transport_Network'));
addpath(here);
outdir = fullfile(here, 'outputs_nm');
if ~isempty(getenv('NM_OUT_DIR')), outdir = getenv('NM_OUT_DIR'); end
if ~isfolder(outdir), mkdir(outdir); end
out_file = fullfile(outdir, sprintf('nm_%s__%s__%s__%s.mat', conv, method, truth_id, param_label));

T = readtable(fullfile(here, 'manifests', 'recovery_truths.csv'), 'TextType', 'string');
T = T(T.convention == conv & T.truth_id == truth_id, :);
assert(height(T) == 1, 'truth %s not found for %s', truth_id, conv);
truth_file = fullfile(here, T.truth_dir, T.truth_case_id + ".mat");
% NM_TRUTH_FILE replaces the truth run (tests, e.g. the surrogate's own run at
% the true values); the true values still come from recovery_truths.csv
if ~isempty(getenv('NM_TRUTH_FILE')), truth_file = getenv('NM_TRUTH_FILE'); end
truth = load(truth_file, 'N', 'meta');
theta_true = [T.gamma1, T.lambda1, T.delta, T.epsilon];
assert(all(ismember(param, param_names)), 'param must be among %s', strjoin(param_names, ', '));
i = arrayfun(@(p) find(param_names == p), param);   % searched parameters

sur = make_surrogate(method, cfg);

is_log = param_log(i);
lo = toU(param_range(i,1)', is_log);  hi = toU(param_range(i,2)', is_log);
span = hi - lo;
u_start = (lo + hi) / 2;
if ~isempty(getenv('NM_START')), u_start = toU(envList('NM_START', numel(i)), is_log); end
% fminsearch searches w, with u = u_start + (w - 1) .* step / 0.05 and w = 1 at
% the start: its first simplex points (w = 1.05 in one coordinate) are then
% u_start + step in that coordinate. By default step is 5% of u_start, which is
% fminsearch's own first simplex in u
step = 0.05 * u_start;
if ~isempty(getenv('NM_STEP')), step = envList('NM_STEP', numel(i)); end
u_of_w = @(w) u_start + (w - 1) .* step / 0.05;
tol_x = tol_x_frac * span;
if ~isempty(getenv('NM_STEP')), tol_x = min(tol_x, abs(step) / 10); end   % finer than the first step

%% 3. Search
trace = struct('u', {}, 'value', {}, 'nmse', {}, 'N', {}, 'runtime_s', {});
best = struct('nmse', Inf);
if isfile(out_file)
    S = load(out_file, 'trace', 'best');
    trace = S.trace;  best = S.best;
    fprintf('[%s] resuming with %d saved simulations\n', param_label, numel(trace));
end

x_start = toX(u_start, is_log);
for k = 1:numel(i)
    fprintf('[%s %s %s] true %s = %g, start %g, range [%g, %g]\n', conv, method, truth_id, param(k), ...
        theta_true(i(k)), x_start(k), param_range(i(k),:));
end
options = optimset('Display', 'iter', 'MaxFunEvals', max_evals, 'MaxIter', max_evals, ...
    'TolX', min(tol_x * 0.05 ./ abs(step)), 'TolFun', tol_fun);   % TolX in w
[w_hat, nmse_hat, exitflag, output] = fminsearch(@objective, ones(1, numel(i)), options);
u_hat = u_of_w(w_hat);

result = struct();
result.truth_id = truth_id;  result.truth_file = string(truth_file);  result.conv = conv;  result.method = method;  result.param = param;
result.true_value = theta_true(i);
result.start_value = x_start;
result.recovered_value = toX(min(max(u_hat, lo), hi), is_log);
result.nmse = nmse_hat;
result.exitflag = exitflag;
result.output = output;
result.n_simulations = numel(trace);
result.settings = struct('max_evals', max_evals, 'tol_x_frac', tol_x_frac, 'tol_fun', tol_fun, ...
    'param_range', param_range(i,:), 'log', is_log, 'first_step', step, 'tol_x', tol_x, ...
    'trange', cfg.trange);
save(out_file, 'trace', 'best', 'result', '-v7.3');
fprintf('[%s] recovered %s (true %s), NMSE %.4g, %d simulations, exitflag %d\n', param_label, ...
    mat2str(result.recovered_value, 6), mat2str(result.true_value, 6), nmse_hat, result.n_simulations, exitflag);

    function loss = objective(w)
        u = u_of_w(w);
        u_in = min(max(u, lo), hi);
        loss = simulateNMSE(u_in) + sum(abs(u - u_in) ./ span);
    end

    function nmse = simulateNMSE(u)
        % a point already simulated (this run or before a restart) is not simulated again
        if ~isempty(trace)
            k = find(all(abs(vertcat(trace.u) - u) <= 1e-12 * span, 2), 1);
            if ~isempty(k)
                nmse = trace(k).nmse;
                return
            end
        end
        theta = theta_true;
        theta(i) = toX(u, is_log);
        p = struct('gamma1', theta(1), 'lambda1', theta(2), 'lambda2', truth.meta.lambda2, ...
            'delta', theta(3), 'epsilon', theta(4));
        if conv == "l2eq", p.lambda2 = theta(2); end
        t0 = tic;
        out = sntm_simulate(cfg, sur, true, p, truth.meta.seed, truth.meta.subset, ...
            truth.meta.init_rescale, numel(trace) + 1);
        runtime_s = toc(t0);
        N = out.Predicted.N;
        assert(isequal(size(N), size(truth.N)), 'simulation and truth have different time points (SNTM_T_END?)');
        nmse = sum((N(:) - truth.N(:)).^2) / sum(truth.N(:).^2);
        if ~isfinite(nmse), nmse = Inf; end   % failed simulation: worst possible loss
        trace(end+1) = struct('u', u, 'value', theta(i), 'nmse', nmse, 'N', N, 'runtime_s', runtime_s);
        if nmse < best.nmse
            best = struct('nmse', nmse, 'value', theta(i), 'N', N, 'F', out.Predicted.F, ...
                'W1', out.Predicted.W1, 'M', out.Predicted.M);
        end
        save(out_file, 'trace', 'best', '-v7.3');
        fprintf('  sim %2d: %s  NMSE = %.5g  (%.0f s)\n', numel(trace), ...
            strjoin(compose('%s = %-10.4g', param, theta(i)), '  '), nmse, runtime_s);
    end
end

function u = toU(x, is_log)
% parameter values to search coordinates (log10 where is_log)
    u = x;
    u(is_log) = log10(x(is_log));
end

function x = toX(u, is_log)
    x = u;
    x(is_log) = 10.^u(is_log);
end

function v = envList(name, n)
% comma-separated numbers of an environment variable, one per searched parameter
    v = str2double(split(string(getenv(name)), ","))';
    assert(numel(v) == n && all(isfinite(v)), '%s needs %d comma-separated numbers', name, n);
end

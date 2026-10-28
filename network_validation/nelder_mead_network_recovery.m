function nelder_mead_network_recovery(truth_id, param, method, conv)
% Single-parameter recovery on the network with Nelder-Mead (fminsearch):
% one parameter is searched, the other three are fixed at their true values.
% Each loss evaluation is one surrogate network simulation, run here, one
% after the other, as fminsearch asks for them.
%
% Loss = NMSE of the surrogate's N against the numeric truth's N over all
% regions and time points, sum (s - n)^2 / sum n^2 (Eq. 7 of the paper), the
% same loss as analyze_network_recovery.m.
% The search is in log10 for gamma1 and lambda1 and linear for delta and
% epsilon, starting from the middle of the training range. fminsearch has no
% bounds: a point outside the range is simulated at the nearest bound and
% its loss gets + (distance outside / range width), which pushes the
% simplex back inside.
%
%   nelder_mead_network_recovery('g0.0005_l0.05_d20_e80', 'lambda1')              % MLP, l2f
%   nelder_mead_network_recovery('g0.0005_l0.05_d20_e80', 'lambda1', 'mlp', 'l2eq')
%
% truth_id  a truth of manifests/recovery_truths.csv
% param     gamma1, lambda1, delta or epsilon
% method    make_surrogate method name (default mlp)
% conv      l2f (lambda2 = 0.01) or l2eq (lambda2 = lambda1); default l2f
%
% Writes outputs_nm/nm_<conv>__<method>__<truth_id>__<param>.mat with
%   trace   every simulation: search coordinate u, parameter value, NMSE, N, runtime
%   best    N, F, W1, M of the lowest-NMSE simulation
%   result  recovered value, NMSE, fminsearch exit flag and output, when finished
% The file is saved after every simulation. Rerunning after a timeout replays
% the search from the saved simulations (fminsearch is deterministic), so only
% the remaining ones are run.

if nargin < 3 || isempty(method), method = 'mlp'; end
if nargin < 4 || isempty(conv), conv = 'l2f'; end
truth_id = string(truth_id);  param = string(param);  method = char(method);  conv = string(conv);

%% 1. Settings
param_names = ["gamma1", "lambda1", "delta", "epsilon"];
param_log = [true, true, false, false];
param_range = [5e-5 5e-2; 0.01 0.1; 10 100; 10 100];   % training range
max_evals = 40;          % simulations at most (about 5 min each for MLP); NM_MAX_EVALS overrides (tests)
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
out_file = fullfile(outdir, sprintf('nm_%s__%s__%s__%s.mat', conv, method, truth_id, param));

T = readtable(fullfile(here, 'manifests', 'recovery_truths.csv'), 'TextType', 'string');
T = T(T.convention == conv & T.truth_id == truth_id, :);
assert(height(T) == 1, 'truth %s not found for %s', truth_id, conv);
truth = load(fullfile(here, T.truth_dir, T.truth_case_id + ".mat"), 'N', 'meta');
theta_true = [T.gamma1, T.lambda1, T.delta, T.epsilon];
i = find(param_names == param);
assert(~isempty(i), 'param must be one of %s', strjoin(param_names, ', '));

sur = make_surrogate(method, cfg);

u_of = @(x) x;  x_of = @(u) u;
if param_log(i), u_of = @(x) log10(x);  x_of = @(u) 10.^u; end
lo = u_of(param_range(i,1));  hi = u_of(param_range(i,2));
span = hi - lo;
u_start = (lo + hi) / 2;

%% 3. Search
trace = struct('u', {}, 'value', {}, 'nmse', {}, 'N', {}, 'runtime_s', {});
best = struct('nmse', Inf);
if isfile(out_file)
    S = load(out_file, 'trace', 'best');
    trace = S.trace;  best = S.best;
    fprintf('[%s] resuming with %d saved simulations\n', param, numel(trace));
end

fprintf('[%s %s %s] true %s = %g, start %g, range [%g, %g]\n', conv, method, truth_id, param, ...
    theta_true(i), x_of(u_start), param_range(i,:));
options = optimset('Display', 'iter', 'MaxFunEvals', max_evals, 'MaxIter', max_evals, ...
    'TolX', tol_x_frac * span, 'TolFun', tol_fun);
[u_hat, nmse_hat, exitflag, output] = fminsearch(@objective, u_start, options);

result = struct();
result.truth_id = truth_id;  result.conv = conv;  result.method = method;  result.param = param;
result.true_value = theta_true(i);
result.start_value = x_of(u_start);
result.recovered_value = x_of(min(max(u_hat, lo), hi));
result.nmse = nmse_hat;
result.exitflag = exitflag;
result.output = output;
result.n_simulations = numel(trace);
result.settings = struct('max_evals', max_evals, 'tol_x_frac', tol_x_frac, 'tol_fun', tol_fun, ...
    'param_range', param_range(i,:), 'log', param_log(i));
save(out_file, 'trace', 'best', 'result', '-v7.3');
fprintf('[%s] recovered %g (true %g), NMSE %.4g, %d simulations, exitflag %d\n', ...
    param, result.recovered_value, result.true_value, nmse_hat, result.n_simulations, exitflag);

    function loss = objective(u)
        u_in = min(max(u, lo), hi);
        loss = simulateNMSE(u_in) + abs(u - u_in) / span;
    end

    function nmse = simulateNMSE(u)
        % a value already simulated (this run or before a restart) is not simulated again
        k = find(abs([trace.u] - u) <= 1e-12 * span, 1);
        if ~isempty(k)
            nmse = trace(k).nmse;
            return
        end
        theta = theta_true;
        theta(i) = x_of(u);
        p = struct('gamma1', theta(1), 'lambda1', theta(2), 'lambda2', truth.meta.lambda2, ...
            'delta', theta(3), 'epsilon', theta(4));
        if conv == "l2eq", p.lambda2 = theta(2); end
        t0 = tic;
        out = sntm_simulate(cfg, sur, true, p, truth.meta.seed, truth.meta.subset, ...
            truth.meta.init_rescale, numel(trace) + 1);
        runtime_s = toc(t0);
        N = out.Predicted.N;
        nmse = sum((N(:) - truth.N(:)).^2) / sum(truth.N(:).^2);
        if ~isfinite(nmse), nmse = Inf; end   % failed simulation: worst possible loss
        trace(end+1) = struct('u', u, 'value', theta(i), 'nmse', nmse, 'N', N, 'runtime_s', runtime_s);
        if nmse < best.nmse
            best = struct('nmse', nmse, 'value', theta(i), 'N', N, 'F', out.Predicted.F, ...
                'W1', out.Predicted.W1, 'M', out.Predicted.M);
        end
        save(out_file, 'trace', 'best', '-v7.3');
        fprintf('  sim %2d: %s = %-10.4g NMSE = %.5g  (%.0f s)\n', numel(trace), param, theta(i), nmse, runtime_s);
    end
end

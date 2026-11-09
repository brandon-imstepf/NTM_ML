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
% NM_RESTARTS (default 1) runs that many searches and keeps the one with the
% lowest loss: the first from the start above, the others from points drawn
% uniformly over the range in the search coordinates (seed = start number, so
% a rerun draws the same points); result.starts has every search.
% NM_LOSS (default nmse) is the loss, see nm_loss.m; other than nmse it is
% added to the file name (__loss-<NM_LOSS>). trace.nmse is always the NMSE.
% Each search's first step is NM_STEP if set, otherwise 5% of its start coordinate.
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
n_starts = 1;
if ~isempty(getenv('NM_RESTARTS')), n_starts = str2double(getenv('NM_RESTARTS')); end
loss_kind = "nmse";
if ~isempty(getenv('NM_LOSS')), loss_kind = string(getenv('NM_LOSS')); end

%% 2. Truth and surrogate
here = fileparts(mfilename('fullpath'));
cfg = sntm_config();
addpath(fullfile(cfg.sim_root, 'Code', 'Transport_Network'));
addpath(here);
outdir = fullfile(here, 'outputs_nm');
if ~isempty(getenv('NM_OUT_DIR')), outdir = getenv('NM_OUT_DIR'); end
if ~isfolder(outdir), mkdir(outdir); end
loss_label = "";
if loss_kind ~= "nmse", loss_label = "__loss-" + loss_kind; end
out_file = fullfile(outdir, sprintf('nm_%s__%s__%s__%s%s.mat', conv, method, truth_id, param_label, loss_label));

T = readtable(fullfile(here, 'manifests', 'recovery_truths.csv'), 'TextType', 'string');
T = T(T.convention == conv & T.truth_id == truth_id, :);
assert(height(T) == 1, 'truth %s not found for %s', truth_id, conv);
truth_file = fullfile(here, T.truth_dir, T.truth_case_id + ".mat");
% NM_TRUTH_FILE replaces the truth run (tests, e.g. the surrogate's own run at
% the true values); the true values still come from recovery_truths.csv
if ~isempty(getenv('NM_TRUTH_FILE')), truth_file = getenv('NM_TRUTH_FILE'); end
truth = load(truth_file, 'N', 'meta');
day = cfg.trange * 180;
nm_loss(truth.N, truth.N, loss_kind, day);   % an unknown NM_LOSS fails here, before any simulation
theta_true = [T.gamma1, T.lambda1, T.delta, T.epsilon];
assert(all(ismember(param, param_names)), 'param must be among %s', strjoin(param_names, ', '));
i = arrayfun(@(p) find(param_names == p), param);   % searched parameters

sur = make_surrogate(method, cfg);

is_log = param_log(i);
lo = toU(param_range(i,1)', is_log);  hi = toU(param_range(i,2)', is_log);
span = hi - lo;
u_start_1 = (lo + hi) / 2;
if ~isempty(getenv('NM_START')), u_start_1 = toU(envList('NM_START', numel(i)), is_log); end

%% 3. Search
trace = struct('u', {}, 'value', {}, 'nmse', {}, 'loss', {}, 'N', {}, 'runtime_s', {});
best = struct('nmse', Inf, 'loss', Inf);
if isfile(out_file)
    S = load(out_file, 'trace', 'best');
    trace = S.trace;  best = S.best;
    if ~isfield(trace, 'loss')   % saved before NM_LOSS: the loss was the NMSE
        assert(loss_kind == "nmse", 'saved file %s has no loss field', out_file);
        [trace.loss] = trace.nmse;  best.loss = best.nmse;
        trace = orderfields(trace, {'u', 'value', 'nmse', 'loss', 'N', 'runtime_s'});
    end
    fprintf('[%s] resuming with %d saved simulations\n', param_label, numel(trace));
end

starts = struct('start_value', {}, 'recovered_value', {}, 'loss', {}, 'nmse', {}, ...
    'exitflag', {}, 'output', {}, 'first_step', {}, 'tol_x', {});
for r = 1:n_starts
    u_start = u_start_1;
    if r > 1
        rs = RandStream('mt19937ar', 'Seed', r);
        u_start = lo + rand(rs, 1, numel(i)) .* span;
    end
    % fminsearch searches w, with u = u_start + (w - 1) .* step / 0.05 and w = 1 at
    % the start: its first simplex points (w = 1.05 in one coordinate) are then
    % u_start + step in that coordinate. By default step is 5% of u_start, which is
    % fminsearch's own first simplex in u
    step = 0.05 * u_start;
    if ~isempty(getenv('NM_STEP')), step = envList('NM_STEP', numel(i)); end
    u_of_w = @(w) u_start + (w - 1) .* step / 0.05;
    tol_x = tol_x_frac * span;
    if ~isempty(getenv('NM_STEP')), tol_x = min(tol_x, abs(step) / 10); end   % finer than the first step

    x_start = toX(u_start, is_log);
    for k = 1:numel(i)
        fprintf('[%s %s %s] start %d of %d: true %s = %g, start %g, range [%g, %g], loss %s\n', conv, method, ...
            truth_id, r, n_starts, param(k), theta_true(i(k)), x_start(k), param_range(i(k),:), loss_kind);
    end
    options = optimset('Display', 'iter', 'MaxFunEvals', max_evals, 'MaxIter', max_evals, ...
        'TolX', min(tol_x * 0.05 ./ abs(step)), 'TolFun', tol_fun);   % TolX in w
    [w_hat, loss_hat, exitflag, output] = fminsearch(@objective, ones(1, numel(i)), options);
    u_hat = min(max(u_of_w(w_hat), lo), hi);
    starts(r) = struct('start_value', x_start, 'recovered_value', toX(u_hat, is_log), 'loss', loss_hat, ...
        'nmse', lookupNMSE(u_hat), 'exitflag', exitflag, 'output', output, 'first_step', step, 'tol_x', tol_x);
end
[~, r_best] = min([starts.loss]);

result = struct();
result.truth_id = truth_id;  result.truth_file = string(truth_file);  result.conv = conv;  result.method = method;  result.param = param;
result.true_value = theta_true(i);
result.start_value = starts(r_best).start_value;
result.recovered_value = starts(r_best).recovered_value;
result.loss_kind = loss_kind;
result.loss = starts(r_best).loss;
result.nmse = starts(r_best).nmse;
result.exitflag = starts(r_best).exitflag;
result.output = starts(r_best).output;
result.best_start = r_best;
result.starts = starts;
result.n_simulations = numel(trace);
result.settings = struct('max_evals', max_evals, 'tol_x_frac', tol_x_frac, 'tol_fun', tol_fun, ...
    'param_range', param_range(i,:), 'log', is_log, 'first_step', starts(r_best).first_step, ...
    'tol_x', starts(r_best).tol_x, 'trange', cfg.trange, 'n_starts', n_starts, 'loss', loss_kind);
save(out_file, 'trace', 'best', 'result', '-v7.3');
for r = 1:n_starts
    fprintf('  start %d: from %s to %s, loss %.4g, exitflag %d\n', r, mat2str(starts(r).start_value, 4), ...
        mat2str(starts(r).recovered_value, 6), starts(r).loss, starts(r).exitflag);
end
fprintf('[%s] recovered %s (true %s), loss %.4g, NMSE %.4g, start %d of %d, %d simulations, exitflag %d\n', ...
    param_label, mat2str(result.recovered_value, 6), mat2str(result.true_value, 6), result.loss, ...
    result.nmse, r_best, n_starts, result.n_simulations, result.exitflag);

    function loss = objective(w)
        u = u_of_w(w);
        u_in = min(max(u, lo), hi);
        loss = simulateLoss(u_in) + sum(abs(u - u_in) ./ span);
    end

    function nmse = lookupNMSE(u)
        % NMSE of a simulated point (the recovered point always is one)
        nmse = NaN;
        k = find(all(abs(vertcat(trace.u) - u) <= 1e-12 * span, 2), 1);
        if ~isempty(k), nmse = trace(k).nmse; end
    end

    function loss = simulateLoss(u)
        % a point already simulated (this run or before a restart) is not simulated again
        if ~isempty(trace)
            k = find(all(abs(vertcat(trace.u) - u) <= 1e-12 * span, 2), 1);
            if ~isempty(k)
                loss = trace(k).loss;
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
        nmse = nm_loss(N, truth.N, "nmse");
        loss = nm_loss(N, truth.N, loss_kind, day);
        if ~isfinite(nmse), nmse = Inf; end   % failed simulation: worst possible loss
        if ~isfinite(loss), loss = Inf; end
        trace(end+1) = struct('u', u, 'value', theta(i), 'nmse', nmse, 'loss', loss, 'N', N, 'runtime_s', runtime_s);
        if loss < best.loss
            best = struct('nmse', nmse, 'loss', loss, 'value', theta(i), 'N', N, 'F', out.Predicted.F, ...
                'W1', out.Predicted.W1, 'M', out.Predicted.M);
        end
        save(out_file, 'trace', 'best', '-v7.3');
        fprintf('  sim %2d: %s  loss = %.5g  NMSE = %.5g  (%.0f s)\n', numel(trace), ...
            strjoin(compose('%s = %-10.4g', param(:), reshape(theta(i), [], 1)), '  '), loss, nmse, runtime_s);
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

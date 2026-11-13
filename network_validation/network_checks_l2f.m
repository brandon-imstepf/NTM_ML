function A = network_checks_l2f(run_set)
% Network-level accuracy and behaviour of each surrogate against the numeric
% truth, for the 8 l2f truths of manifests/recovery_truths.csv (lambda2 = 0.01 in
% the numeric runs, as in the surrogates' training data; 30 regions, 180 days),
% for the reviewer responses (R1.m1, R2.M4). Seed ECl from outputs_recovery/;
% seeds CA1 and RSPd (manifests/seeds_l2f_*.csv) from outputs_seeds/, where present.
%   nmse         NMSE (Eq. 7, nmse_eq7.m): sum (s - n)^2 / sum (n - mean n)^2 over
%                all regions and time points
%   r2           1 - nmse (R^2 against y = x over the same points)
%   nmse_last30  nmse over days 150-180
%   sse_rel, sse_rel_last30  sum (s - n)^2 / sum n^2 (sse_rel.m), the metric
%                reported before 2026-11-10
%   min_N        smallest N of the run (negative = loss of non-negativity)
%   n_negative   number of N values below 0
%   n_clamped    number of N values equal to eps: NetworkTransportModel replaces
%                negative N by eps after each step, so this counts the values
%                that went negative
%   finite       all N finite
%   rel_change_last  max over regions of |N(180) - N(171)| / max N, the change
%                over the last time step (long-term stability)
%   runtime_s    MATLAB tic/toc of the run
% Time-step sensitivity of the MLP: NMSE between the 60-day MLP runs at the
% default grid and at a 0.225-day maximum step (outputs_t60_*), at the shared
% time points up to day 59.
% Writes results/network_checks_l2f.csv and results/network_dt_mlp_l2f.csv.
% run_set (l2f_runs.m): "original" (default) or "consistent", which reads
% outputs_consistent/, writes results/network_checks_l2f_consistent.csv and
% skips the time-step check (run in the original convention only).
%
%   A = network_checks_l2f();
%   A = network_checks_l2f("consistent");

if nargin < 1, run_set = "original"; end
here = fileparts(mfilename('fullpath'));
T = readtable(fullfile(here, 'manifests', 'recovery_truths.csv'), 'TextType', 'string');
T = T(T.convention == "l2f", :);
methods = ["numeric", "mlp", "sr", "gp", "lr"];
file = @(dir, m, id) fullfile(here, dir, "rec_l2f__" + m + "__ECl__" + id + ".mat");

seeds = ["ECl", "CA1", "RSPd"];
[seed_dirs, tag] = l2f_runs(run_set);
seed_file = @(k, m, id) fullfile(here, seed_dirs(k), "rec_l2f__" + m + "__" + seeds(k) + "__" + id + ".mat");

rows = {};
for k = 1:numel(seeds)
for s = 1:height(T)
    if ~all(arrayfun(@(m) isfile(seed_file(k, m, T.truth_id(s))), methods)), continue, end
    num = load(seed_file(k, "numeric", T.truth_id(s)), 'N', 'trange');
    day = num.trange * 180;
    last30 = day >= 150;
    for m = methods
        r = load(seed_file(k, m, T.truth_id(s)), 'N', 'runtime_s');
        N = r.N;  n = num.N;
        nmse = nmse_eq7(N, n);
        r2 = 1 - nmse;
        nmse_last = nmse_eq7(N(:, last30), n(:, last30));
        rel_change = max(abs(N(:, end) - N(:, end-1))) / max(N(:));
        rows(end+1, :) = {seeds(k), T.truth_id(s), m, nmse, r2, nmse_last, min(N(:)), sum(N(:) < 0), ...
            sum(N(:) == eps), all(isfinite(N(:))), rel_change, r.runtime_s, sse_rel(N, n), ...
            sse_rel(N(:, last30), n(:, last30))}; %#ok<AGROW>
    end
end
end
A = cell2table(rows, 'VariableNames', {'seed', 'truth_id', 'method', 'nmse', 'r2', 'nmse_last30', 'min_N', ...
    'n_negative', 'n_clamped', 'finite', 'rel_change_last', 'runtime_s', 'sse_rel', 'sse_rel_last30'});
writetable(A, fullfile(here, 'results', "network_checks_l2f" + tag + ".csv"));

fprintf('== l2f truths, 180 days, seeds %s: median [min, max] over %d simulations per method\n', ...
    strjoin(unique(A.seed, 'stable'), ', '), sum(A.method == "numeric"));
fprintf('%-8s %-28s %-28s %-28s %8s %8s %-22s\n', 'method', 'NMSE', 'R^2', 'NMSE days 150-180', 'neg. N', 'finite', 'runtime (s)');
for m = methods
    a = A(A.method == m, :);
    fprintf('%-8s %-28s %-28s %-28s %8d %5d/%d %-22s\n', m, mmm(a.nmse), mmm(a.r2), mmm(a.nmse_last30), ...
        sum(a.n_negative), sum(a.finite), height(a), mmm(a.runtime_s));
end
fprintf('N clamped from negative to eps (simulations with any / most values in one simulation): %s\n', ...
    strjoin(arrayfun(@(m) sprintf('%s %d/%d, %d', m, nnz(A.n_clamped(A.method == m) > 0), ...
    nnz(A.method == m), max(A.n_clamped(A.method == m))), methods, 'UniformOutput', false), '; '));
fprintf('largest change over the last step (fraction of max N): %s\n', ...
    strjoin(compose('%s %.2g', methods', arrayfun(@(m) max(A.rel_change_last(A.method == m)), methods')), ', '));

%% Time-step sensitivity of the MLP (60 days)
if run_set ~= "original", return, end
rows = {};
for s = 1:height(T)
    b = load(file('outputs_t60_base', "mlp", T.truth_id(s)), 'N', 'trange');
    f = load(file('outputs_t60_dt0.225', "mlp", T.truth_id(s)), 'N', 'trange');
    nb = load(file('outputs_t60_base', "numeric", T.truth_id(s)), 'N');
    nf = load(file('outputs_t60_dt0.225', "numeric", T.truth_id(s)), 'N');
    [~, ib, jf] = intersect(round(b.trange * 180, 6), round(f.trange * 180, 6));
    keep = b.trange(ib) * 180 <= 59;
    ib = ib(keep);  jf = jf(keep);
    d = @nmse_eq7;
    rows(end+1, :) = {T.truth_id(s), d(b.N(:, ib), f.N(:, jf)), d(nb.N(:, ib), nf.N(:, jf)), ...
        d(f.N, nf.N)}; %#ok<AGROW>
end
D = cell2table(rows, 'VariableNames', {'truth_id', 'mlp_base_vs_dt0225', 'numeric_base_vs_dt0225', 'mlp_vs_numeric_dt0225'});
writetable(D, fullfile(here, 'results', 'network_dt_mlp_l2f.csv'));
fprintf('\n== Time step (60 days, default grid vs 0.225-day steps, shared points to day 59)\n');
disp(D)
end

function s = mmm(x)
% non-finite values (a failed run) are left out
s = sprintf('%.3g [%.3g, %.3g]', median(x, 'omitnan'), min(x), max(x));
end

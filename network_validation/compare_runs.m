function s = compare_runs(file_a, file_b, t_end_days)
% Compare two run_sntm_cases outputs (a = reference, b = the run to check).
% Runs on different time grids are compared at the time points they share;
% t_end_days (optional) keeps only points up to that day.
% Per field (N, M, W1, F, and EdgeMass when both have it) prints:
%   isequal   1 = bit-identical
%   max abs   max |a - b|
%   max rel   max |a - b| / max |a|
%   NMSE      sum (a - b)^2 / sum a^2
%   first     first shared time point (day) where a and b differ
% and the runtime, MATLAB version, host and thread count of each run, since
% the last bits depend on the thread count (and possibly the machine).
% With an output argument nothing is printed; s is a table with one row per
% field (isequal, max_abs, max_rel, nmse, first_day) and s.Properties.UserData
% holds the number of shared time points and the last shared day.
%
%   compare_runs('outputs_t60_base/rec_l2f__numeric__ECl__g0.0005_l0.02_d80_e20.mat', ...
%                'outputs_t60_dt0.45/rec_l2f__numeric__ECl__g0.0005_l0.02_d80_e20.mat', 60)

if nargin < 3, t_end_days = inf; end
a = load(file_a);
b = load(file_b);
show = nargout == 0;

% Time points the two runs share (trange is in units of 180 days)
[in_b, loc] = ismembertol(a.trange, b.trange, 1e-9, 'DataScale', 1);
in_b = in_b & a.trange * 180 <= t_end_days + 1e-9;
ia = find(in_b);
ib = loc(in_b);
days = a.trange(ia) * 180;
if show
    fprintf('a: %s\nb: %s\n', file_a, file_b);
    fprintf('time points: a %d, b %d, shared %d (day %g to %g)\n', ...
        numel(a.trange), numel(b.trange), numel(ia), days(1), days(end));
end

fields = {'N', 'M', 'W1', 'F', 'EdgeMass'};
s = table('Size', [0 6], 'VariableTypes', {'string', 'logical', 'double', 'double', 'double', 'double'}, ...
    'VariableNames', {'field', 'isequal', 'max_abs', 'max_rel', 'nmse', 'first_day'});
for k = 1:numel(fields)
    f = fields{k};
    if ~isfield(a, f) || ~isfield(b, f), continue, end
    xa = at_times(a.(f), ia);
    xb = at_times(b.(f), ib);
    d = abs(xa - xb);
    differs = find(any(d > 0, 1), 1);
    first_day = NaN;
    if ~isempty(differs), first_day = days(differs); end
    s(end+1, :) = {f, isequal(xa, xb), max(d(:)), max(d(:)) / max(abs(xa(:))), ...
        sum(d(:).^2) / sum(xa(:).^2), first_day}; %#ok<AGROW>
    if show
        first = 'none';
        if ~isnan(first_day), first = sprintf('day %g', first_day); end
        fprintf('%-8s isequal %d  max abs %.3g  max rel %.3g  NMSE %.3g  first %s\n', f, ...
            s.isequal(end), s.max_abs(end), s.max_rel(end), s.nmse(end), first);
    end
end
s.Properties.UserData = struct('n_shared', numel(ia), 'last_day', days(end));
if ~show, return, end

if isfield(a, 'env') && isfield(b, 'env')
    fprintf('runtime %.1f s vs %.1f s (ratio %.1f)\n', a.runtime_s, b.runtime_s, a.runtime_s / b.runtime_s);
    fprintf('matlab  %s | %s\nhost    %s | %s\nthreads %d | %d\n', a.env.matlab, b.env.matlab, ...
        a.env.host, b.env.host, a.env.max_threads, b.env.max_threads);
end
if isfield(a, 'cfg') && isfield(b, 'cfg')
    fprintf('sim     %s | %s\n', a.cfg.sim_root, b.cfg.sim_root);
end
end

function x = at_times(x, idx)
% Columns idx of the last (time) dimension, one column per time point
sz = size(x);
x = reshape(x, [], sz(end));
x = x(:, idx);
end

function aggregate_sntm(out_dir)
% Compare every surrogate run with its numeric run (same setting_id) and
% summarise input-space coverage of the numeric runs.
%
% Writes to <out_dir>/summary/:
%   accuracy_by_time.csv  one row per (surrogate case, time point)
%   accuracy_by_case.csv  one row per surrogate case (worst / median over time)
%   coverage.csv          per numeric case: N reached vs. training distribution
%
% Metrics compare regional N across regions at each time point:
%   r2_identity  1 - SS_res/SS_tot against the y = x line (penalises bias)
%   slope0       least-squares slope of surrogate on numeric through the origin
%   r2_origin    R^2 of that through-origin fit
%   nmse         sum (s - n)^2 / sum n^2   (Eq. 7 of the paper)
%   max_rel_err  max_i |s_i - n_i| / max(n)
%   num_cv       std/mean of numeric N across regions (how much spatial
%                pattern there is to predict; R^2 at one time point is
%                unstable when this is near zero)
% Per case:
%   r2_pooled, nmse_pooled   all regions and time points at once
%   r2_ss_median             median over time points from SS_DAYS on of r2_identity
%   r2_ss_pooled, nmse_ss    all regions and time points from SS_DAYS on at once

if nargin < 1, cfg = sntm_config(); out_dir = cfg.out_dir; else, cfg = sntm_config(); end
sum_dir = fullfile(out_dir, 'summary');
if ~isfolder(sum_dir), mkdir(sum_dir); end
SS_DAYS = 120; % "steady state": from month 4 on

files = dir(fullfile(out_dir, '*.mat'));
meta = cell(numel(files), 1);
for k = 1:numel(files)
    S = load(fullfile(files(k).folder, files(k).name), 'meta');
    meta{k} = S.meta;
end
meta = [meta{:}];
methods = string({meta.method});
is_num = methods == "numeric";
num_ids = string({meta(is_num).setting_id});
num_files = files(is_num);

rows_t = {};
rows_c = {};
for k = find(~is_num)
    j = find(num_ids == string(meta(k).setting_id), 1);
    if isempty(j), continue, end
    s = load(fullfile(files(k).folder, files(k).name), 'N', 'trange', 'runtime_s');
    n = load(fullfile(num_files(j).folder, num_files(j).name), 'N', 'runtime_s');
    nt = min(size(s.N, 2), size(n.N, 2));
    m = per_time_metrics(s.N(:, 2:nt), n.N(:, 2:nt)); % skip t=0 (identical by construction)
    days = s.trange(2:nt).' * 180;
    base = struct2table(meta(k), 'AsArray', true);
    rows_t{end+1} = [repmat(base, numel(days), 1), table(days), struct2table(m)]; %#ok<AGROW>
    finite = isfinite(m.r2_identity);
    sv = s.N(:, 2:nt); nv = n.N(:, 2:nt);
    ss = days >= SS_DAYS;
    [r2_pooled, nmse_pooled] = pooled_metrics(sv, nv);
    [r2_ss_pooled, nmse_ss] = pooled_metrics(sv(:, ss), nv(:, ss));
    rows_c{end+1} = [base, table( ...
        r2_pooled, nmse_pooled, median(m.r2_identity(ss)), r2_ss_pooled, nmse_ss, ...
        min(m.r2_identity), median(m.r2_identity(finite)), m.r2_identity(end), ...
        median(m.slope0), max(m.nmse), m.nmse(end), any(~isfinite(s.N(:))), ...
        s.runtime_s, n.runtime_s, n.runtime_s / s.runtime_s, ...
        'VariableNames', {'r2_pooled', 'nmse_pooled', 'r2_ss_median', 'r2_ss_pooled', 'nmse_ss', ...
        'r2_min', 'r2_median', 'r2_final', 'slope0_median', 'nmse_max', ...
        'nmse_final', 'nonfinite', 'runtime_s', 'numeric_runtime_s', 'speedup'})]; %#ok<AGROW>
end
if ~isempty(rows_t)
    writetable(vertcat(rows_t{:}), fullfile(sum_dir, 'accuracy_by_time.csv'));
    writetable(vertcat(rows_c{:}), fullfile(sum_dir, 'accuracy_by_case.csv'));
end

% Coverage of edge inputs (N_i, N_j on edges with tau present) vs. training data
T = readmatrix(cfg.train_superset_file);
Ntrain = [T(:,5); T(:,6)];
Ns = sort(Ntrain); q = Ns(ceil([0.5 0.99 0.999] * numel(Ns))).'; % no Statistics Toolbox needed
rows_v = {};
for j = 1:numel(num_files)
    S = load(fullfile(num_files(j).folder, num_files(j).name), 'N', 'C', 'meta');
    [ei, ej] = find(S.C > 0);
    N1 = S.N(ei, :); N2 = S.N(ej, :);
    used = N1 > 0 | N2 > 0;
    v = [N1(used); N2(used)];
    rows_v{end+1} = [struct2table(S.meta, 'AsArray', true), table( ...
        numel(v), max(v), mean(v > q(2)), mean(v > q(3)), mean(v > max(Ntrain)), ...
        'VariableNames', {'n_edge_inputs', 'N_max', 'frac_above_train_q99', ...
        'frac_above_train_q999', 'frac_above_train_max'})]; %#ok<AGROW>
end
if ~isempty(rows_v)
    writetable(vertcat(rows_v{:}), fullfile(sum_dir, 'coverage.csv'));
end
fprintf('training N: median %.3g, q99 %.3g, q99.9 %.3g, max %.3g\n', q, max(Ntrain));
fprintf('wrote %s\n', sum_dir);
end

function m = per_time_metrics(S, N)
nt = size(N, 2);
[m.r2_identity, m.slope0, m.r2_origin, m.nmse, m.max_rel_err, m.num_cv] = deal(nan(nt, 1));
for t = 1:nt
    s = S(:, t); n = N(:, t);
    m.r2_identity(t) = 1 - sum((s - n).^2) / sum((n - mean(n)).^2);
    m.slope0(t) = (n' * s) / (n' * n);
    m.r2_origin(t) = 1 - sum((s - m.slope0(t) * n).^2) / sum(s.^2);
    m.nmse(t) = sum((s - n).^2) / sum(n.^2);
    m.max_rel_err(t) = max(abs(s - n)) / max(n);
    m.num_cv(t) = std(n) / mean(n);
end
end

function [r2, nmse] = pooled_metrics(S, N)
% R^2 against y = x and NMSE (Eq. 7) over all entries of S and N together
r2 = 1 - sum((S(:) - N(:)).^2) / sum((N(:) - mean(N(:))).^2);
nmse = sum((S(:) - N(:)).^2) / sum(N(:).^2);
end

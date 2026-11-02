function R = timing_table3(out_dir)
% Run times for Table 3 (R1.m4): the repeated runs of manifests/timing_fig8.csv
% (Fig. 8 setting, original simulator code), mean and standard error of the
% mean per method, with the machines they ran on. runtime_s is the tic/toc of
% the simulation in run_sntm_cases.m, without MATLAB start-up and surrogate loading.
% Writes results/timing_table3.csv (one row per run).
%
%   R = timing_table3();                   % outputs_timing/
%   R = timing_table3('outputs_timing');

here = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(out_dir), out_dir = fullfile(here, 'outputs_timing'); end
methods = ["numeric", "mlp", "sr", "gp", "lr"];
files = dir(fullfile(out_dir, 'timing__*.mat'));
rows = {};
for f = files'
    r = load(fullfile(f.folder, f.name), 'meta', 'runtime_s', 'env');
    rows(end+1, :) = {string(r.meta.method), string(r.meta.case_id), r.runtime_s, string(r.env.host), ...
        string(r.env.cpu), string(r.env.matlab), r.env.max_threads}; %#ok<AGROW>
end
R = cell2table(rows, 'VariableNames', {'method', 'case_id', 'runtime_s', 'host', 'cpu', 'matlab', 'threads'});
writetable(R, fullfile(here, 'results', 'timing_table3.csv'));

fprintf('== Table 3 timing, %s\n', out_dir);
fprintf('%-8s %3s %14s %10s %8s\n', 'method', 'n', 'mean (s)', 's.e.m.', 'speedup');
t_num = mean(R.runtime_s(R.method == "numeric"));
for m = methods
    t = R.runtime_s(R.method == m);
    if isempty(t), continue, end
    fprintf('%-8s %3d %14.1f %10.1f %8.1f\n', m, numel(t), mean(t), std(t) / sqrt(numel(t)), t_num / mean(t));
end
fprintf('machines: %s\ncpu: %s\nMATLAB: %s; threads: %s\n', strjoin(unique(R.host), ', '), ...
    strjoin(unique(R.cpu), '; '), strjoin(unique(R.matlab), ', '), mat2str(unique(R.threads)'));
end

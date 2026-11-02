function check_t60_runs()
% Checks on the 60-day recovery truths copied from the cluster
% (manifests/recovery_truth_mlp_l2f.csv, 16 cases):
%   1. what each folder holds: time points, last day, step limit, MATLAB, hosts, threads
%   2. repeats: the three 180-day default-grid arrays in outputs_t60_wrong_180d
%      (same settings) against each other
%   3. truncation: a 180-day run cut to 60 days against outputs_t60_base, up to
%      day 59 (point 83 of 84): NetworkTransportModel fills W1 at a point only
%      when it steps on from it, so W1 at the last point of every run is 0
%   4. time step: default grid vs 0.45 vs 0.225 days, at the points they share
% Uses compare_runs.m; see CHECKING.md.

folders = ["outputs_t60_base", "outputs_t60_dt0.45", "outputs_t60_dt0.225", ...
    "outputs_t60_wrong_180d/base", "outputs_t60_wrong_180d/dt0.45", "outputs_t60_wrong_180d/dt0.225"];
M = readtable('manifests/recovery_truth_mlp_l2f.csv', 'TextType', 'string', 'Delimiter', ',');
cases = M.case_id;
file = @(folder, c) fullfile(folder, c + ".mat");

fprintf('== 1. Folders\n');
for d = folders
    pts = []; last = []; dt = []; ml = strings(0); host = strings(0); thr = [];
    for c = cases'
        r = load(file(d, c), 'trange', 'cfg', 'env');
        pts(end+1) = numel(r.trange); last(end+1) = r.trange(end) * 180; %#ok<AGROW>
        dt(end+1) = step_limit(r.cfg); ml(end+1) = r.env.matlab; %#ok<AGROW>
        host(end+1) = r.env.host; thr(end+1) = r.env.max_threads; %#ok<AGROW>
    end
    fprintf('%-30s %d files  points %s  last day %s  dt_max %s (0 = default grid)  threads %s\n   matlab %s\n   hosts %s\n', d, ...
        numel(cases), num2str(unique(pts)), num2str(unique(round(last, 4))), num2str(unique(dt)), ...
        num2str(unique(thr)), strjoin(unique(ml), ', '), strjoin(unique(host), ', '));
end

fprintf('\n== 2. Repeats (180-day default grid, three arrays); N, M, W1, F, EdgeMass all isequal?\n');
reps = folders(4:6);
for c = cases'
    s12 = compare_runs(file(reps(1), c), file(reps(2), c));
    s13 = compare_runs(file(reps(1), c), file(reps(3), c));
    fprintf('%-46s rep1=rep2 %d  rep1=rep3 %d  (max rel %.2g, %.2g)\n', c, all(s12.isequal), ...
        all(s13.isequal), max(s12.max_rel), max(s13.max_rel));
end

fprintf('\n== 3. Truncation: 180-day run cut to day 59 vs outputs_t60_base; all isequal?\n');
for c = cases'
    s = compare_runs(file("outputs_t60_base", c), file(reps(1), c), 59);
    fprintf('%-46s %d  (%d points to day %g, max rel %.2g)\n', c, all(s.isequal), ...
        s.Properties.UserData.n_shared, s.Properties.UserData.last_day, max(s.max_rel));
end

fprintf('\n== 4. Time step, days 0-60, at the points all three grids share\n');
fprintf('   difference to the 0.225-day run; N: NMSE / max rel, F: NMSE\n');
fprintf('%-46s %24s %24s\n', '', 'default grid vs 0.225', '0.45 vs 0.225');
for c = cases'
    sb = compare_runs(file("outputs_t60_dt0.225", c), file("outputs_t60_base", c));
    sh = compare_runs(file("outputs_t60_dt0.225", c), file("outputs_t60_dt0.45", c));
    fprintf('%-46s %9.2g /%6.3f  F %7.2g %9.2g /%6.3f  F %7.2g\n', c, ...
        sb.nmse(1), sb.max_rel(1), sb.nmse(sb.field == "F"), ...
        sh.nmse(1), sh.max_rel(1), sh.nmse(sh.field == "F"));
end
end

function dt = step_limit(cfg)
% SNTM_DT_MAX of the run; 0 = default grid (also for runs made before the setting existed)
dt = 0;
if isfield(cfg, 'dt_max_days') && ~isempty(cfg.dt_max_days), dt = cfg.dt_max_days; end
end

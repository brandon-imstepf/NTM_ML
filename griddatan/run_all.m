% run_all.m
%
% Milestone 11: single entry point that regenerates every artifact in
% saved_models/, figures/, and results/ from the raw input data alone, in
% dependency order, and proves this project reproduces cleanly end to end.
%
% Design note: every s##_*.m script in this project starts with
% `clc; clear;`. Chaining them as plain script calls inside one shared
% MATLAB workspace (e.g. `s00_data_and_scaling_parity; s01_toy_validation; ...`)
% would have each script's own `clear` wipe out THIS script's loop state
% too, since a script-calling-a-script shares one workspace in MATLAB
% (unlike a function call). Spawning a fresh `matlab -batch` subprocess per
% script avoids that entirely, and is a stronger "clean-session" test
% besides -- each script is verified to run correctly starting from
% nothing, not just to not crash when chained after another script already
% populated the workspace. This also matches how every script in this
% project was individually validated throughout development.
%
% Known limitation, stated rather than silently accepted: `system()` here
% blocks with no timeout, so a child script that hangs would hang this
% script forever with no detection. Every script currently completes in
% well under 15 minutes (see the summary table this produces), so this has
% not been worth the added complexity of a real subprocess-timeout
% mechanism -- revisit if a future script's runtime becomes unpredictable.

clc; clear;
here = fileparts(mfilename('fullpath'));

scripts = {
    's00_data_and_scaling_parity'
    's01_toy_validation'
    's02_fast_predictor'
    's03_train_validate_models'
    's04_build_comparison_table'
    's05_compare_nmse_boxplot'
    's06_parameter_recovery_gamma'
    's07_parameter_recovery_joint'
    's08_final_report'
};
% Expected total wall time: a few minutes, dominated by s01_toy_validation
% (griddatan('linear') on the toy problem takes ~35-200s per call,
% confirmed at Milestone 4) -- everything else completes in seconds.

log_dir = fullfile(here, 'results', 'run_all_logs');
if ~isfolder(log_dir); mkdir(log_dir); end

n = numel(scripts);
script_names = strings(n,1);
exit_codes = zeros(n,1);
wall_times = zeros(n,1);

overall_tic = tic;
for i = 1:n
    name = scripts{i};
    fprintf('\n=== [%d/%d] Running %s ===\n', i, n, name);
    log_file = fullfile(log_dir, [name '.log']);
    cmd = sprintf('matlab -batch "%s"', name);

    t0 = tic;
    old_dir = cd(here); % run from this folder so `matlab -batch "name"` resolves the script; each script computes its own absolute paths from mfilename, so this doesn't affect what it reads/writes
    try
        % Two-output system() captures stdout+stderr directly into
        % `cmdout`, rather than relying on shell redirection
        % (`> log 2>&1`) around a nested `matlab -batch` call -- found at
        % Milestone 11's adversarial review to silently fail to produce a
        % log file at all (the nested MATLAB's output bubbled up to this
        % process's own console instead), while leaving exit-code capture
        % unaffected. This two-output form sidesteps that: it works
        % regardless of what the shell-level redirect would have done.
        [status, cmdout] = system(cmd);
        fid = fopen(log_file, 'w');
        fprintf(fid, '%s', cmdout);
        fclose(fid);
    catch ME
        % try/catch here only to guarantee cd(old_dir) below still runs
        % and this project's own working directory isn't left wherever
        % the loop happened to be if system()/fopen ever throws --
        % re-thrown immediately after, not silently swallowed.
        cd(old_dir);
        rethrow(ME);
    end
    cd(old_dir);
    t_elapsed = toc(t0);

    script_names(i) = name;
    exit_codes(i) = status;
    wall_times(i) = t_elapsed;

    if status == 0
        fprintf('  PASSED in %.1f s\n', t_elapsed);
    else
        fprintf('  FAILED (exit code %d) in %.1f s -- see %s\n', status, t_elapsed, log_file);
    end
end
total_time = toc(overall_tic);

results = table(script_names, exit_codes, wall_times, 'VariableNames', {'script','exit_code','wall_time_s'});

fprintf('\n=== run_all.m summary ===\n');
disp(results);
fprintf('Total wall time: %.1f s\n', total_time);

n_failed = sum(exit_codes ~= 0);
if n_failed > 0
    fprintf('\n%d/%d scripts FAILED. See logs in %s for details. NOT a clean reproducibility pass.\n', ...
        n_failed, n, log_dir);
else
    fprintf('\nAll %d scripts completed successfully end-to-end from a clean state -- reproducibility confirmed.\n', n);
end

folder_results = fullfile(here, 'results');
if ~isfolder(folder_results); mkdir(folder_results); end
writetable(results, fullfile(folder_results, 'run_all_summary.csv'));

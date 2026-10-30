function forward_cut_truths(t_end, out_dir)
% MLP and numeric sNTM runs at the recovery truth (g0.0005_l0.02_d80_e20, l2f,
% ECl seed, Hippocampus+PC+RSP) up to t_end days, cut from the 180-day runs in
% outputs_nm/local_tests_10-29/truth, and a figure of tau in each region with
% the numeric run (the truth), the MLP at the true values, and the MLP at the
% gamma1 of the 180-day Nelder-Mead search fitted to the numeric run (testB of
% nm_local_tests.m, its lowest-loss simulation).
% NetworkTransportModel steps forward Euler on the output grid, so the first
% days of a longer run are the shorter run itself; this is checked against the
% 50-day MLP run of nm_failfast_tests.m when that run exists.
%
%   forward_cut_truths(60)   % outputs_nm/forward_60d

here = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(t_end), t_end = 60; end
if nargin < 2 || isempty(out_dir), out_dir = fullfile(here, 'outputs_nm', sprintf('forward_%gd', t_end)); end
if ~isfolder(out_dir), mkdir(out_dir); end
full_dir = fullfile(here, 'outputs_nm', 'local_tests_10-29', 'truth');
case_name = 'ECl__g0.0005_l0.02_d80_e20.mat';

runs = struct();
for m = ["mlp", "numeric"]
    S = load(fullfile(full_dir, "rec_l2f__" + m + "__" + case_name));
    keep = S.trange <= t_end / 180 + 1e-12;
    for f = ["N", "M", "F", "W1", "C", "EdgeMass"]   % fields with time as the last dimension
        if isfield(S, f) && size(S.(f), ndims(S.(f))) == numel(S.trange)
            idx = repmat({':'}, 1, ndims(S.(f)) - 1);
            S.(f) = S.(f)(idx{:}, keep);
        end
    end
    S.trange = S.trange(keep);
    save(fullfile(out_dir, "rec_l2f__" + m + "__" + case_name), '-struct', 'S', '-v7.3');
    runs.(m) = S;
end

% a shorter run is the start of a longer one (forward Euler on the output grid)
short = fullfile(here, 'outputs_nm', 'failfast_50d', 'truth', "rec_l2f__mlp__" + case_name);
if isfile(short)
    R = load(short, 'N');
    n = size(R.N, 2);
    fprintf('50-day MLP run vs the 180-day run cut to 50 days: max |difference| %.3g (max N %.3g)\n', ...
        max(abs(R.N - runs.mlp.N(:, 1:n)), [], 'all'), max(runs.mlp.N, [], 'all'));
end

% MLP at the best gamma1 of the 180-day search fitted to the numeric run
B = load(fullfile(here, 'outputs_nm', 'local_tests_10-29', 'testB', ...
    'nm_l2f__mlp__g0.0005_l0.02_d80_e20__gamma1.mat'), 'trace');
[~, b] = min([B.trace.nmse]);
best_gamma1 = B.trace(b).value;
best_N = B.trace(b).N(:, 1:numel(runs.numeric.trange));
save(fullfile(out_dir, "best_gamma1__mlp__" + case_name), 'best_gamma1', 'best_N');

% tau in each region: numeric, MLP at the truth, MLP at the best gamma1
days = runs.numeric.trange * 180;
names = string(runs.numeric.region_names(:))';
n_reg = numel(names);
n_cols = 6;  n_rows = ceil(n_reg / n_cols);
fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 6*n_cols 5*n_rows + 3]);
tl = tiledlayout(fig, n_rows, n_cols, 'TileSpacing', 'compact', 'Padding', 'compact');
m = runs.numeric.meta;
title(tl, {sprintf('sNTM forward solve to %g days, MLP and numeric', round(days(end), 1)), ...
    sprintf('truth:  \\gamma_1 = %g,  \\lambda_1 = %g,  \\lambda_2 = %g,  \\delta = %g,  \\epsilon = %g,  seed %s', ...
    m.gamma1, m.lambda1, m.lambda2, m.delta, m.epsilon, m.seed)}, 'FontSize', 13);
for r = 1:n_reg
    ax = nexttile(tl); hold(ax, 'on');
    h = [plot(ax, days, runs.numeric.N(r,:), '-k', 'LineWidth', 2), ...
         plot(ax, days, runs.mlp.N(r,:), '--', 'Color', [0.15 0.3 0.6], 'LineWidth', 1.5), ...
         plot(ax, days, best_N(r,:), ':', 'Color', [0.8 0.1 0.1], 'LineWidth', 2)];
    set(ax, 'Box', 'off', 'TickDir', 'out', 'FontSize', 8);
    grid(ax, 'on');
    title(ax, wrapName(names(r)), 'FontWeight', 'normal', 'FontSize', 9, 'Interpreter', 'none');
    if r > n_reg - n_cols, xlabel(ax, 'days'); end
    if mod(r - 1, n_cols) == 0, ylabel(ax, '\tau'); end
    hold(ax, 'off');
end
lg = legend(h, {'numeric at the truth', 'MLP at the truth', ...
    sprintf('MLP at the best fit to numeric (\\gamma_1 = %.3g, 180-day search)', best_gamma1)}, ...
    'Orientation', 'horizontal', 'Box', 'off', 'FontSize', 11);
lg.Layout.Tile = 'south';
exportgraphics(fig, fullfile(out_dir, sprintf('forward_%gd_regions.png', t_end)), 'Resolution', 150);
close(fig);
fprintf('saved to %s (%d time points, last at %.1f days)\n', out_dir, numel(days), days(end));
end

function lines = wrapName(name)
% region names on two lines, split at the first comma (as in plotNelderMeadRecovery.m),
% with an empty second line otherwise so every title has the same height
    lines = [name, ""];
    k = strfind(name, ", ");
    if ~isempty(k)
        lines = [extractBefore(name, k(1) + 1), extractAfter(name, k(1) + 1)];
    end
end

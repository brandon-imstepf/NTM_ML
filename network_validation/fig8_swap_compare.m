function T = fig8_swap_compare()
% Total soluble tau sum_i N_i of each Fig. 8 run divided by the numeric run, at
% selected days, and the per-month slope of sNTM against numeric regional N, for
% the swap test (FIG8_OFFSET.md, step 3). Runs that do not exist yet are skipped.
% Writes results/fig8_swap_compare.csv.
%
%   T = fig8_swap_compare();

here = fileparts(mfilename('fullpath'));
f = @(m) fullfile(here, 'outputs_fig8', "fig8__" + m + "__ECl__g5e-05_l0.01_d10_e50.mat");
num = load(f("numeric"), 'N', 'trange');
day = num.trange * 180;
runs = ["mlp", "gp", "swapA_numflux_mlpw1", "swapB_mlpflux_numw1"];
days = [0.45 2.25 10 30 60 150];
rows = {};
for r = runs
    if ~isfile(f(r)), fprintf('%s: not available yet\n', r); continue; end
    s = load(f(r), 'N');
    ratio = arrayfun(@(t) sum(s.N(:, find(day >= t, 1))) / sum(num.N(:, find(day >= t, 1))), days);
    slope = arrayfun(@(mo) subsref(polyfit(num.N(:, find(day >= 30 * mo, 1)), s.N(:, find(day >= 30 * mo, 1)), 1), ...
        struct('type', '()', 'subs', {{1}})), 1:5);
    rows(end+1, :) = [{r}, num2cell(ratio), num2cell(slope)]; %#ok<AGROW>
end
T = cell2table(rows, 'VariableNames', [{'run'}, compose('ratio_day%g', days), compose('slope_month%d', 1:5)]);
writetable(T, fullfile(here, 'results', 'fig8_swap_compare.csv'));
disp(T)
end

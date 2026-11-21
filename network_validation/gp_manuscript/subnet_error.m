function [T, M] = subnet_error(src_fmt, tag)
% Network-scale test (R2.M1): error of each sNTM against the numeric NTM for the
% runs of make_scale_manifest.py (Primary Parameter Set, seed BLA, 180 days).
% (11/18) The first cluster runs (outputs_subnet/) were made with a simulator copy
% without the 'Amygdala+CTXsp' case, so connectome_subset fell through to the
% whole connectome (426 regions, 9,934 directed edges) in both adjacency_cached.m
% and NetworkTransportModel.m: these are whole-brain runs.
%   T: per method, NMSE (Eq. 7, nmse_eq7.m) over all regions and time points,
%      finite, runtime (s) and speedup over the numeric run.
%   M: per method and month 1-5, over the regions: NMSE, slope b of the
%      least-squares fit through the origin and R^2 = squared Pearson correlation
%      (as in main text Fig. 10, plot_fig8_gp.m).
% Writes results/subnet_error<tag>.csv and results/subnet_error_by_month<tag>.csv.
%
%   [T, M] = subnet_error();     % outputs_subnet/, tag "_wholebrain"

here = fileparts(mfilename('fullpath'));
nv = fileparts(here);
addpath(nv);
if nargin < 1 || isempty(src_fmt), src_fmt = 'outputs_subnet/subnet__%s__BLA__g0.005_l0.02_d20_e80__amyg.mat'; end
if nargin < 2, tag = "_wholebrain"; end
methods = ["lr", "sr", "interp", "gp", "mlp"];   % paper order
num = load(fullfile(nv, sprintf(src_fmt, 'numeric')), 'N', 'C', 'trange', 'runtime_s');
day = num.trange * 180;
hm = arrayfun(@(m) find(abs(day - 30 * m) == min(abs(day - 30 * m)), 1), 1:5);
fprintf('%d regions, %d directed edges, %d time points; numeric runtime %.0f s\n', ...
    size(num.N, 1), nnz(num.C), numel(day), num.runtime_s);
rows = {};  mrows = {};
for m = methods
    s = load(fullfile(nv, sprintf(src_fmt, m)), 'N', 'trange', 'runtime_s');
    assert(isequal(size(s.N), size(num.N)) && max(abs(s.trange - num.trange)) < 1e-12, 'time grid differs: %s', m);
    fin = all(isfinite(s.N(:)));
    e = NaN;  if fin, e = nmse_eq7(s.N, num.N); end
    rows(end+1, :) = {m, e, fin, s.runtime_s, num.runtime_s / s.runtime_s}; %#ok<AGROW>
    for k = 1:5
        x = num.N(:, hm(k));  y = s.N(:, hm(k));
        if all(isfinite(y))
            R = corrcoef(x, y);  b = sum(x .* y) / sum(x .^ 2);
            mrows(end+1, :) = {m, k, day(hm(k)), nmse_eq7(y, x), b, R(1, 2)^2}; %#ok<AGROW>
        else
            mrows(end+1, :) = {m, k, day(hm(k)), NaN, NaN, NaN}; %#ok<AGROW>
        end
    end
end
T = cell2table(rows, 'VariableNames', {'method', 'nmse', 'finite', 'runtime_s', 'speedup'});
M = cell2table(mrows, 'VariableNames', {'method', 'month', 'day', 'nmse', 'slope_origin', 'r2_corr'});
writetable(T, fullfile(nv, 'results', "subnet_error" + tag + ".csv"));
writetable(M, fullfile(nv, 'results', "subnet_error_by_month" + tag + ".csv"));
disp(T);  disp(M);
end

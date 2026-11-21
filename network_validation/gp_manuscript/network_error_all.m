function T = network_error_all()
% Network error of every sNTM against the numeric NTM in all network settings
% of R1.m1 (fig8_cases.m): "interior" (8 settings x 3 seeds), "edge" (8 settings
% on the corners of the training range) and "ood" (8 settings outside it), ECl
% seed for edge and ood; 180 days. NMSE (Eq. 7, nmse_eq7.m) over all regions and
% time points, as Supplementary Table S2 (network_checks_l2f.m); R^2 = 1 - NMSE.
% The nearest-neighbor sNTM (interp; log gamma/lambda, 10^6 lookup table) runs
% are those of manifests/consistent_l2f_interp.csv and
% consistent_edge_ood_interp.csv in outputs_consistent/.
% A run that produced non-finite concentrations has NMSE NaN, finite 0.
% Writes results/network_error_all.csv and results/network_error_by_month.csv (NMSE
% and least-squares R^2 over the 30 regions at months 1-5) and prints median [min, max] per group
% and method.
%
%   T = network_error_all();

here = fileparts(mfilename('fullpath'));
nv = fileparts(here);
addpath(nv);
methods = ["lr", "sr", "interp", "gp", "mlp"];   % paper order (11/17)
cases = fig8_cases(["interior", "edge", "ood"]);
rows = {};  mrows = {};
for c = cases
    num_path = c.num_file;  if isempty(num_path), num_path = sprintf(c.src_fmt, 'numeric'); end
    num = load(fullfile(nv, num_path), 'N', 'trange');
    day = num.trange * 180;
    hm = arrayfun(@(m) find(abs(day - 30 * m) == min(abs(day - 30 * m)), 1), 1:5);
    for m = methods
        if m == "interp" && c.group ~= "interior"
            f = fullfile(nv, 'outputs_consistent', c.group + "__" + c.setting + "__" + c.seed + "__interp.mat");
        else
            f = fullfile(nv, sprintf(c.src_fmt, m));
        end
        if ~isfile(f), fprintf('missing %s\n', f); continue, end
        s = load(f, 'N', 'trange');
        assert(isequal(size(s.N), size(num.N)) && max(abs(s.trange - num.trange)) < 1e-12, 'time grid differs: %s', f);
        fin = all(isfinite(s.N(:)));
        e = NaN;  if fin, e = nmse_eq7(s.N, num.N); end
        v = sscanf(char(c.setting), 'g%f_l%f_d%f_e%f')';
        rows(end+1, :) = {c.group, c.seed, c.setting, v(1), v(2), v(3), v(4), m, e, 1 - e, fin}; %#ok<AGROW>
        % by month, over the 30 regions (as in the network figure of the main text):
        % NMSE (Eq. 7) and R^2 of the least-squares fit; NaN once the run is non-finite.
        % (11/17) also the fit through the origin of the network figure (plot_fig8_gp.m):
        % slope b = sum(x y) / sum(x^2) and its R^2 = 1 - sum (y - b x)^2 / sum (y - mean y)^2
        for k = 1:5
            x = num.N(:, hm(k));  y = s.N(:, hm(k));
            if all(isfinite(y))
                R = corrcoef(x, y);  q = polyfit(x, y, 1);
                b = sum(x .* y) / sum(x .^ 2);  r2o = 1 - sum((y - b * x) .^ 2) / sum((y - mean(y)) .^ 2);
                mrows(end+1, :) = {c.group, c.seed, c.setting, m, k, day(hm(k)), nmse_eq7(y, x), R(1, 2)^2, ...
                    q(1), q(2) * 1e3, b, r2o, 1}; %#ok<AGROW>
            else
                mrows(end+1, :) = {c.group, c.seed, c.setting, m, k, day(hm(k)), NaN, NaN, NaN, NaN, NaN, NaN, 0}; %#ok<AGROW>
            end
        end
    end
end
M = cell2table(mrows, 'VariableNames', {'group', 'seed', 'setting', 'method', 'month', 'day', 'nmse', 'r2_fit', ...
    'slope', 'intercept_nM', 'slope_origin', 'r2_origin', 'finite'});
writetable(M, fullfile(nv, 'results', 'network_error_by_month.csv'));
T = cell2table(rows, 'VariableNames', {'group', 'seed', 'setting', 'gamma', 'lambda', 'delta', 'epsilon', ...
    'method', 'nmse', 'r2', 'finite'});
writetable(T, fullfile(nv, 'results', 'network_error_all.csv'));
for g = ["interior", "edge", "ood"]
    fprintf('== %s\n', g);
    for m = methods
        a = T(T.group == g & T.method == m, :);
        x = a.nmse(a.finite == 1);
        fprintf('%-7s n %2d, non-finite %d: NMSE median %.3g [%.3g, %.3g]; NMSE < 0.01 in %d, < 0.1 in %d\n', m, ...
            height(a), nnz(~a.finite), median(x), min(x), max(x), nnz(x < 0.01), nnz(x < 0.1));
    end
end
end

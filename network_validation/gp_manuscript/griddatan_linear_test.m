function T = griddatan_linear_test(sizes, n_eval)
% griddatan 'linear' (barycentric interpolation in the Delaunay simplex that
% contains the query; delaunayn + tsearchn, what griddatan does, triangulated
% once) against 'nearest', on the same scaled inputs as the nearest-neighbor sNTM
% (log10 gamma1, log10 lambda1, then all six z-scored). Queries outside the
% convex hull of the training inputs take the nearest value (fraction reported).
% Training rows: randperm(10^5, n) of data_bias_e5_1.csv after rng(42).
% Evaluated on the 10^3-sample validation set (bias_val_e3.csv) and on n_eval
% samples of the held-out set (data_bias_e5_11.csv, rng(1)); NMSE = Eq. 7.
% Records triangulation size and time, and query time.
% Writes results/griddatan_linear_test.csv.
%
%   T = griddatan_linear_test([1000 2000 5000], 2000);

if nargin < 1 || isempty(sizes), sizes = [1000 2000 5000]; end
if nargin < 2 || isempty(n_eval), n_eval = 2000; end
here = fileparts(mfilename('fullpath'));
nv = fileparts(here);
addpath(nv);
cfg = sntm_config();
super = readmatrix(cfg.train_superset_file);
val = readmatrix(fullfile(cfg.ntm_ml_root, 'data', 'data_baseline_models', 'e3_data', 'bias_val_e3.csv'));
held = readmatrix(fullfile(cfg.ntm_ml_root, 'data', 'data_neural_network_model', 'training_data_e6', 'data_bias_e5_11.csv'));
rng(1);  held = held(randperm(size(held, 1), n_eval), :);
lg = @(X) [log10(X(:, 1:2)), X(:, 3:6)];
rows = {};
for n = sizes
    rng(42);  D = super(randperm(size(super, 1), n), :);
    mu = mean(lg(D(:, 1:6)));  sd = std(lg(D(:, 1:6)));
    tf = @(X) (lg(X) - mu) ./ sd;
    Xt = tf(D(:, 1:6));  Y = D(:, 7:8);
    ns = createns(Xt, 'NSMethod', 'kdtree');
    t0 = tic;  tri = delaunayn(Xt);  t_tri = toc(t0);
    fprintf('n %d: %d simplices, delaunayn %.0f s\n', n, size(tri, 1), t_tri);
    for set = ["val", "heldout"]
        if set == "val", E = val; else, E = held; end
        Xq = tf(E(:, 1:6));
        Pn = Y(knnsearch(ns, Xq), :);
        t0 = tic;  [k, w] = tsearchn(Xt, tri, Xq);  t_q = toc(t0) / size(Xq, 1);
        Pl = Pn;  in = ~isnan(k);
        for c = 1:2
            y = Y(:, c);
            Pl(in, c) = sum(y(tri(k(in), :)) .* w(in, :), 2);
        end
        rows(end+1, :) = {n, set, size(tri, 1), t_tri, t_q, mean(~in), ...
            nmse_eq7(Pl(:, 1), E(:, 7)), nmse_eq7(Pl(:, 2), E(:, 8)), ...
            nmse_eq7(Pn(:, 1), E(:, 7)), nmse_eq7(Pn(:, 2), E(:, 8)), ...
            nmse_eq7(Pl(in, 1), E(in, 7)), nmse_eq7(Pl(in, 2), E(in, 8))}; %#ok<AGROW>
        fprintf('  %-7s outside hull %.1f%%, query %.3g s: linear NMSE %.3g / %.3g (inside hull %.3g / %.3g), nearest %.3g / %.3g\n', ...
            set, 100 * mean(~in), t_q, rows{end, 7}, rows{end, 8}, rows{end, 11}, rows{end, 12}, rows{end, 9}, rows{end, 10});
    end
end
T = cell2table(rows, 'VariableNames', {'n_train', 'eval_set', 'n_simplices', 'triangulation_s', 'query_s', ...
    'frac_outside_hull', 'linear_flux_nmse', 'linear_emc_nmse', 'nearest_flux_nmse', 'nearest_emc_nmse', ...
    'linear_inside_flux_nmse', 'linear_inside_emc_nmse'});
writetable(T, fullfile(nv, 'results', 'griddatan_linear_test.csv'));
end

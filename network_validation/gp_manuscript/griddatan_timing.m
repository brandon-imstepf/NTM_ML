function T = griddatan_timing(n_lin, n_near, n_query)
% Computational cost of griddatan 'linear' against 'nearest' on the same
% training samples and the same queries, one CPU thread: set-up (Delaunay
% triangulation, delaunayn; k-d tree, createns) and time per query (tsearchn plus
% the barycentric interpolation; knnsearch), against the number of training
% samples n. Inputs scaled as the nearest-neighbor sNTM (log10 gamma1, log10
% lambda1, then all z-scored). Training rows: randperm(10^5, n) of
% data_bias_e5_1.csv after rng(42), and the 10^6 MLP training samples
% (interp_training_data.m) for n > 10^5. Queries: n_query held-out samples
% (data_bias_e5_11.csv, rng(1)). 'linear' is limited to the sizes that can be
% triangulated (the number of simplices grows rapidly in six dimensions);
% plot_griddatan_timing.py extrapolates it.
% Writes results/griddatan_timing.csv.
%
%   T = griddatan_timing([500 1000 2000 3000 5000 7000], [500 1e3 2e3 5e3 1e4 1e5 1e6], 200);

if nargin < 1 || isempty(n_lin), n_lin = [500 1000 2000 3000 5000 7000]; end
if nargin < 2 || isempty(n_near), n_near = [500 1000 2000 3000 5000 7000 1e4 3e4 1e5 3e5 1e6]; end
if nargin < 3 || isempty(n_query), n_query = 200; end
maxNumCompThreads(1);
here = fileparts(mfilename('fullpath'));
nv = fileparts(here);
addpath(nv);
cfg = sntm_config();
super = readmatrix(cfg.train_superset_file);
full = interp_training_data(cfg);
held = readmatrix(fullfile(cfg.ntm_ml_root, 'data', 'data_neural_network_model', 'training_data_e6', 'data_bias_e5_11.csv'));
rng(1);  Q = held(randperm(size(held, 1), n_query), 1:6);
lg = @(X) [log10(X(:, 1:2)), X(:, 3:6)];
rows = {};
for n = unique([n_lin n_near])
    rng(42);
    if n <= 1e5, D = super(randperm(size(super, 1), n), :); else, D = full(randperm(size(full, 1), n), :); end
    mu = mean(lg(D(:, 1:6)));  sd = std(lg(D(:, 1:6)));
    Xt = (lg(D(:, 1:6)) - mu) ./ sd;  Xq = (lg(Q) - mu) ./ sd;  Y = D(:, 7);
    if ismember(n, n_near)
        t0 = tic;  ns = createns(Xt, 'NSMethod', 'kdtree');  t_set = toc(t0);
        reps = max(1, round(2e4 / n_query));        % repeat the queries so the timing is measurable
        t0 = tic;  for r = 1:reps, k = knnsearch(ns, Xq); end;  t_q = toc(t0) / (reps * n_query); %#ok<NASGU>
        rows(end+1, :) = {"nearest", n, t_set, t_q, NaN}; %#ok<AGROW>
        fprintf('nearest n %7d: set-up %.3g s, query %.3g s\n', n, t_set, t_q);
    end
    if ismember(n, n_lin)
        t0 = tic;  tri = delaunayn(Xt);  t_set = toc(t0);
        t0 = tic;  [k, w] = tsearchn(Xt, tri, Xq);
        in = ~isnan(k);  p = nan(n_query, 1);  p(in) = sum(Y(tri(k(in), :)) .* w(in, :), 2); %#ok<NASGU>
        t_q = toc(t0) / n_query;
        rows(end+1, :) = {"linear", n, t_set, t_q, size(tri, 1)}; %#ok<AGROW>
        fprintf('linear  n %7d: %d simplices, set-up %.3g s, query %.3g s\n', n, size(tri, 1), t_set, t_q);
        clear tri
    end
    T = cell2table(rows, 'VariableNames', {'method', 'n_train', 'setup_s', 'query_s', 'n_simplices'});
    writetable(T, fullfile(nv, 'results', 'griddatan_timing.csv'));
end
end

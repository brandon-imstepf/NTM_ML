function T = supp_gp_griddatan(steps)
% Supplement (advisor meeting 11/10): the GP trained on the same 10^3-sample
% dataset as linear and symbolic regression, a GP training-size curve, and
% scattered-data interpolation (griddatan) as a lookup-table baseline, all in
% the format of Table 1 of the paper.
%
% Models
%   gp_lrsr_1000   exact GP, linear basis, ARD Matern 3/2, fminunc (the configuration
%                  selected for the paper's GP), trained on bias_train_e3.csv, the
%                  training set of LR and SR
%   gp_<n>         the same configuration trained on n = 500, 1000, 2000, 3000, 4000
%                  rows of data_bias_e5_1.csv (rows randperm(1e5, n) after
%                  rng(42), the sampling of gaussian_process.m); n = 5000 is the
%                  paper's GP (models/, not refitted)
%   griddatan_nearest_raw_<n>     griddatan 'nearest' on the raw inputs, as in
%                  griddatan/s03 (the 1.83 of the response, there with n = 3000)
%   griddatan_nearest_scaled_<n>  'nearest' on scaled inputs: log10 gamma1,
%                  log10 lambda1, then each input z-scored with the training mean
%                  and s.d.
%   griddatan_linear_scaled_1000  'linear' on the same scaled inputs; queries
%                  outside the convex hull of the training inputs (NaN) take the
%                  'nearest' value; the fraction outside is reported. Evaluated on
%                  the test set only: point location (tsearchn) takes about 0.16 s
%                  per query in 6-D with 1000 training points (6e5 simplices), and
%                  the 5000-point triangulation has 4.2e6 simplices
%   griddatan for n = 1000 (bias_train_e3.csv) and n = 5000 (the GP's 5,000 rows;
%   'nearest' only)
% Targets: GP on z-scored targets (training mean and s.d.), as gaussian_process.m;
% griddatan interpolates the raw targets (equivalent for both methods).
% NMSE = Eq. 7 (nmse_eq7.m); mean +- s.d. over 10^3 bootstrap resamples of 10^3
% samples with replacement (rng(1) before each set), from
%   train    the model's own training rows
%   test     bias_val_e3.csv (10^3; the test set of Table 1 for LR and SR, and the
%            set on which the paper's GP configuration was selected)
%   heldout  data_bias_e5_11.csv (10^5; never used for fitting or selection)
% plus the NMSE over the whole set (_all).
% steps: "fit" (GP fits, cached in gp_manuscript/models_size/), "eval", or both.
% Writes results/supp_gp_griddatan.csv and results/supp_gp_size_curve.csv.
%
%   T = supp_gp_griddatan();

if nargin < 1, steps = ["fit", "eval"]; end
steps = string(steps);
here = fileparts(mfilename('fullpath'));
nv = fileparts(here);
addpath(nv);
cfg = sntm_config();
root = cfg.ntm_ml_root;
model_dir = fullfile(here, 'models_size');
if ~isfolder(model_dir), mkdir(model_dir); end

super = readmatrix(fullfile(root, 'data', 'data_neural_network_model', 'training_data_e6', 'data_bias_e5_1.csv'));
lrsr = readmatrix(fullfile(root, 'data', 'data_baseline_models', 'e3_data', 'bias_train_e3.csv'));
test = readmatrix(fullfile(root, 'data', 'data_baseline_models', 'e3_data', 'bias_val_e3.csv'));
heldout = readmatrix(fullfile(root, 'data', 'data_neural_network_model', 'training_data_e6', 'data_bias_e5_11.csv'));
G5 = load(cfg.gp_model_file);
sizes = [500 1000 2000 3000 4000];
super_rows = @(n) super(randperm_seed(42, 1e5, n), :);
assert(isequal(super_rows(G5.dataset_size), super_rows(5000)) || G5.dataset_size == 5000);

%% GP fits (cached)
train_sets = struct('gp_lrsr_1000', lrsr);
for n = sizes, train_sets.(sprintf('gp_%d', n)) = super_rows(n); end
if any(steps == "fit")
    for f = string(fieldnames(train_sets))'
        out = fullfile(model_dir, f + ".mat");
        if isfile(out), continue; end
        D = train_sets.(f);
        t0 = tic;
        M = fit_gp(D);
        M.fit_time_s = toc(t0);
        save(out, '-struct', 'M');
        fprintf('%s: fitted on %d rows in %.0f s\n', f, size(D, 1), M.fit_time_s);
    end
end
if ~any(steps == "eval"), T = []; return; end

%% Predictors
pred = struct();  train_of = struct();
for f = string(fieldnames(train_sets))'
    M = load(fullfile(model_dir, f + ".mat"));
    pred.(f) = @(X) gp_predict(M, X);  train_of.(f) = train_sets.(f);
end
pred.gp_5000 = @(X) [predict(G5.stzd_flux_gprModel, X) * G5.flux_std + G5.flux_mean, ...
                     predict(G5.stzd_emc_gprModel, X) * G5.emc_std + G5.emc_mean];
train_of.gp_5000 = super_rows(5000);
gd_sets = struct('n1000', lrsr, 'n5000', super_rows(5000));
inside = struct();
for g = ["n1000", "n5000"]
    D = gd_sets.(g);  n = extractAfter(g, "n");
    pred.("griddatan_nearest_raw_" + n) = make_gd(D, 'nearest', false);
    pred.("griddatan_nearest_scaled_" + n) = make_gd(D, 'nearest', true);
    kinds = ["nearest_raw", "nearest_scaled"];
    if n == "1000"
        [pred.("griddatan_linear_scaled_" + n), inside.("griddatan_linear_scaled_" + n)] = make_gd(D, 'linear', true);
        kinds(end+1) = "linear_scaled"; %#ok<AGROW>
    end
    for kind = kinds
        train_of.("griddatan_" + kind + "_" + n) = D;
    end
end

%% NMSE
names = string(fieldnames(pred))';
rows = {};
for f = names
    P = struct();
    sets = struct('train', train_of.(f), 'test', test, 'heldout', heldout);
    out_frac = NaN;
    for s = ["train", "test", "heldout"]
        D = sets.(s);
        if contains(f, "linear") && s ~= "test"   % test set only (see above)
            P.(s) = NaN(1, 6);
            continue
        end
        p = chunked(pred.(f), D(:, 1:6));
        if contains(f, "linear")
            out_frac = mean(~inside.(f)(D(:, 1:6)));
        end
        [mf, sf] = boot(p(:, 1), D(:, 7));  [me, se] = boot(p(:, 2), D(:, 8));
        P.(s) = [mf, sf, me, se, nmse_eq7(p(:, 1), D(:, 7)), nmse_eq7(p(:, 2), D(:, 8))];
    end
    rows(end+1, :) = {f, size(train_of.(f), 1), P.train, P.test, P.heldout, out_frac}; %#ok<AGROW>
    fprintf('%-30s n=%5d  test flux %.3g emc %.3g | heldout flux %.3g emc %.3g\n', f, size(train_of.(f), 1), ...
        P.test(1), P.test(3), P.heldout(1), P.heldout(3));
end
T = table();
T.model = string(rows(:, 1));
T.n_train = cell2mat(rows(:, 2));
lab = ["flux_nmse_mean", "flux_nmse_sd", "emc_nmse_mean", "emc_nmse_sd", "flux_nmse_all", "emc_nmse_all"];
for s = ["train", "test", "heldout"]
    k = find(["train", "test", "heldout"] == s);
    V = cell2mat(rows(:, 2 + k));
    for j = 1:numel(lab), T.(s + "_" + lab(j)) = V(:, j); end
end
T.test_frac_outside_hull = cell2mat(rows(:, 6));
writetable(T, fullfile(nv, 'results', 'supp_gp_griddatan.csv'));

% training-size curve: fit time and held-out NMSE of the GP by training size
C = T(startsWith(T.model, "gp_") & T.model ~= "gp_lrsr_1000", {'model', 'n_train', 'heldout_flux_nmse_mean', ...
    'heldout_flux_nmse_sd', 'heldout_emc_nmse_mean', 'heldout_emc_nmse_sd'});
C = sortrows(C, 'n_train');
writetable(C, fullfile(nv, 'results', 'supp_gp_size_curve.csv'));
disp(T(:, {'model', 'n_train', 'test_flux_nmse_mean', 'test_emc_nmse_mean', 'heldout_flux_nmse_mean', 'heldout_emc_nmse_mean', 'test_frac_outside_hull'}))
end

function idx = randperm_seed(seed, n, k)
rng(seed);  idx = randperm(n, k);
end

function M = fit_gp(D)
% the configuration of the paper's GP, fitted as in gaussian_process.m
X = D(:, 1:6);
M.flux_mean = mean(D(:, 7));  M.flux_std = std(D(:, 7));
M.emc_mean = mean(D(:, 8));   M.emc_std = std(D(:, 8));
opts = {'Standardize', true, 'FitMethod', 'exact', 'BasisFunction', 'linear', ...
    'KernelFunction', 'ardmatern32', 'Optimizer', 'fminunc'};
M.flux = fitrgp(X, (D(:, 7) - M.flux_mean) / M.flux_std, opts{:});
M.emc = fitrgp(X, (D(:, 8) - M.emc_mean) / M.emc_std, opts{:});
M.n_train = size(D, 1);
end

function P = gp_predict(M, X)
P = [predict(M.flux, X) * M.flux_std + M.flux_mean, predict(M.emc, X) * M.emc_std + M.emc_mean];
end

function Z = scale_inputs(X, ref)
X(:, 1:2) = log10(X(:, 1:2));  ref(:, 1:2) = log10(ref(:, 1:2));
Z = (X - mean(ref)) ./ std(ref);
end

function [f, inside] = make_gd(D, method, scaled)
% griddatan as a reusable predictor: 'nearest' = the nearest training input
% (k-d tree; identical to griddatan 'nearest', griddatan/README.md, milestone 5);
% 'linear' = barycentric interpolation in the Delaunay simplex containing the
% query (delaunayn + tsearchn, what griddatan 'linear' does), triangulated once;
% outside the convex hull, the nearest value.
Xt = D(:, 1:6);
tf = @(X) X;
if scaled, tf = @(X) scale_inputs(X, D(:, 1:6)); Xt = tf(Xt); end
Y = D(:, 7:8);
ns = createns(Xt, 'NSMethod', 'kdtree');
nearest = @(Xq) Y(knnsearch(ns, Xq), :);
if method == "nearest"
    f = @(X) nearest(tf(X));
    inside = @(X) true(size(X, 1), 1);
    return
end
t0 = tic;
tri = delaunayn(Xt);
fprintf('delaunayn: %d points, %d simplices, %.0f s\n', size(Xt, 1), size(tri, 1), toc(t0));
f = @(X) linear_pred(Xt, tri, Y, tf(X), nearest);
inside = @(X) ~isnan(tsearchn(Xt, tri, tf(X)));
end

function P = linear_pred(Xt, tri, Y, Xq, nearest)
[k, w] = tsearchn(Xt, tri, Xq);
P = nearest(Xq);
in = ~isnan(k);
for c = 1:2
    y = Y(:, c);
    P(in, c) = sum(y(tri(k(in), :)) .* w(in, :), 2);
end
end

function y = chunked(f, X)
y = zeros(size(X, 1), 2);
for a = 1:20000:size(X, 1)
    b = min(a + 19999, size(X, 1));
    y(a:b, :) = f(X(a:b, :));
end
end

function [mu, sd] = boot(p, y)
rng(1);
v = zeros(1000, 1);
for i = 1:1000
    j = randi(numel(y), 1000, 1);
    v(i) = nmse_eq7(p(j), y(j));
end
mu = mean(v);  sd = std(v);
end

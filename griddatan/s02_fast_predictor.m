% s02_fast_predictor.m
%
% Milestone 5: a reusable "fit once, query many times" nearest-neighbor
% predictor that exactly reproduces griddatan(...,'nearest') but without
% rebuilding anything per query. This is the engineering fix the
% parameter-recovery stage requires to be tractable at all: fminbnd and the
% joint-recovery grid sweep both make hundreds to thousands of repeated
% queries per dataset row, and re-running griddatan from scratch for each
% one (~35-45 s per call at toy scale, Milestone 4) would make recovery
% take hours to days.
%
% Design note (informed by Milestone 4): griddatan('nearest',...) does not
% use a Delaunay/Qhull triangulation internally -- its timing showed none of
% the interior/exterior asymmetry that 'linear' (genuinely Qhull-based)
% displayed. So the reusable predictor here is a k-d tree
% (MATLAB's createns/knnsearch), not the delaunayn+tsearchn wrapper
% originally planned in README.md's file-layout table -- simpler and
% cheaper, and the correct match for the method Milestone 4 actually chose.

clc; clear;
here = fileparts(mfilename('fullpath'));
gp_dir = fullfile(here, '..', 'gaussian_process-interpolant');

%% Parity test 1 -- toy problem (same setup as Milestone 4, s01)
rng(7);
f = @(X) sin(X(:,1)) + 0.5*cos(2*X(:,2)) + X(:,3).*X(:,4) - 0.3*X(:,5).^2 + log1p(X(:,6).^2);
n_train = 2000; n_test = 500;
X_train_toy = rand(n_train, 6) * 2*pi;
y_train_toy = f(X_train_toy);
X_test_interior = rand(n_test, 6) * 2*pi;
X_test_exterior = 2*pi + rand(n_test, 6) * 2*pi;

model_toy = fit_nearest_predictor(X_train_toy, y_train_toy);

y_griddatan_int = griddatan(X_train_toy, y_train_toy, X_test_interior, 'nearest');
y_fast_int = predict_nearest(model_toy, X_test_interior);
d_int = max(abs(y_griddatan_int - y_fast_int));
fprintf('[Toy, interior] griddatan vs. fast predictor, max abs diff: %g\n', d_int);
assert(d_int == 0, 'Fast predictor does not exactly reproduce griddatan on interior toy queries.');

y_griddatan_ext = griddatan(X_train_toy, y_train_toy, X_test_exterior, 'nearest');
y_fast_ext = predict_nearest(model_toy, X_test_exterior);
d_ext = max(abs(y_griddatan_ext - y_fast_ext));
fprintf('[Toy, exterior] griddatan vs. fast predictor, max abs diff: %g\n', d_ext);
assert(d_ext == 0, 'Fast predictor does not exactly reproduce griddatan on exterior toy queries.');

%% Parity test 2 -- real data, using the Milestone 2-verified 3000-row corpus
% Same reconstruction as s00_data_and_scaling_parity.m (rng(42);
% randperm(100000,3000)), proven there to match the real GP training corpus
% exactly. This is the more relevant test: real feature scales, real
% duplicate/near-duplicate rows, real validation set -- not just a clean
% synthetic function.
tr_filepath = fullfile(here, '..', 'data', 'data_neural_network_model', 'training_data_e6', 'data_bias_e5_1.csv');
va_filepath = fullfile(here, '..', 'data', 'data_baseline_models', 'e3_data', 'bias_val_e3.csv');
training_superset = readmatrix(tr_filepath);
validation_dataset = readmatrix(va_filepath);

rng(42); idx_3000 = randperm(100000, 3000);
subset_3000 = training_superset(idx_3000, :);
X_train_real = subset_3000(:, 1:6);
flux_train_real = subset_3000(:, 7);
val_features_real = validation_dataset(:, 1:6);

tic;
model_real = fit_nearest_predictor(X_train_real, flux_train_real);
t_fit_real = toc;

tic;
y_griddatan_real = griddatan(X_train_real, flux_train_real, val_features_real, 'nearest');
t_griddatan_real = toc;

tic;
y_fast_real = predict_nearest(model_real, val_features_real);
t_fast_real = toc;

d_real = max(abs(y_griddatan_real - y_fast_real));
fprintf('[Real data] griddatan vs. fast predictor (flux, full validation set), max abs diff: %g\n', d_real);
assert(d_real == 0, 'Fast predictor does not exactly reproduce griddatan on real validation data.');

fprintf('\nOne-time fit cost (3000 training rows): %.3f s\n', t_fit_real);
fprintf('griddatan(''nearest'',...) cost for the WHOLE %d-row validation batch: %.3f s\n', ...
    size(val_features_real,1), t_griddatan_real);
fprintf('Fast predictor cost for the SAME whole batch (after one-time fit): %.3f s\n', t_fast_real);

%% Repeated small-query benchmark -- the actual recovery-loop workload shape
% The recovery loop (fminbnd, or the joint grid sweep) doesn't query in one
% big batch -- it queries ONE point at a time, hundreds to thousands of
% times per dataset row. That's the workload this benchmark reproduces.
n_reps = 500;
rng(1);
single_queries = val_features_real(randi(size(val_features_real,1), n_reps, 1), :);

tic;
for i = 1:n_reps
    predict_nearest(model_real, single_queries(i,:));
end
t_fast_many = toc;

fprintf('\n%d repeated single-point queries via the fast predictor: %.3f s total (%.4f s/query)\n', ...
    n_reps, t_fast_many, t_fast_many/n_reps);
fprintf(['Calling griddatan fresh for each of those %d single-point queries would cost ' ...
    'roughly %d x the %.3f s one-shot-call time measured above (~%.0f s) -- the ' ...
    'reusable predictor turns an intractable recovery loop into a practical one.\n'], ...
    n_reps, n_reps, t_griddatan_real, n_reps * t_griddatan_real);

fprintf('\nMilestone 5 PASSED: fast predictor is exact (diff=0) against raw griddatan, both on the toy problem and on real training/validation data.\n');

%% ---- Local functions ----

function model = fit_nearest_predictor(X_train, y_train)
    % One-time cost: build a k-d tree over the training features. Reused by
    % every subsequent predict_nearest() call -- no re-triangulation, no
    % re-scan of the training set.
    %
    % Known gap vs. griddatan (confirmed at Milestone 5 close-out, not just
    % suspected): griddatan('nearest',...) always triangulates via
    % delaunayn and locates the nearest SIMPLEX VERTEX via dsearchn, not a
    % global nearest-neighbor search -- and it pre-averages exact-duplicate
    % training rows before triangulating. This k-d tree does neither, so it
    % is only proven to match griddatan when no two training rows are exact
    % duplicates (true for the 3000-row corpus this project uses -- min
    % pairwise distance 0.00994, verified at Milestone 5 close-out). Warn
    % rather than silently diverge if that ever stops being true.
    if size(unique(X_train, 'rows'), 1) < size(X_train, 1)
        warning('fit_nearest_predictor:DuplicateRows', ...
            ['Training data contains duplicate feature rows. griddatan ' ...
             'averages duplicates'' targets before use; this k-d tree ' ...
             'predictor does not. Re-run the Milestone 5 parity checks ' ...
             'before trusting predictions from this model.']);
    end
    model.ns = createns(X_train, 'NSMethod', 'kdtree');
    model.y = y_train;
end

function y_pred = predict_nearest(model, X_query)
    % Exact ties (query equidistant from 2+ training points): confirmed at
    % Milestone 5 close-out to disagree with griddatan roughly half the
    % time when forced (griddatan locates the nearest vertex of the
    % enclosing Delaunay simplex via dsearchn; knnsearch does a global
    % distance search -- the two have no reason to agree on a tie). Not
    % observed on any toy or real query tested here (diff=0 throughout),
    % but not something a k-d tree can be made to match in general. Re-run
    % the parity checks in this file if a new dataset is substituted in.
    idx = knnsearch(model.ns, X_query);
    y_pred = model.y(idx);
end

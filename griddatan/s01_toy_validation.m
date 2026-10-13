% s01_toy_validation.m
%
% Milestone 4: sanity-check griddatan itself on a synthetic function with a
% known closed-form answer, in a dimensionality comparable to the real
% problem (6-D), before trusting it on real data. Also independently
% re-verifies the one finding already on record from the old draft
% (interpolant.m): that 'linear' produces NaN for query points outside the
% training data's convex hull.

clc; clear;
rng(7);   % arbitrary, fixed -- this script's own data is synthetic, not
          % the shared training corpus, so it doesn't need to match GP's seed

%% A known 6-D function, cheap to sample exactly
% Smooth, non-trivial (not separable, not purely linear). NOT dimension-
% balanced, though: a variance decomposition shows x3*x4 and x5^2 alone
% explain ~98% of f's variance, so the NMSE numbers below are softer
% evidence than they look -- they mainly tell you "nearest-neighbor
% distance works" on the ~3 dimensions that actually move the output. The
% convex-hull / NaN-rate finding this script exists to check does NOT
% depend on f's value structure (it's pure point-cloud geometry using all 6
% dimensions), so that finding is unaffected by this imbalance. See
% STATUS.md, Milestone 4, for the full caveat.
f = @(X) sin(X(:,1)) + 0.5*cos(2*X(:,2)) + X(:,3).*X(:,4) - 0.3*X(:,5).^2 + log1p(X(:,6).^2);

n_train = 2000;    % same order of magnitude as the real 6-D training sizes used elsewhere
n_test  = 500;

X_train = rand(n_train, 6) * 2*pi;   % uniform over a bounded box: [0, 2*pi]^6
y_train = f(X_train);

% Interior test points: sampled strictly inside the training box, so they
% are guaranteed inside (or very near) the training convex hull.
X_test_interior = rand(n_test, 6) * 2*pi;
y_test_interior = f(X_test_interior);

% Exterior test points: deliberately outside the training box in every
% dimension, to reproduce the convex-hull NaN issue interpolant.m found.
X_test_exterior = 2*pi + rand(n_test, 6) * 2*pi;   % shifted fully outside [0,2*pi]^6

%% Nearest-neighbor griddatan -- the method the old draft kept
tic;
y_pred_nearest = griddatan(X_train, y_train, X_test_interior, 'nearest');
t_nearest = toc;
nmse_nearest = mean((y_test_interior - y_pred_nearest).^2) / var(y_test_interior);
fprintf('[nearest, interior] NMSE = %g, wall time = %.1f s (%d train x %d query, 6-D)\n', ...
    nmse_nearest, t_nearest, n_train, n_test);
assert(nmse_nearest < 1, ...
    'nearest-neighbor interpolation should beat predicting the mean (NMSE < 1) on a smooth toy function -- it does not.');

%% Linear griddatan -- the method the old draft dropped
tic;
y_pred_linear_interior = griddatan(X_train, y_train, X_test_interior, 'linear');
t_linear_interior = toc;
n_nan_interior = sum(isnan(y_pred_linear_interior));
fprintf('[linear, interior] NaN count: %d / %d (%.1f%%), wall time = %.1f s\n', ...
    n_nan_interior, n_test, 100*n_nan_interior/n_test, t_linear_interior);

nmse_linear = mean((y_test_interior(~isnan(y_pred_linear_interior)) - ...
    y_pred_linear_interior(~isnan(y_pred_linear_interior))).^2) / var(y_test_interior);
fprintf('[linear, interior] NMSE (excluding NaN): %g\n', nmse_linear);

tic;
y_pred_linear_exterior = griddatan(X_train, y_train, X_test_exterior, 'linear');
t_linear_exterior = toc;
n_nan_exterior = sum(isnan(y_pred_linear_exterior));
fprintf('[linear, exterior] NaN count: %d / %d (%.1f%%), wall time = %.1f s\n', ...
    n_nan_exterior, n_test, 100*n_nan_exterior/n_test, t_linear_exterior);

%% Verdict
% Thresholds (50% / 10%) are arbitrary but not fragile in practice: observed
% exterior NaN rate is always ~100% (vs. a 50% bar) and observed interior
% NaN rate ranges ~30-50% across seeds/training sizes (vs. a 10% bar) -- see
% STATUS.md, Milestone 4, for the sensitivity check.
fprintf('\n--- Milestone 4 verdict ---\n');
if n_nan_exterior > 0.5 * n_test
    fprintf(['CONFIRMED: ''linear'' griddatan poisons output with NaN for ' ...
        'query points outside the training convex hull in 6-D, matching ' ...
        'interpolant.m''s original finding. ''nearest'' has no such failure ' ...
        'mode (always returns a value) and is the only method carried ' ...
        'forward past this milestone.\n']);
else
    fprintf(['UNEXPECTED: exterior NaN rate (%d/%d) is lower than ' ...
        'interpolant.m''s finding would suggest -- re-examine before ' ...
        'assuming the old finding still holds.\n'], n_nan_exterior, n_test);
end
if n_nan_interior > 0.1 * n_test
    fprintf(['ADDITIONAL FINDING: ''linear'' also fails for a substantial ' ...
        'fraction (%.1f%%) of points sampled INSIDE the training bounding ' ...
        'box -- in 6-D, a random point cloud''s convex hull covers much ' ...
        'less volume than its bounding box, so "interior to the box" does ' ...
        'not imply "interior to the hull." This is expected to transfer ' ...
        'directly to the real 6-D validation set at Milestone 6.\n'], ...
        100*n_nan_interior/n_test);
end
fprintf(['Total griddatan wall time for this toy problem (%d train points): ' ...
    '%.1f s across 3 calls. Scaling behavior as training size grows is the ' ...
    'subject of Milestone 5 (the reusable triangulation-based predictor).\n'], ...
    n_train, t_nearest + t_linear_interior + t_linear_exterior);

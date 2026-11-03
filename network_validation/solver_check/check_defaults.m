% The edited simulator with flux_per_edge and lambda2_w1 at their defaults must
% reproduce the stored runs: a short numeric run (SNTM_SMOKE, the first 5 time
% points of the full grid) against the stored 180-day run at the same points.
% Run with SNTM_SMOKE=1, SNTM_OUT_DIR set to a scratch folder, and neither
% SNTM_FLUX_PER_EDGE nor SNTM_LAMBDA2_W1 set.
here = fileparts(mfilename('fullpath'));
nv = fileparts(here);
cfg = sntm_config();
assert(~cfg.flux_per_edge && ~cfg.lambda2_w1_is_lambda1, 'unset SNTM_FLUX_PER_EDGE and SNTM_LAMBDA2_W1');
run_sntm_cases(fullfile(nv, 'manifests', 'recovery_truth_l2f.csv'), 1, 1);
id = 'rec_l2f__numeric__ECl__g0.0005_l0.02_d80_e20.mat';
s = load(fullfile(cfg.out_dir, id), 'N', 'F', 'W1');
r = load(fullfile(nv, 'outputs_recovery', id), 'N', 'F', 'W1');
k = size(s.N, 2);
fprintf('max |N - stored| over %d time points: %g; F: %g; W1 (first %d steps): %g\n', k, ...
    max(abs(s.N - r.N(:, 1:k)), [], 'all'), max(abs(s.F - r.F(:, :, 1:k)), [], 'all'), k - 1, ...
    max(abs(s.W1(:, :, 1:k-1) - r.W1(:, :, 1:k-1)), [], 'all'));

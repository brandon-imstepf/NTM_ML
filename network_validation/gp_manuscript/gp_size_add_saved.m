function row = gp_size_add_saved(model_file, n)
% Adds a GP fitted elsewhere (gaussian_process-interpolant/gaussian_process.m,
% same configuration as the paper's GP: exact, linear basis, ARD Matern 3/2,
% fminunc, on randperm(1e5, n) rows of data_bias_e5_1.csv after rng(42)) to the
% GP training-size curve: training, test (bias_val_e3.csv) and held-out
% (data_bias_e5_11.csv) NMSE, as supp_gp_griddatan.m computes them (Eq. 7; mean
% and s.d. over 10^3 bootstrap resamples of 10^3 samples, rng(1); _all = whole
% set). Appends or replaces the row gp_<n> of results/supp_gp_griddatan.csv and
% prints the fitted noise s.d. and ARD length scales.
%
%   gp_size_add_saved(fullfile(ntm_ml_root, 'gaussian_process-interpolant', ...
%       'saved_models_superset_7500', 'FM_exact__BF_linear__KF_ARDMatern32__OPT_fminunc.mat'), 7500)

here = fileparts(mfilename('fullpath'));
nv = fileparts(here);
addpath(nv);
cfg = sntm_config();
root = cfg.ntm_ml_root;
super = readmatrix(fullfile(root, 'data', 'data_neural_network_model', 'training_data_e6', 'data_bias_e5_1.csv'));
test = readmatrix(fullfile(root, 'data', 'data_baseline_models', 'e3_data', 'bias_val_e3.csv'));
heldout = readmatrix(fullfile(root, 'data', 'data_neural_network_model', 'training_data_e6', 'data_bias_e5_11.csv'));
rng(42);  train = super(randperm(1e5, n), :);

M = load(model_file, 'stzd_flux_gprModel', 'stzd_emc_gprModel', 'flux_mean', 'flux_std', 'emc_mean', 'emc_std');
assert(M.stzd_flux_gprModel.NumObservations == n, 'model has %d observations, not %d', ...
    M.stzd_flux_gprModel.NumObservations, n);
assert(abs(M.flux_mean - mean(train(:, 7))) < 1e-12 * abs(M.flux_mean) + eps, 'training rows do not match the model');
f = @(X) [predict(M.stzd_flux_gprModel, X) * M.flux_std + M.flux_mean, ...
          predict(M.stzd_emc_gprModel, X) * M.emc_std + M.emc_mean];

row = table("gp_" + n, n, 'VariableNames', ["model", "n_train"]);
sets = struct('train', train, 'test', test, 'heldout', heldout);
for s = ["train", "test", "heldout"]
    D = sets.(s);  P = f(D(:, 1:6));
    tg = ["flux", "emc"];
    for t = 1:2
        [m, sd] = boot(P(:, t), D(:, 6 + t));
        row.(s + "_" + tg(t) + "_nmse_mean") = m;
        row.(s + "_" + tg(t) + "_nmse_sd") = sd;
        row.(s + "_" + tg(t) + "_nmse_all") = nmse_eq7(P(:, t), D(:, 6 + t));
    end
    fprintf('%-8s NMSE flux %.4g (all %.4g), emc %.4g (all %.4g)\n', s, row.(s + "_flux_nmse_mean"), ...
        row.(s + "_flux_nmse_all"), row.(s + "_emc_nmse_mean"), row.(s + "_emc_nmse_all"));
end
row.test_frac_outside_hull = NaN;
for t = ["flux", "emc"]
    G = M.("stzd_" + t + "_gprModel");
    fprintf('%s: noise s.d. %.4g, length scales %s\n', t, G.Sigma, mat2str(G.KernelInformation.KernelParameters(1:6)', 3));
end

file = fullfile(nv, 'results', 'supp_gp_griddatan.csv');
T = readtable(file, 'TextType', 'string');
T(T.model == row.model, :) = [];
T = [T; row(:, T.Properties.VariableNames)];
writetable(T, file);
fprintf('wrote %s (row %s)\n', file, row.model);
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

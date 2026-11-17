function T = gp_sigma_curve(sigmas)
% Complexity curve of the GP surrogate: the paper's GP configuration (exact,
% linear basis, ARD Matern 3/2, fminunc; 5,000 rows randperm(1e5, 5000) of
% data_bias_e5_1.csv after rng(42)) refitted with the noise standard deviation
% fixed ('ConstantSigma') at each value of sigmas (units of the Z-scored target).
% Small noise: the GP nearly interpolates its training data (flexible); large
% noise: it smooths (the analogue of the number of features for linear
% regression, the complexity for symbolic regression and the epochs for the MLP).
% The paper's GP is the point sigma = 0.01: its fitted noise ended at that value,
% the default lower bound of fitrgp, and it is used as is (not refitted).
% Training, test (bias_val_e3.csv) and held-out (data_bias_e5_11.csv) NMSE as in
% supp_gp_griddatan.m: mean and s.d. over 10^3 bootstrap resamples of 10^3.
% Fits are cached in gp_manuscript/models_size/gp_sigma_<s>.mat. Writes
% results/gp_sigma_curve.csv (rows added as fits finish).
%
%   T = gp_sigma_curve([1e-4 1e-3 1e-2 3e-2 0.1 0.3 1]);

if nargin < 1 || isempty(sigmas), sigmas = [1e-4 1e-3 1e-2 3e-2 0.1 0.3 1]; end
here = fileparts(mfilename('fullpath'));
nv = fileparts(here);
addpath(nv);
cfg = sntm_config();
root = cfg.ntm_ml_root;
super = readmatrix(cfg.train_superset_file);
test = readmatrix(fullfile(root, 'data', 'data_baseline_models', 'e3_data', 'bias_val_e3.csv'));
heldout = readmatrix(fullfile(root, 'data', 'data_neural_network_model', 'training_data_e6', 'data_bias_e5_11.csv'));
rng(42);  D = super(randperm(1e5, 5000), :);
G = load(cfg.gp_model_file);
out_csv = fullfile(nv, 'results', 'gp_sigma_curve.csv');
rows = {};
for s = sigmas
    if abs(s - 0.01) < 1e-12
        M = struct('flux', G.stzd_flux_gprModel, 'emc', G.stzd_emc_gprModel, 'flux_mean', G.flux_mean, ...
            'flux_std', G.flux_std, 'emc_mean', G.emc_mean, 'emc_std', G.emc_std, 'fit_time_s', NaN);
    else
        f = fullfile(here, 'models_size', sprintf('gp_sigma_%g.mat', s));
        if isfile(f)
            M = load(f);
        else
            M.flux_mean = mean(D(:, 7));  M.flux_std = std(D(:, 7));
            M.emc_mean = mean(D(:, 8));   M.emc_std = std(D(:, 8));
            opts = {'Standardize', true, 'FitMethod', 'exact', 'BasisFunction', 'linear', ...
                'KernelFunction', 'ardmatern32', 'Optimizer', 'fminunc', ...
                'Sigma', s, 'ConstantSigma', true, 'SigmaLowerBound', s / 10};
            t0 = tic;
            M.flux = compact(fitrgp(D(:, 1:6), (D(:, 7) - M.flux_mean) / M.flux_std, opts{:}));
            M.emc = compact(fitrgp(D(:, 1:6), (D(:, 8) - M.emc_mean) / M.emc_std, opts{:}));
            M.fit_time_s = toc(t0);
            save(f, '-struct', 'M');
        end
    end
    pf = @(X) [predict(M.flux, X) * M.flux_std + M.flux_mean, predict(M.emc, X) * M.emc_std + M.emc_mean];
    row = {s, M.fit_time_s};
    sets = {D, test, heldout};
    for k = 1:3
        P = pf(sets{k}(:, 1:6));
        for t = 1:2
            [m, sd] = boot(P(:, t), sets{k}(:, 6 + t));
            row = [row, {m, sd}]; %#ok<AGROW>
        end
    end
    rows(end+1, :) = row; %#ok<AGROW>
    T = cell2table(rows, 'VariableNames', {'sigma', 'fit_time_s', 'train_flux_mean', 'train_flux_sd', ...
        'train_emc_mean', 'train_emc_sd', 'test_flux_mean', 'test_flux_sd', 'test_emc_mean', 'test_emc_sd', ...
        'heldout_flux_mean', 'heldout_flux_sd', 'heldout_emc_mean', 'heldout_emc_sd'});
    writetable(T, out_csv);
    fprintf('sigma %g: train %.3g / %.3g, held-out %.4g / %.4g (fit %.0f s)\n', s, row{3}, row{5}, row{11}, row{13}, row{2});
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

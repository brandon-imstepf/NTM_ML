function T = gp_single_edge_table()
% Single-edge accuracy of the GP surrogate for Table 1 of GP_manuscript.tex,
% and its ARD length scales (feature relevance).
% NMSE = MSE / Var(target) (Eq. 7), as mean +- s.d. over 10^3 bootstrap
% resamples of 10^3 samples drawn with replacement (rng 1) from
%   train: the GP's 5,000 training rows (rows of data_bias_e5_1.csv chosen with
%          rng(training_seed), randperm(1e5, dataset_size), as in
%          gaussian_process-interpolant/gaussian_process.m)
%   test:  bias_val_e3.csv (10^3 held-out samples, as for LR and SR). This is
%          also the set on which gaussian_process.m ranked the GP configurations
%   heldout_e5_11: data_bias_e5_11.csv (10^5 samples, never used for fitting or
%          selection), to check the test NMSE for selection bias
% Length scales of the ARD Matern 3/2 kernels are for standardized inputs
% (fitrgp 'Standardize', true); a shorter length scale = a more relevant input.
% Writes results/gp_single_edge_table.csv and results/gp_ard_lengthscales.csv.
%
%   T = gp_single_edge_table();

here = fileparts(mfilename('fullpath'));
nv = fileparts(here);
addpath(nv);
cfg = sntm_config();
G = load(cfg.gp_model_file);
root = cfg.ntm_ml_root;

S = readmatrix(fullfile(root, 'data', 'data_neural_network_model', 'training_data_e6', 'data_bias_e5_1.csv'));
rng(G.training_seed);
train = S(randperm(1e5, G.dataset_size), :);
test = readmatrix(fullfile(root, 'data', 'data_baseline_models', 'e3_data', 'bias_val_e3.csv'));
heldout = readmatrix(fullfile(root, 'data', 'data_neural_network_model', 'training_data_e6', 'data_bias_e5_11.csv'));

predict_flux = @(X) predict(G.stzd_flux_gprModel, X) * G.flux_std + G.flux_mean;
predict_emc = @(X) predict(G.stzd_emc_gprModel, X) * G.emc_std + G.emc_mean;
nmse = @(p, y) mean((p - y).^2) / var(y);

rows = {};
sets = {train, 'train'; test, 'test'; heldout, 'heldout_e5_11'};
for k = 1:3
    D = sets{k, 1};
    pf = chunked(predict_flux, D(:, 1:6));  pe = chunked(predict_emc, D(:, 1:6));
    rng(1);
    b = zeros(1000, 2);
    for i = 1:1000
        j = randi(size(D, 1), 1000, 1);
        b(i, :) = [nmse(pf(j), D(j, 7)), nmse(pe(j), D(j, 8))];
    end
    rows(end+1, :) = {sets{k, 2}, size(D, 1), mean(b(:, 1)), std(b(:, 1)), mean(b(:, 2)), std(b(:, 2)), ...
        nmse(pf, D(:, 7)), nmse(pe, D(:, 8))}; %#ok<AGROW>
end
T = cell2table(rows, 'VariableNames', {'set', 'n', 'flux_nmse_mean', 'flux_nmse_sd', 'emc_nmse_mean', 'emc_nmse_sd', ...
    'flux_nmse_all', 'emc_nmse_all'});
writetable(T, fullfile(nv, 'results', 'gp_single_edge_table.csv'));
disp(T)

names = {'gamma', 'lambda', 'delta', 'epsilon', 'N1', 'N2'};
lf = G.stzd_flux_gprModel.KernelInformation.KernelParameters(1:6);
le = G.stzd_emc_gprModel.KernelInformation.KernelParameters(1:6);
L = table(names', lf(:), le(:), 'VariableNames', {'feature', 'flux_lengthscale', 'emc_lengthscale'});
writetable(L, fullfile(nv, 'results', 'gp_ard_lengthscales.csv'));
disp(L)
fprintf('kernel %s, basis %s, fit %s\n', G.stzd_flux_gprModel.KernelInformation.Name, ...
    G.stzd_flux_gprModel.BasisFunction, G.stzd_flux_gprModel.FitMethod);
end

function y = chunked(f, X)
% Predict in blocks of 10^4 rows (an exact GP's kernel against 5,000 training rows)
y = zeros(size(X, 1), 1);
for a = 1:10000:size(X, 1)
    b = min(a + 9999, size(X, 1));
    y(a:b) = f(X(a:b, :));
end
end

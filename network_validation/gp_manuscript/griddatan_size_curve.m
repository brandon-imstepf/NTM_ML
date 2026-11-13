function T = griddatan_size_curve()
% Held-out error of griddatan 'nearest' (the nearest training input; k-d tree,
% identical to griddatan(..., 'nearest'), griddatan/README.md milestone 5)
% against the number of training samples, for three scalings of the inputs:
%   raw     the inputs as they are (griddatan/s03; the 1.83 of the response)
%   zscore  each input z-scored with the training mean and s.d. (the "interp"
%           sNTM before 11/13, n = 10^4)
%   logz    log10 gamma1 and log10 lambda1, then z-scored (the interp sNTM, n = 10^6)
% Training rows: for n <= 10^5, randperm(10^5, n) of data_bias_e5_1.csv after
% rng(42) (the GP's sampling; n = 10^5 is the whole file); for n > 10^5,
% randperm(10^6, n) after rng(42) of the MLP's training data data_bias_e5_1..10
% (n = 10^6 is all of it; interp_training_data.m). Held out: data_bias_e5_11.csv
% (10^5 samples never used for fitting or selection). Validation: the 10^3-sample
% bias_val_e3.csv, on which the scaling is chosen (val_*_nmse).
% NMSE = Eq. 7; mean and s.d. over 10^3 bootstrap resamples of 10^3 held-out
% samples (rng(1)), as supp_gp_griddatan.m, and over all 10^5 (_all).
% Training NMSE is 0 by construction (a training input is its own nearest
% neighbour; the training inputs are distinct) and is not tabulated.
% Writes results/griddatan_size_curve.csv.
%
%   T = griddatan_size_curve();

here = fileparts(mfilename('fullpath'));
nv = fileparts(here);
addpath(nv);
cfg = sntm_config();
super = readmatrix(cfg.train_superset_file);
held = readmatrix(fullfile(cfg.ntm_ml_root, 'data', 'data_neural_network_model', 'training_data_e6', 'data_bias_e5_11.csv'));
val = readmatrix(fullfile(cfg.ntm_ml_root, 'data', 'data_baseline_models', 'e3_data', 'bias_val_e3.csv'));
full = cell(10, 1);
for k = 1:10, full{k} = readmatrix(cfg.interp_train_files{k}); end
full = vertcat(full{:});
sizes = [500 1000 2000 5000 10000 20000 50000 100000 200000 500000 1000000];
scalings = ["raw", "zscore", "logz"];
rows = {};
for n = sizes
    rng(42);
    if n <= 1e5, D = super(randperm(size(super, 1), n), :); else, D = full(randperm(size(full, 1), n), :); end
    for sc = scalings
        tf = scaler(sc, D(:, 1:6));
        t0 = tic;  ns = createns(tf(D(:, 1:6)), 'NSMethod', 'kdtree');  t_fit = toc(t0);
        t0 = tic;  P = D(knnsearch(ns, tf(held(:, 1:6))), 7:8);  t_query = toc(t0) / size(held, 1);
        [mf, sf] = boot(P(:, 1), held(:, 7));  [me, se] = boot(P(:, 2), held(:, 8));
        Pv = D(knnsearch(ns, tf(val(:, 1:6))), 7:8);
        rows(end+1, :) = {sc, n, mf, sf, me, se, nmse_eq7(P(:, 1), held(:, 7)), nmse_eq7(P(:, 2), held(:, 8)), ...
            t_fit, t_query, nmse_eq7(Pv(:, 1), val(:, 7)), nmse_eq7(Pv(:, 2), val(:, 8))}; %#ok<AGROW>
        fprintf('%-6s n %6d: held-out NMSE flux %.3g, emc %.3g\n', sc, n, mf, me);
    end
end
T = cell2table(rows, 'VariableNames', {'scaling', 'n_train', 'heldout_flux_nmse_mean', 'heldout_flux_nmse_sd', ...
    'heldout_emc_nmse_mean', 'heldout_emc_nmse_sd', 'heldout_flux_nmse_all', 'heldout_emc_nmse_all', ...
    'fit_time_s', 'query_time_s', 'val_flux_nmse', 'val_emc_nmse'});
writetable(T, fullfile(nv, 'results', 'griddatan_size_curve.csv'));
end

function tf = scaler(sc, ref)
switch sc
    case "raw",    tf = @(X) X;
    case "zscore", mu = mean(ref);  sd = std(ref);  tf = @(X) (X - mu) ./ sd;
    case "logz"
        lg = @(X) [log10(X(:, 1:2)), X(:, 3:end)];
        mu = mean(lg(ref));  sd = std(lg(ref));  tf = @(X) (lg(X) - mu) ./ sd;
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

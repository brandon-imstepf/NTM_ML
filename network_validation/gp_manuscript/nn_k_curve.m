function T = nn_k_curve(ks)
% Complexity curve of the nearest-neighbor surrogate: training and held-out NMSE
% (Eq. 7) of k-nearest-neighbor averaging (the mean of phi_tau, Delta_tau over the
% k nearest training samples; k = 1 is the nearest-neighbor sNTM) against k, on the
% sNTM's inputs (log10 gamma1, log10 lambda1, then z-scored) and lookup table (the
% 10^6 MLP training samples, interp_training_data.m). Small k = flexible (k = 1
% reproduces the training data), large k = smooth; the analogue of the number of
% features (LR), complexity (SR) or epochs (MLP).
% Training NMSE on 10^4 of the training samples (each is its own first neighbor);
% held-out NMSE on data_bias_e5_11.csv (10^5). Mean and s.d. over 10^3 bootstrap
% resamples of 10^3 samples (rng(1)), as for the other surrogates.
% Writes results/nn_k_curve.csv.
%
%   T = nn_k_curve([1 2 3 5 8 12 20 30 50 80 120]);

if nargin < 1 || isempty(ks), ks = [1 2 3 4 5 6 8 10 15 20 30 50 80 120]; end
here = fileparts(mfilename('fullpath'));
nv = fileparts(here);
addpath(nv);
cfg = sntm_config();
D = interp_training_data(cfg);
held = readmatrix(fullfile(cfg.ntm_ml_root, 'data', 'data_neural_network_model', 'training_data_e6', 'data_bias_e5_11.csv'));
L = cfg.interp_log_inputs;
lg = @(X) [X(:, 1:6)] .* 1;  %#ok<NASGU>
tf0 = @(X) logc(X(:, 1:6), L);
Z = tf0(D);  mu = mean(Z);  sd = std(Z);
tf = @(X) (tf0(X) - mu) ./ sd;
ns = createns(tf(D), 'NSMethod', 'kdtree');
rng(7);  tr = D(randperm(size(D, 1), 1e4), :);
K = max(ks);
idx_h = knnsearch(ns, tf(held), 'K', K);
idx_t = knnsearch(ns, tf(tr), 'K', K);
rows = {};
for k = ks
    for t = 1:2
        Yv = D(:, 6 + t);
        ph = mean(Yv(idx_h(:, 1:k)), 2);  pt = mean(Yv(idx_t(:, 1:k)), 2);
        [mh, sh] = boot(ph, held(:, 6 + t));  [mt, st] = boot(pt, tr(:, 6 + t));
        rows(end+1, :) = {k, string(t), mt, st, mh, sh, nmse_eq7(pt, tr(:, 6 + t)), nmse_eq7(ph, held(:, 6 + t))}; %#ok<AGROW>
    end
    fprintf('k %3d: held-out NMSE flux %.4f emc %.4f, training flux %.4g emc %.4g\n', k, rows{end - 1, 5}, ...
        rows{end, 5}, rows{end - 1, 3}, rows{end, 3});
end
T = cell2table(rows, 'VariableNames', {'k', 'target', 'train_nmse_mean', 'train_nmse_sd', 'heldout_nmse_mean', ...
    'heldout_nmse_sd', 'train_nmse_all', 'heldout_nmse_all'});
T.target = replace(T.target, ["1", "2"], ["flux", "emc"]);
writetable(T, fullfile(nv, 'results', 'nn_k_curve.csv'));
end

function X = logc(X, L)
X(:, L) = log10(X(:, L));
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

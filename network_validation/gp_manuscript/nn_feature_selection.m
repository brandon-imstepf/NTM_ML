function T = nn_feature_selection()
% Forward input selection for the nearest-neighbor surrogate, the analogue of the
% forward feature selection of linear regression (Fig. 3 of the main text): for
% phi_tau and Delta_tau separately, start from no inputs and at each step add the
% input that most lowers the leave-one-out NMSE on the training data, until all
% six inputs are used (the nearest-neighbor sNTM).
% Model: 1-nearest neighbor on the sNTM's inputs (log10 gamma1, log10 lambda1, the
% others as they are; each z-scored with the training mean and s.d.) restricted to
% the selected inputs; lookup table: the 10^6 MLP training samples
% (interp_training_data.m).
% Leave-one-out NMSE: 10^4 training samples (rng(7)), each predicted by its
% nearest *other* training sample (second neighbor). Training NMSE is 0 for
% every input set (each training sample is its own nearest neighbor).
% Held-out NMSE: data_bias_e5_11.csv (10^5). NMSE = Eq. 7; mean and s.d. over
% 10^3 bootstrap resamples of 10^3 samples (rng(1)).
% Writes results/nn_feature_selection.csv (one row per target and step).
%
%   T = nn_feature_selection();

here = fileparts(mfilename('fullpath'));
nv = fileparts(here);
addpath(nv);
cfg = sntm_config();
D = interp_training_data(cfg);
held = readmatrix(fullfile(cfg.ntm_ml_root, 'data', 'data_neural_network_model', 'training_data_e6', 'data_bias_e5_11.csv'));
names = ["gamma1", "lambda1", "delta", "epsilon", "N1", "N2"];
L = cfg.interp_log_inputs;
Z = logc(D(:, 1:6), L);  mu = mean(Z);  sd = std(Z);
tf = @(X) (logc(X(:, 1:6), L) - mu) ./ sd;
ZD = tf(D);  ZH = tf(held);
rng(7);  itr = randperm(size(D, 1), 1e4);
targets = ["flux", "emc"];
rows = {};
for t = 1:2
    y = D(:, 6 + t);
    sel = [];  rest = 1:6;
    for step = 1:6
        best = Inf;
        for c = rest
            cols = [sel c];
            ns = createns(ZD(:, cols), 'NSMethod', 'kdtree');
            idx = knnsearch(ns, ZD(itr, cols), 'K', 2);
            p = loo_pred(idx, itr, y);
            e = nmse_eq7(p, y(itr));
            if e < best, best = e;  bc = c;  end
        end
        sel = [sel bc];  rest = setdiff(rest, bc);
        ns = createns(ZD(:, sel), 'NSMethod', 'kdtree');
        idx = knnsearch(ns, ZD(itr, sel), 'K', 2);
        p_loo = loo_pred(idx, itr, y);
        p_h = y(knnsearch(ns, ZH(:, sel)));
        [ml, sl] = boot(p_loo, y(itr));  [mh, sh] = boot(p_h, held(:, 6 + t));
        rows(end+1, :) = {targets(t), step, names(bc), strjoin(names(sel), "+"), ml, sl, mh, sh, ...
            nmse_eq7(p_loo, y(itr)), nmse_eq7(p_h, held(:, 6 + t))}; %#ok<AGROW>
        fprintf('%s step %d: + %-8s leave-one-out %.4g, held-out %.4g\n', targets(t), step, names(bc), ml, mh);
    end
end
T = cell2table(rows, 'VariableNames', {'target', 'n_inputs', 'added', 'inputs', 'loo_nmse_mean', 'loo_nmse_sd', ...
    'heldout_nmse_mean', 'heldout_nmse_sd', 'loo_nmse_all', 'heldout_nmse_all'});
writetable(T, fullfile(nv, 'results', 'nn_feature_selection.csv'));
end

function p = loo_pred(idx, itr, y)
% nearest other sample: the first neighbor unless it is the query itself
nb = idx(:, 1);
self = nb == itr(:);
nb(self) = idx(self, 2);
p = y(nb);
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

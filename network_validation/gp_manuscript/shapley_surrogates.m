function shapley_surrogates(models, n_query, n_back)
% Shapley values of the inputs of the GP and of the nearest-neighbour
% interpolant (griddatan 'nearest'; the "interp" sNTM of make_surrogate.m),
% for phi_tau and Delta_tau.
% Exact interventional Shapley values: for 6 inputs all 2^6 = 64 coalitions are
% evaluated, the inputs outside a coalition taken from n_back background
% samples (the value of a coalition is the mean prediction over them):
%   phi_i(x) = sum over S not containing i of |S|! (6-|S|-1)! / 6! [v(S+i) - v(S)]
% Query points: n_query samples of the held-out set data_bias_e5_11.csv (never
% used for fitting or selection). Background: n_back samples of the model's own
% training rows. Both drawn after rng(1). Shapley values are in the units of the
% target (phi_tau, Delta_tau as in the training data); their sum over the inputs
% is f(x) minus the mean prediction over the background.
% Models
%   gp      the paper's GP (exact, linear basis, ARD Matern 3/2; 5000 rows),
%           sntm_config gp_model_file
%   mlp     the MLP (make_surrogate.m); background from data_bias_e5_1.csv
%   lr, sr  the feature-selected linear regressions and the symbolic regression
%           expressions (make_surrogate.m), trained on bias_train_e3.csv (10^3)
%   interp  nearest neighbour over 10^6 training rows on log10 gamma1, log10
%           lambda1 and the other inputs, all z-scored (make_surrogate.m,
%           cfg.interp_log_inputs)
% Writes results/shapley_<model>_<target>.csv: one row per query point, the 6
% input values and the 6 Shapley values (phi_<input>), f(x) and the base value;
% and the s.d. of each target over the model's training samples in
% results/shapley_target_std.csv (to put the values in Z-scored units).
%
%   shapley_surrogates();                         % gp and interp, 300 x 50
%   shapley_surrogates("interp", 2000, 100);

if nargin < 1 || isempty(models), models = ["gp", "interp"]; end
if nargin < 2 || isempty(n_query), n_query = 300; end
if nargin < 3 || isempty(n_back), n_back = 50; end
here = fileparts(mfilename('fullpath'));
nv = fileparts(here);
addpath(nv);
cfg = sntm_config();
names = ["gamma1", "lambda1", "delta", "epsilon", "N1", "N2"];
heldout = readmatrix(fullfile(cfg.ntm_ml_root, 'data', 'data_neural_network_model', 'training_data_e6', 'data_bias_e5_11.csv'));
super = readmatrix(cfg.train_superset_file);

for m = string(models)
    switch m
        case "gp"
            G = load(cfg.gp_model_file);
            rng(G.training_seed);  train = super(randperm(size(super, 1), G.dataset_size), :);
            f = {@(X) predict(G.stzd_flux_gprModel, X) * G.flux_std + G.flux_mean, ...
                 @(X) predict(G.stzd_emc_gprModel, X) * G.emc_std + G.emc_mean};
        case "interp"
            train = interp_training_data(cfg);          % all 10^6 MLP training samples
            S = make_surrogate('interp', cfg);           % log gamma1, log lambda1, z-scored
            f = {S.predict_flux, S.predict_w1};
        case "mlp"
            % the MLP of the network runs (make_surrogate.m); background from the
            % 10^5 samples data_bias_e5_1.csv, the first 10^5 of its training data
            train = super;
            S = make_surrogate('mlp', cfg);
            f = {S.predict_flux, S.predict_w1};
        case {"lr", "sr"}
            % the feature-selected linear regressions and the PySR expressions of the
            % network runs, both fit on the 10^3-sample training set bias_train_e3.csv
            train = readmatrix(fullfile(cfg.ntm_ml_root, 'data', 'data_baseline_models', 'e3_data', 'bias_train_e3.csv'));
            S = make_surrogate(char(m), cfg);
            f = {S.predict_flux, S.predict_w1};
    end
    rng(1);
    Xq = heldout(randperm(size(heldout, 1), n_query), 1:6);
    Xb = train(randperm(size(train, 1), n_back), 1:6);
    targets = ["flux", "emc"];
    for t = 1:2
        tic;
        [phi, fx, base] = shapley_exact(f{t}, Xq, Xb);
        T = [array2table(Xq, 'VariableNames', names), array2table(phi, 'VariableNames', "phi_" + names), ...
            table(fx, repmat(base, n_query, 1), 'VariableNames', ["f_x", "base"])];
        out = fullfile(nv, 'results', "shapley_" + m + "_" + targets(t) + ".csv");
        writetable(T, out);
        save_target_std(fullfile(nv, 'results', 'shapley_target_std.csv'), m, targets(t), std(train(:, 6 + t)));
        fprintf('%s %s: %.0f s, mean |phi|: %s; additivity error %.2g\n', m, targets(t), toc, ...
            strjoin(compose('%s %.3g', names', mean(abs(phi))'), ', '), max(abs(sum(phi, 2) + base - fx)));
    end
end
end

function [phi, fx, base] = shapley_exact(f, Xq, Xb)
% Exact interventional Shapley values of f (n x 6 -> n x 1) at the rows of Xq
% with background Xb.
[nq, p] = size(Xq);  nb = size(Xb, 1);
S = dec2bin(0:2^p - 1) == '1';                 % coalitions (2^p x p), true = from x
sz = sum(S, 2);
phi = zeros(nq, p);  fx = zeros(nq, 1);
base = mean(f(Xb));
chunk = max(1, floor(2e4 / (2^p * nb)));       % query points per call of f
for q0 = 1:chunk:nq
    qs = q0:min(nq, q0 + chunk - 1);
    X = zeros(numel(qs) * 2^p * nb, p);
    r = 0;
    for q = qs
        for s = 1:2^p
            blk = Xb;  blk(:, S(s, :)) = repmat(Xq(q, S(s, :)), nb, 1);
            X(r + (1:nb), :) = blk;  r = r + nb;
        end
    end
    v = reshape(mean(reshape(f(X), nb, []), 1), 2^p, numel(qs));   % v(S) per query
    for i = 1:p
        without = find(~S(:, i));
        with = arrayfun(@(s) find(all(S == (S(s, :) | (1:p) == i), 2)), without);
        w = factorial(sz(without)) .* factorial(p - sz(without) - 1) / factorial(p);
        phi(qs, i) = (w' * (v(with, :) - v(without, :)))';
    end
    fx(qs) = v(end, :)';
end
end

function save_target_std(file, model, target, v)
if isfile(file), S = readtable(file, 'TextType', 'string');
else, S = table(strings(0, 1), strings(0, 1), zeros(0, 1), 'VariableNames', ["model", "target", "train_std"]); end
S(S.model == model & S.target == target, :) = [];
S(end+1, :) = {string(model), string(target), v};
writetable(S, file);
end

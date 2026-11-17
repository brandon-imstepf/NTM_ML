function T = table1_bootstrap(n_boot, n_samp)
% Table 1 of the paper with one bootstrap procedure for every surrogate
% (R2.M3): NMSE (Eq. 7, nmse_eq7.m) on n_boot resamples (with replacement,
% rng(b) for resample b) of n_samp samples, mean and standard deviation, for
% the training and the test data of each selected surrogate (make_surrogate.m,
% the models of the network simulations):
%   LR (feature selected), SR: train bias_train_e3.csv, test bias_val_e3.csv (10^3 each)
%   GP: its 5,000 training rows (data_bias_e5_1.csv, randperm with the GP's seed),
%       test data_bias_e5_11.csv (10^5 held-out)
%   1-NN interpolation: train data_bias_e5_1.csv (10^5 of its 10^6 lookup-table
%       rows; zero by construction), test data_bias_e5_11.csv
%   MLP: train data_bias_e5_1.csv (10^5 of its 10^6 training rows), test
%       data_bias_e5_11.csv
% The same resample indices are used for every surrogate on a given data set.
% Also gives the NMSE over the whole set (no resampling).
% Writes results/table1_bootstrap.csv.
%
%   T = table1_bootstrap();          % 1000 resamples of 1000

if nargin < 1, n_boot = 1000; end
if nargin < 2, n_samp = 1000; end
here = fileparts(mfilename('fullpath'));
nv = fileparts(here);
addpath(nv);
cfg = sntm_config();
root = cfg.ntm_ml_root;
e3 = fullfile(root, 'data', 'data_baseline_models', 'e3_data');
e6 = fullfile(root, 'data', 'data_neural_network_model', 'training_data_e6');
lrsr_train = readmatrix(fullfile(e3, 'bias_train_e3.csv'));
lrsr_test = readmatrix(fullfile(e3, 'bias_val_e3.csv'));
super = readmatrix(fullfile(e6, 'data_bias_e5_1.csv'));
heldout = readmatrix(fullfile(e6, 'data_bias_e5_11.csv'));
G = load(cfg.gp_model_file);
rng(G.training_seed);
gp_train = super(randperm(1e5, G.dataset_size), :);
methods = ["lr", "sr", "gp", "interp", "mlp"];
sets = struct('lr', {{lrsr_train, lrsr_test}}, 'sr', {{lrsr_train, lrsr_test}}, ...
    'gp', {{gp_train, heldout}}, 'interp', {{super, heldout}}, 'mlp', {{super, heldout}});
setname = ["train", "test"];

% resample indices per data-set size
idx = containers.Map('KeyType', 'double', 'ValueType', 'any');
rows = {};
for m = methods
    sur = make_surrogate(char(m), cfg);
    for k = 1:2
        D = sets.(m){k};
        n = size(D, 1);
        if ~isKey(idx, n)
            I = zeros(n_samp, n_boot);
            for b = 1:n_boot
                rng(b);
                I(:, b) = randi(n, n_samp, 1);
            end
            idx(n) = I;
        end
        I = idx(n);
        p = [chunked(sur.predict_flux, D(:, 1:6)), chunked(sur.predict_w1, D(:, 1:6))];
        y = D(:, 7:8);
        for t = 1:2
            e = zeros(n_boot, 1);
            for b = 1:n_boot
                e(b) = nmse_eq7(p(I(:, b), t), y(I(:, b), t));
            end
            rows(end+1, :) = {m, setname(k), t, n, mean(e), std(e), median(e), nmse_eq7(p(:, t), y(:, t))}; %#ok<AGROW>
        end
    end
end
T = cell2table(rows, 'VariableNames', {'method', 'set', 'target', 'n', 'nmse_mean', 'nmse_sd', 'nmse_median', 'nmse_full'});
T.target = categorical(T.target, [1 2], {'phi_tau', 'Delta_tau'});
writetable(T, fullfile(nv, 'results', 'table1_bootstrap.csv'));
disp(T)
end

function y = chunked(f, X)
y = zeros(size(X, 1), 1);
for a = 1:20000:size(X, 1)
    b = min(a + 19999, size(X, 1));
    y(a:b) = f(X(a:b, :));
end
end

function R = domain_edge_error(n_boot)
% Single-edge error of the four surrogates near the edge of the training
% domain (advisor/reviewer question: does the GP, or LR, fail out of domain
% or where tau concentrates?).
% Held-out data: data_bias_e5_11.csv (10^5 rows, same sampling as the
% training data). No surrogate was fit on it; nn-hyp-tuning.py only records
% the MLP's error on it (early stopping and the saved model use the held-out
% fold). The GP was fit on 5,000 rows of data_bias_e5_1.csv, LR and SR on
% bias_train_e3.csv.
% Position in the domain: for each input k, u_k = empirical CDF of the
% sampling distribution (data_bias_e5_1.csv) at x_k. Edge depth
%   d = min_k min(u_k, 1 - u_k)  in [0, 0.5],
% the tail mass beyond the most extreme input; d < 0.01 means at least one
% input lies in the outer 1% of its range. Rows are binned by d, and, per
% input, by its lower/upper 1% tail.
% Error in a bin: NMSE_bin = mean (p - y)^2 / Var(y) with Var over the bin
% (Eq. 7 restricted to the bin), as mean +- s.d. over n_boot bootstrap
% resamples within the bin; also the MSE relative to the interior bin.
% GP predictive s.d. (flux, phi_tau) per bin.
% Writes results/domain_edge_error.csv (depth bins) and
% results/domain_edge_tails.csv (per-input tails).
%
%   R = domain_edge_error();       % n_boot = 200

if nargin < 1, n_boot = 200; end
here = fileparts(mfilename('fullpath'));
cfg = sntm_config();
dd = fullfile(cfg.ntm_ml_root, 'data', 'data_neural_network_model', 'training_data_e6');
ref = readmatrix(fullfile(dd, 'data_bias_e5_1.csv'));
H = readmatrix(fullfile(dd, 'data_bias_e5_11.csv'));
X = H(:, 1:6);  Y = H(:, 7:8);
inputs = ["gamma", "lambda", "delta", "epsilon", "N1", "N2"];
targets = ["flux", "w1"];
methods = ["mlp", "gp", "sr", "lr"];

% Position of each held-out input in the sampling distribution
U = zeros(size(X));
for k = 1:6
    r = sort(ref(:, k));
    U(:, k) = (discretize_rank(r, X(:, k)) - 0.5) / numel(r);
end
tail = min(U, 1 - U);
d = min(tail, [], 2);
edges = [0, 0.01, 0.05, 0.15, 0.5];
bin_names = ["d<0.01", "0.01-0.05", "0.05-0.15", "interior d>=0.15"];
bin = discretize(d, edges, 'IncludedEdge', 'left');
bin(d >= 0.5) = 4;

P = struct();
for m = methods
    s = make_surrogate(char(m), cfg);
    P.(m) = [chunked(s.predict_flux, X), chunked(s.predict_w1, X)];
end
G = load(cfg.gp_model_file);
gp_sd = zeros(size(X, 1), 1);
for a = 1:10000:size(X, 1)
    b = min(a + 9999, size(X, 1));
    [~, gp_sd(a:b)] = predict(G.stzd_flux_gprModel, X(a:b, :));
end
gp_sd = gp_sd * G.flux_std;

rng(1);
rows = {};
for m = methods
    for t = 1:2
        mse_int = mean((P.(m)(bin == 4, t) - Y(bin == 4, t)).^2);
        for b = 1:4
            idx = find(bin == b);
            [mu, sd] = boot_nmse(P.(m)(idx, t), Y(idx, t), n_boot);
            mse = mean((P.(m)(idx, t) - Y(idx, t)).^2);
            rows(end+1, :) = {m, targets(t), bin_names(b), numel(idx), mu, sd, mse / mse_int, ...
                median(gp_sd(idx))}; %#ok<AGROW>
        end
    end
end
R.depth = cell2table(rows, 'VariableNames', {'method', 'target', 'bin', 'n', 'nmse_mean', 'nmse_sd', ...
    'mse_over_interior', 'gp_flux_sd_median'});

rows = {};
for m = methods
    for t = 1:2
        mse_mid = mean((P.(m)(d >= 0.15, t) - Y(d >= 0.15, t)).^2);
        for k = 1:6
            for side = ["low", "high"]
                if side == "low", idx = find(U(:, k) < 0.01); else, idx = find(U(:, k) > 0.99); end
                [mu, sd] = boot_nmse(P.(m)(idx, t), Y(idx, t), n_boot);
                mse = mean((P.(m)(idx, t) - Y(idx, t)).^2);
                rows(end+1, :) = {m, targets(t), inputs(k), side, numel(idx), mu, sd, mse / mse_mid, ...
                    median(gp_sd(idx))}; %#ok<AGROW>
            end
        end
    end
end
R.tails = cell2table(rows, 'VariableNames', {'method', 'target', 'input', 'tail', 'n', 'nmse_mean', 'nmse_sd', ...
    'mse_over_interior', 'gp_flux_sd_median'});

writetable(R.depth, fullfile(here, 'results', 'domain_edge_error.csv'));
writetable(R.tails, fullfile(here, 'results', 'domain_edge_tails.csv'));
fprintf('held-out rows %d; rows per depth bin: %s\n', size(X, 1), mat2str(accumarray(bin, 1)'));
disp(R.depth)
disp(R.tails)
end

function r = discretize_rank(sorted_ref, x)
% Mid-rank of x in sorted_ref (number of reference values below x, plus half the ties)
[ru, first] = unique(sorted_ref, 'first');
[~, last] = unique(sorted_ref, 'last');
j = discretize(x, [-Inf; ru; Inf]) - 1;     % ru(j) <= x < ru(j+1), j = 0 below all
lo = zeros(size(x));  hi = zeros(size(x));
in = j > 0;
hi(in) = last(j(in));
lo(in) = hi(in);
tie = in & x == ru(max(j, 1));
lo(tie) = first(j(tie)) - 1;
r = (lo + hi) / 2 + 0.5;
end

function y = chunked(f, X)
y = zeros(size(X, 1), 1);
for a = 1:10000:size(X, 1)
    b = min(a + 9999, size(X, 1));
    y(a:b) = f(X(a:b, :));
end
end

function [mu, sd] = boot_nmse(p, y, n_boot)
v = zeros(n_boot, 1);
n = numel(y);
for i = 1:n_boot
    j = randi(n, n, 1);
    v(i) = mean((p(j) - y(j)).^2) / var(y(j));
end
mu = mean(v);  sd = std(v);
end

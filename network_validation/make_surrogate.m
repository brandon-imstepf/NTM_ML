function sur = make_surrogate(method, cfg)
% Build phi_tau / Delta_tau surrogates as function handles that plug into
% NetworkTransportModel's existing 'sr_fun_flux' / 'sr_fun_w1' inputs, so the
% simulator code itself never needs to change.
%
%   sur.flux_fun(g,l,d,e,N1,N2)     -> phi_tau,   same size as N1
%   sur.w1_fun(g,l,d,e,N1,N2,flux)  -> Delta_tau, same size as N1
%   sur.em_fun                      -> zeros (edge mass is not used by the N update)
%
% g,l,d,e are scalars and N1,N2 are nroi x nroi matrices (N1(i,j) = N_i,
% N2(i,j) = N_j). Every model sees raw physical inputs and returns raw outputs.
% Flux is forced to zero on edges where both nodes are empty, matching the
% MLP branch inside NetworkFluxCalculator.
%
% method: 'numeric' | 'mlp' | 'sr' | 'gp' | 'interp' | 'lr'

sur.method = method;
sur.em_fun = @(varargin) zeros(size(varargin{5}));

switch method
    case 'numeric'
        sur.flux_fun = [];
        sur.w1_fun = [];
        return

    case 'sr'
        fs = load(fullfile(cfg.sr_dir, 'flux_scaling_params.mat'));
        es = load(fullfile(cfg.sr_dir, 'emc_scaling_params.mat'));
        zf = @(x, k) (x - fs.feature_means(k)) ./ fs.feature_stds(k);
        ze = @(x, k) (x - es.feature_means(k)) ./ es.feature_stds(k);
        % Final PySR expressions (synthetic_parameter_recovery_flux_only.m)
        flux_z = @(g,l,d,e,N1,N2) -0.487*zf(N2,6) - 0.812 + ...
            ((zf(N1,5) + 3.38).^(-2.64 .* zf(g,1))) ./ (4.96 - zf(d,3));
        emc_z = @(g,l,d,e,N1,N2) ...
            (0.0453*ze(d,3) - ze(g,1) + 0.0453*ze(l,2) - 0.509) ./ (ze(g,1) + 0.691);
        flux_raw = @(g,l,d,e,N1,N2) flux_z(g,l,d,e,N1,N2) * fs.target_std + fs.target_mean;
        % The Delta_tau expression has no N dependence, so broadcast it to matrix size
        w1_raw = @(g,l,d,e,N1,N2) (emc_z(g,l,d,e,N1,N2) * es.target_std + es.target_mean) + zeros(size(N1));
        sur.flux_fun = @(g,l,d,e,N1,N2) gate_flux(flux_raw(g,l,d,e,N1,N2), N1, N2);
        sur.w1_fun = @(g,l,d,e,N1,N2,~) w1_raw(g,l,d,e,N1,N2);

    case 'mlp'
        W = load(cfg.mlp_weights_file);
        predict_f = @(X) mlp_forward(W.flux, X);
        predict_w = @(X) mlp_forward(W.w1, X);
        sur.flux_fun = @(g,l,d,e,N1,N2) gate_flux(apply_rowwise(predict_f,g,l,d,e,N1,N2), N1, N2);
        sur.w1_fun = @(g,l,d,e,N1,N2,~) apply_rowwise(predict_w,g,l,d,e,N1,N2);

    case 'gp'
        % Models predict z-scored targets; the file stores the scaling used in fitting
        G = load(cfg.gp_model_file);
        predict_f = @(X) predict(G.stzd_flux_gprModel, X) * G.flux_std + G.flux_mean;
        predict_w = @(X) predict(G.stzd_emc_gprModel, X) * G.emc_std + G.emc_mean;
        sur.flux_fun = @(g,l,d,e,N1,N2) gate_flux(apply_rowwise(predict_f,g,l,d,e,N1,N2), N1, N2);
        sur.w1_fun = @(g,l,d,e,N1,N2,~) apply_rowwise(predict_w,g,l,d,e,N1,N2);

    case 'interp'
        % Nearest-neighbour lookup on the same 10k subset as interpolant.m.
        % Unlike griddatan on raw inputs, features are z-scored first so that
        % delta/epsilon (~1e1) do not swamp N1/N2 (~1e-3) in the distance.
        T = training_subset(cfg, cfg.interp_train_size);
        mu = mean(T(:,1:6)); sd = std(T(:,1:6));
        tree = KDTreeSearcher((T(:,1:6) - mu) ./ sd);
        yf = T(:,7); yw = T(:,8);
        predict_f = @(X) yf(knnsearch(tree, (X - mu) ./ sd));
        predict_w = @(X) yw(knnsearch(tree, (X - mu) ./ sd));
        sur.flux_fun = @(g,l,d,e,N1,N2) gate_flux(apply_rowwise(predict_f,g,l,d,e,N1,N2), N1, N2);
        sur.w1_fun = @(g,l,d,e,N1,N2,~) apply_rowwise(predict_w,g,l,d,e,N1,N2);

    case 'lr'
        % Feature-selected linear regressions (fit_lr_surrogate.py), fit on raw inputs
        L = load(cfg.lr_weights_file);
        predict_f = @(X) lr_predict(L.flux, X);
        predict_w = @(X) lr_predict(L.w1, X);
        sur.flux_fun = @(g,l,d,e,N1,N2) gate_flux(apply_rowwise(predict_f,g,l,d,e,N1,N2), N1, N2);
        sur.w1_fun = @(g,l,d,e,N1,N2,~) apply_rowwise(predict_w,g,l,d,e,N1,N2);

    otherwise
        error('make_surrogate:method', 'Unknown method "%s".', method);
end
end

function Y = apply_rowwise(predict_fun, g, l, d, e, N1, N2)
% Evaluate a model that takes an n x 6 feature matrix on every (i,j) pair.
n = numel(N1);
X = [repmat([g l d e], n, 1), N1(:), N2(:)];
Y = reshape(predict_fun(X), size(N1));
end

function y = lr_predict(model, X)
% Features are named by input index: 'x1' is X(:,1), 'x1x5' is X(:,1).*X(:,5).
y = model.intercept * ones(size(X,1), 1);
for k = 1:numel(model.features)
    idx = sscanf(strrep(char(model.features{k}), 'x', ' '), '%d');
    y = y + model.coef(k) * prod(X(:, idx), 2);
end
end

function F = gate_flux(F, N1, N2)
F(N1 <= 0 & N2 <= 0) = 0;
end

function T = training_subset(cfg, n)
persistent superset
if isempty(superset)
    superset = readmatrix(cfg.train_superset_file);
end
rng(cfg.train_seed);
idx = randperm(size(superset,1), n);
T = superset(idx, :);
end

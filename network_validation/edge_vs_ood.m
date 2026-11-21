function [S, E] = edge_vs_ood()
% Why are some sNTMs more accurate at the out-of-domain (OOD) settings than at
% the edge settings (main text Table 3b)? For the 8 edge and 8 OOD network
% settings (fig8_cases.m; ECl seed, 180 days), per setting:
%   S (one row per setting and sNTM): the network NMSE (Eq. 7, over all regions
%     and time points, as network_error_all.m) split into its numerator and
%     denominator: MSE of the sNTM concentrations and variance of the numeric
%     ones; finite or not; and the single-edge NMSE of the surrogate's phi_tau
%     and Delta_tau evaluated at the numeric simulation's own states (active
%     edges, all time points), which measures the surrogate's error at these
%     inputs without feedback through the network.
%   E (one row per setting): where the numeric simulation's single-edge inputs
%     fall relative to the training data (N columns of data_bias_e5_*.csv): the
%     fraction of active edge-time inputs with one node empty (N = 0, absent from
%     the training data), and of the nonzero N below the training minimum, below
%     its 1st percentile and above its maximum; plus the numeric total mass at
%     day 180 relative to day 0 and the variance of the numeric N.
% Writes results/edge_vs_ood_nmse.csv and results/edge_vs_ood_inputs.csv.

here = fileparts(mfilename('fullpath'));
addpath(here);  addpath(fullfile(here, 'gp_manuscript'));
cfg = sntm_config();
methods = ["lr", "sr", "interp", "gp", "mlp"];
cases = fig8_cases(["edge", "ood"]);

files = dir(fullfile(cfg.ntm_ml_root, 'data', 'data_neural_network_model', 'training_data_e6', 'data_bias_e5_*.csv'));
Ntrain = [];
for f = files'
    D = readmatrix(fullfile(f.folder, f.name));
    Ntrain = [Ntrain; D(:, 5); D(:, 6)]; %#ok<AGROW>
end
lo = min(Ntrain);  hi = max(Ntrain);  p1 = prctile(Ntrain, 1);
fprintf('training N: %d values, %d zeros, min %.3g, 1st pct %.3g, max %.3g\n', numel(Ntrain), nnz(Ntrain == 0), lo, p1, hi);

sur = struct();
for m = methods, sur.(m) = make_surrogate(char(m), cfg); end
srows = {};  erows = {};
for c = cases
    num = load(fullfile(here, c.num_file), 'N', 'F', 'W1', 'C', 'meta');
    p = num.meta;
    xn = num.N(:);  v = var(xn, 1);
    % single-edge inputs of the numeric simulation (edges the simulator solves)
    [i, j] = find(num.C > 0);
    nt = size(num.N, 2);
    Ni = num.N(i, :);  Nj = num.N(j, :);
    active = Ni > 0 | Nj > 0;
    a = Ni(active);  b = Nj(active);
    x = [a(a > 0); b(b > 0)];
    erows(end+1, :) = {c.group, c.setting, p.gamma1, p.lambda1, p.delta, p.epsilon, nnz(active), ...
        mean(a == 0 | b == 0), mean(x < lo), mean(x < p1), mean(x > hi), median(x), v}; %#ok<AGROW>
    lin = sub2ind(size(num.C), repmat(i, 1, nt), repmat(j, 1, nt));
    tt = repmat(1:nt, numel(i), 1);
    Fn = reshape(num.F, [], nt);  Wn = reshape(num.W1, [], nt);
    Fn = Fn(sub2ind(size(Fn), lin(active), tt(active)));
    Wn = Wn(sub2ind(size(Wn), lin(active), tt(active)));
    w_ok = tt(active) < nt;       % W1 at the last time point is not computed
    for m = methods
        if m == "interp"
            f = fullfile(here, 'outputs_consistent', c.group + "__" + c.setting + "__" + c.seed + "__interp.mat");
        else
            f = fullfile(here, sprintf(c.src_fmt, m));
        end
        s = load(f, 'N');
        fin = all(isfinite(s.N(:)));
        if fin
            mse = mean((s.N(:) - xn).^2);  nm = nmse_eq7(s.N(:), xn);
        else
            mse = NaN;  nm = NaN;
        end
        Fs = sur.(m).flux_fun(p.gamma1, p.lambda1, p.delta, p.epsilon, a, b);
        Ws = sur.(m).w1_fun(p.gamma1, p.lambda1, p.delta, p.epsilon, a, b, []);
        srows(end+1, :) = {c.group, c.setting, m, fin, nm, mse, v, nmse_eq7(Fs, Fn), nmse_eq7(Ws(w_ok), Wn(w_ok))}; %#ok<AGROW>
    end
end
S = cell2table(srows, 'VariableNames', ["group", "setting", "method", "finite", "nmse", "mse", "var_numeric", ...
    "edge_nmse_flux", "edge_nmse_w1"]);
E = cell2table(erows, 'VariableNames', ["group", "setting", "gamma1", "lambda1", "delta", "epsilon", "n_active", ...
    "frac_one_empty", "frac_below_min", "frac_below_p1", "frac_above_max", "median_N", "var_numeric"]);
writetable(S, fullfile(here, 'results', 'edge_vs_ood_nmse.csv'));
writetable(E, fullfile(here, 'results', 'edge_vs_ood_inputs.csv'));
disp(E);
for m = methods
    for g = ["edge", "ood"]
        r = S(S.method == m & S.group == g, :);
        fprintf('%-6s %-4s finite %d/8  median NMSE %.3g  median MSE %.3g  median var %.3g  median single-edge NMSE flux %.3g  W1 %.3g\n', ...
            m, g, nnz(r.finite), median(r.nmse, 'omitnan'), median(r.mse, 'omitnan'), median(r.var_numeric), ...
            median(r.edge_nmse_flux), median(r.edge_nmse_w1));
    end
end
end

function S = single_edge_error_locked()
% Prediction error of each surrogate on the 25 locked single-edge datasets of
% parameter_recovery/locked_edge (numeric phi_tau and Delta_tau on a 32 x 32
% grid of N_1, N_2 in [0, 0.004] at fixed gamma, lambda, delta, epsilon):
% 8 interior sets, 8 on the corners of the training ranges (edge), 8 outside
% them (ood), and the Fig. 8 parameters (paper). For R2.M1.
%   flux_nmse, w1_nmse   sum (s - n)^2 / sum n^2 over the 1024 grid points
% Writes results/single_edge_error_locked.csv.
%
%   S = single_edge_error_locked();

here = fileparts(mfilename('fullpath'));
cfg = sntm_config();
locked = fullfile(cfg.ntm_ml_root, 'parameter_recovery', 'locked_edge');
sets = readtable(fullfile(locked, 'locked_sets.csv'), 'TextType', 'string');
methods = ["mlp", "sr", "gp", "lr"];
sur = struct();
for m = methods, sur.(m) = make_surrogate(char(m), cfg); end

rows = {};
for k = 1:height(sets)
    s = sets(k, :);
    D = readtable(fullfile(locked, 'outputs', s.set_id + ".csv"));
    group = extractBefore(s.set_id, "_");
    if ismissing(group), group = s.set_id; end
    for m = methods
        F = sur.(m).flux_fun(s.gamma1, s.lambda1, s.delta, s.epsilon, D.NRow, D.NCol);
        W = sur.(m).w1_fun(s.gamma1, s.lambda1, s.delta, s.epsilon, D.NRow, D.NCol, []);
        rows(end+1, :) = {s.set_id, group, m, nmse(F, D.FValue), nmse(W, D.W1)}; %#ok<AGROW>
    end
end
S = cell2table(rows, 'VariableNames', {'set_id', 'group', 'method', 'flux_nmse', 'w1_nmse'});
writetable(S, fullfile(here, 'results', 'single_edge_error_locked.csv'));

fprintf('== Single-edge error on the locked sets, median [min, max] per group\n');
for g = ["interior", "edge", "ood", "paper"]
    for m = methods
        r = S(S.group == g & S.method == m, :);
        fprintf('%-8s %-4s flux %-32s W1 %s\n', g, m, mmm(r.flux_nmse), mmm(r.w1_nmse));
    end
end
end

function v = nmse(s, n)
v = sum((s(:) - n(:)).^2) / sum(n(:).^2);
end

function s = mmm(x)
s = sprintf('%.3g [%.3g, %.3g]', median(x, 'omitnan'), min(x), max(x));
end

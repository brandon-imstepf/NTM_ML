function P = gp_trajectory_l2f(every, run_set)
% GP predictive standard deviation along each GP sNTM simulation's own
% trajectory (8 l2f truths, seeds ECl, CA1 and RSPd where present, 180 days),
% for R1.M1: does the
% variance rise when the simulation leaves the training distribution?
% At every 'every'-th time point: the largest N over the regions relative to
% the largest N of the MLP/GP training data, the median and maximum GP flux
% s.d. over the active edges, and the NMSE (Eq. 7, nmse_eq7.m) of the GP run
% against the numeric run up to that time point.
% Writes results/gp_trajectory_l2f.csv; with run_set = "consistent" (l2f_runs.m)
% the runs in outputs_consistent/, into results/gp_trajectory_l2f_consistent.csv.
%
%   P = gp_trajectory_l2f();      % every = 5, run_set "original"
%   P = gp_trajectory_l2f(1);
%   P = gp_trajectory_l2f(5, "consistent");

if nargin < 1, every = 5; end
if nargin < 2, run_set = "original"; end
here = fileparts(mfilename('fullpath'));
cfg = sntm_config();
T = readtable(fullfile(here, 'manifests', 'recovery_truths.csv'), 'TextType', 'string');
T = T(T.convention == "l2f", :);
[seed_dirs, tag] = l2f_runs(run_set);
seeds = ["ECl", "CA1", "RSPd"];
seed_file = @(k, m, id) fullfile(here, seed_dirs(k), "rec_l2f__" + m + "__" + seeds(k) + "__" + id + ".mat");
G = load(cfg.gp_model_file);
Ntrain_max = 0;
for f = dir(fullfile(cfg.ntm_ml_root, 'data', 'data_neural_network_model', 'training_data_e6', 'data_bias_e5_*.csv'))'
    D = readmatrix(fullfile(f.folder, f.name));
    Ntrain_max = max([Ntrain_max; D(:, 5); D(:, 6)]);
end

rows = {};
for k = 1:numel(seeds)
for s = 1:height(T)
    id = T.truth_id(s);
    if ~isfile(seed_file(k, "gp", id)) || ~isfile(seed_file(k, "numeric", id)), continue, end
    g = load(seed_file(k, "gp", id), 'N', 'C', 'trange', 'meta');
    n = load(seed_file(k, "numeric", id), 'N');
    p = g.meta;
    [i, j] = find(g.C > 0);
    for h = 1:every:size(g.N, 2)
        Ni = g.N(i, h);  Nj = g.N(j, h);
        on = (Ni > 0 | Nj > 0) & isfinite(Ni) & isfinite(Nj);
        if any(on)
            X = [repmat([p.gamma1, p.lambda1, p.delta, p.epsilon], nnz(on), 1), Ni(on), Nj(on)];
            [~, sd] = predict(G.stzd_flux_gprModel, X);
            sd = sd * G.flux_std;
        else
            sd = NaN;
        end
        a = g.N(:, 1:h);  b = n.N(:, 1:h);
        rows(end+1, :) = {seeds(k), id, g.trange(h) * 180, max(g.N(:, h)) / Ntrain_max, median(sd), max(sd), ...
            nmse_eq7(a, b)}; %#ok<AGROW>
    end
end
end
P = cell2table(rows, 'VariableNames', {'seed', 'truth_id', 'day', 'maxN_over_train_max', 'sd_median', 'sd_max', 'nmse_to_day'});
writetable(P, fullfile(here, 'results', "gp_trajectory_l2f" + tag + ".csv"));

fprintf('== GP sNTM: s.d. at day 0 and its largest value, largest N / training max, final NMSE\n');
U = unique(P(:, {'seed', 'truth_id'}), 'stable');
for u = 1:height(U)
    q = P(P.seed == U.seed(u) & P.truth_id == U.truth_id(u), :);
    [sd_peak, k] = max(q.sd_max);
    fprintf('%-5s %-24s sd_max day0 %.3g, peak %.3g (day %.0f, x%.1f); max N/train max %.3g; NMSE %.3g\n', ...
        U.seed(u), U.truth_id(u), ...
        q.sd_max(1), sd_peak, q.day(k), sd_peak / q.sd_max(1), max(q.maxN_over_train_max), q.nmse_to_day(end));
end
end

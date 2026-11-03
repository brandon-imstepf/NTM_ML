function B = mass_balance_l2f()
% Tau mass balance of the network simulations of the 8 l2f truths (numeric NTM
% and each sNTM, 180 days; seeds as in network_checks_l2f), for R2.M4.
% NetworkTransportModel updates N by
%   (V_i (1 + m_i) + v1_i + v2_i) (N_i(h+1) - N_i(h)) = dt * (net flux into i),
% with m_i = dM_i/dN_i and v1 + v2 the change of the mass held on the edges of
% node i. The net fluxes sum to zero over the nodes for any flux matrix
% (sum_i [sum_k C_ki F_ki - sum_k C_ik F_ik] = 0), so node mass plus edge mass
% changes only through the first-order step in M: the drift of total mass up to
% step h is
%   sum_{h' <= h} sum_i V_i [ (N + M)_i(h'+1) - (N + M)_i(h') - (1 + m_i(h')) (N_i(h'+1) - N_i(h')) ].
%   drift_max      max over time of |drift| / initial node mass
%   node_mass_end  node mass sum_i V_i (N_i + M_i) at 180 days / initial node mass
%   node_mass_diff max over time of |node mass (sNTM) - node mass (numeric)| /
%                  initial node mass (sNTMs only)
% Writes results/mass_balance_l2f.csv.
%
%   B = mass_balance_l2f();

here = fileparts(mfilename('fullpath'));
cfg = sntm_config();
T = readtable(fullfile(here, 'manifests', 'recovery_truths.csv'), 'TextType', 'string');
T = T(T.convention == "l2f", :);
methods = ["numeric", "mlp", "sr", "gp", "lr"];
seeds = ["ECl", "CA1", "RSPd"];
seed_dirs = ["outputs_recovery", "outputs_seeds", "outputs_seeds"];
seed_file = @(k, m, id) fullfile(here, seed_dirs(k), "rec_l2f__" + m + "__" + seeds(k) + "__" + id + ".mat");

% Region volumes of the 30-region subnetwork, selected as in NetworkTransportModel
load(fullfile(cfg.matdir, 'DefaultAtlas.mat'), 'DefaultAtlas');
load(fullfile(cfg.matdir, 'CCF_labels.mat'), 'CCF_labels');
inds = ismember(CCF_labels(:,3), 'Hippocampus') | ismember(CCF_labels(:,1), 'Piriform area') | ...
    ismember(CCF_labels(:,3), 'Retrosplenial Area');
V = DefaultAtlas.volumes(inds);
V = V(:);

rows = {};
for k = 1:numel(seeds)
for s = 1:height(T)
    id = T.truth_id(s);
    if ~all(arrayfun(@(m) isfile(seed_file(k, m, id)), methods)), continue, end
    node_num = [];
    for m = methods
        r = load(seed_file(k, m, id), 'N', 'M', 'meta');
        g1 = r.meta.gamma1;
        mt = 2 * g1 * r.N / cfg.beta;            % dM/dN for gamma2 = 0
        node = sum(V .* (r.N + r.M), 1);
        step = sum(V .* (diff(r.N + r.M, 1, 2) - (1 + mt(:, 1:end-1)) .* diff(r.N, 1, 2)), 1);
        drift = cumsum(step) / node(1);
        if m == "numeric", node_num = node; end
        rows(end+1, :) = {seeds(k), id, m, max(abs(drift)), node(end) / node(1), ...
            max(abs(node - node_num)) / node_num(1)}; %#ok<AGROW>
    end
end
end
B = cell2table(rows, 'VariableNames', {'seed', 'truth_id', 'method', 'drift_max', 'node_mass_end', 'node_mass_diff'});
writetable(B, fullfile(here, 'results', 'mass_balance_l2f.csv'));

fprintf('== Mass balance, median [min, max] over %d simulations per method (fractions of initial node mass)\n', ...
    sum(B.method == "numeric"));
for m = methods
    b = B(B.method == m, :);
    fprintf('%-8s drift %s   node mass at 180 d %s   |node mass - numeric| %s\n', m, mmm(b.drift_max), ...
        mmm(b.node_mass_end), mmm(b.node_mass_diff));
end
end

function s = mmm(x)
s = sprintf('%.3g [%.3g, %.3g]', median(x, 'omitnan'), min(x), max(x));
end

function export_handoff(setting, seed, out_dir)
% Hand-off package of one network setting (all models): copies the run files
% (numeric NTM and the MLP, GP, SR and LR sNTMs, outputs_consistent/) and writes,
% per model, CSVs of regional soluble tau N and insoluble tau M (uM; rows =
% regions, columns = days), plus the Fig. 8-style and trajectory figures of that
% setting if they exist. Writes a manifest of the files with their settings.
% The nearest-neighbour interpolant sNTM ("interp", make_surrogate.m; log10
% gamma1 and lambda1, then all inputs z-scored) is included; its model, the
% 10k-row lookup table and the scaling, goes to model/.
%
%   export_handoff("g0.005_l0.02_d20_e80", "CA1", out_dir)

here = fileparts(mfilename('fullpath'));
nv = fileparts(here);
methods = ["numeric", "mlp", "gp", "sr", "lr", "interp"];
if ~isfolder(out_dir), mkdir(out_dir); end
mat_dir = fullfile(out_dir, 'mat');  csv_dir = fullfile(out_dir, 'csv');  fig_dir = fullfile(out_dir, 'figures');
for d = {mat_dir, csv_dir, fig_dir}, if ~isfolder(d{1}), mkdir(d{1}); end, end
rows = {};
for m = methods
    src = fullfile(nv, 'outputs_consistent', "rec_l2f__" + m + "__" + seed + "__" + setting + ".mat");
    [~, fname, ext] = fileparts(src);
    copyfile(src, mat_dir);
    s = load(src, 'N', 'M', 'trange', 'region_names', 'meta', 'runtime_s');
    day = s.trange(:)' * 180;
    head = ["region", compose("day_%g", day)];
    for q = ["N", "M"]
        T = [table(string(s.region_names(:)), 'VariableNames', "region"), ...
            array2table(s.(q), 'VariableNames', head(2:end))];
        writetable(T, fullfile(csv_dir, q + "_" + m + ".csv"));
    end
    rows(end+1, :) = {m, string(fname) + ext, s.meta.gamma1, s.meta.lambda1, ...
        s.meta.delta, s.meta.epsilon, lam2(s), s.meta.init_rescale, numel(day), max(day), s.runtime_s}; %#ok<AGROW>
end
writetable(cell2table(rows, 'VariableNames', {'model', 'file', 'gamma', 'lambda', 'delta', 'epsilon', ...
    'lambda2', 'init_rescale', 'n_timepoints', 'days', 'runtime_s'}), fullfile(out_dir, 'manifest.csv'));
figs = ["fig8_sets", "scatter_months1-5"; "trajectory_sets", "trajectories_days0-60"];
for k = 1:size(figs, 1)
    f = fullfile(here, 'figures', figs(k, 1), 'interior', setting + "__" + seed + ".png");
    if isfile(f), copyfile(f, fullfile(fig_dir, figs(k, 2) + "__" + setting + "__" + seed + ".png")); end
end
export_interp_model(fullfile(out_dir, 'model'));
fprintf('wrote %s\n', out_dir);
end

function v = lam2(s)
% lambda2 of the run (not stored for the interp run; surrogate runs do not depend on it)
v = NaN;  if isfield(s.meta, 'lambda2'), v = s.meta.lambda2; end
end

function export_interp_model(model_dir)
% The interp sNTM as make_surrogate.m builds it: its training samples
% (interp_training_data.m: all 10^6 samples of the MLP's training data,
% data_bias_e5_1..10.csv; too large to copy, so model/README_lookup_table.txt
% names them) and the scaling of the six inputs used for the nearest-neighbour
% distance: the inputs in cfg.interp_log_inputs on log10 scale, then each
% z-scored with the mean and s.d. written here (of the log10 values for those
% inputs).
if ~isfolder(model_dir), mkdir(model_dir); end
old = fullfile(model_dir, 'interp_lookup_table.csv');
if isfile(old), delete(old); end                     % the 10^4-row table of before 11/13
cfg = sntm_config();
T = interp_training_data(cfg);
X = T(:, 1:6);  L = cfg.interp_log_inputs;
X(:, L) = log10(X(:, L));
nm = ["gamma1"; "lambda1"; "delta"; "epsilon"; "NRow"; "NCol"];
logged = false(6, 1);  logged(L) = true;
S = table(nm, logged, mean(X)', std(X)', 'VariableNames', ["input", "log10", "mean", "std"]);
writetable(S, fullfile(model_dir, 'interp_input_scaling.csv'));
[~, f, e] = cellfun(@fileparts, cfg.interp_train_files, 'UniformOutput', false);
fid = fopen(fullfile(model_dir, 'README_lookup_table.txt'), 'w');
fprintf(fid, ['The lookup table of the nearest-neighbor sNTM is all %d training samples of the MLP,\n' ...
    'the files below (repository NTM_ML, data/data_neural_network_model/training_data_e6/), in this order.\n' ...
    'Columns: gamma1, lambda1, delta, epsilon, NRow (N1, uM), NCol (N2, uM), FValue (phi_tau), W1 (Delta_tau).\n\n'], size(T, 1));
fprintf(fid, '%s\n', string(f) + string(e));
fclose(fid);
end

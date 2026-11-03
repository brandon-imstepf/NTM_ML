function R = compare_fixedflux()
% Numeric runs with per-edge flux solves (../outputs_fixedflux, SNTM_FLUX_PER_EDGE=1)
% against the original numeric runs (joint solves), and each sNTM against both.
% NMSE = sum (s - n)^2 / sum n^2 over all regions and time points, as in
% network_checks_l2f. Also counts N values clamped from negative to eps.
% Writes ../results/fixedflux_compare.csv.
%
%   R = compare_fixedflux();

here = fileparts(mfilename('fullpath'));
nv = fileparts(here);
nmse = @(s, n) sum((s(:) - n(:)).^2) / sum(n(:).^2);
files = dir(fullfile(nv, 'outputs_fixedflux', 'rec_l2f__numeric__*.mat'));
methods = ["mlp", "sr", "gp", "lr"];
rows = {};
for f = files'
    parts = split(erase(string(f.name), ".mat"), "__");      % rec_l2f, numeric, seed, truth_id
    seed = parts(3);  id = parts(4);
    src = "outputs_recovery";  if seed ~= "ECl", src = "outputs_seeds"; end
    fix = load(fullfile(f.folder, f.name), 'N', 'runtime_s');
    old = load(fullfile(nv, src, "rec_l2f__numeric__" + seed + "__" + id + ".mat"), 'N');
    rows(end+1, :) = {seed, id, "numeric_joint", nmse(old.N, fix.N), NaN, nnz(old.N == eps), nnz(fix.N == eps), fix.runtime_s}; %#ok<AGROW>
    for m = methods
        s = load(fullfile(nv, src, "rec_l2f__" + m + "__" + seed + "__" + id + ".mat"), 'N');
        rows(end+1, :) = {seed, id, m, nmse(s.N, fix.N), nmse(s.N, old.N), nnz(s.N == eps), NaN, NaN}; %#ok<AGROW>
    end
end
R = cell2table(rows, 'VariableNames', {'seed', 'truth_id', 'method', 'nmse_vs_fixed', 'nmse_vs_original', ...
    'n_clamped', 'n_clamped_fixed', 'runtime_fixed_s'});
writetable(R, fullfile(nv, 'results', 'fixedflux_compare.csv'));
disp(R)
end

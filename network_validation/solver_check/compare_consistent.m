function R = compare_consistent()
% MLP sNTM against numeric NTM runs with each fix, seed ECl:
%   original    joint flux solves, W1 at lambda2 = 0.01   (outputs_recovery)
%   fixedflux   per-edge flux solves, W1 at lambda2 = 0.01 (outputs_fixedflux)
%   consistent  per-edge flux solves, W1 at lambda2 = lambda1, the convention of the
%               surrogates' training data (outputs_consistent)
% The MLP runs: as in the paper (W1 path at lambda2 = 0.01, outputs_recovery) and
% with the mass balance at lambda2 = lambda1 (outputs_consistent).
% NMSE = sum (s - n)^2 / sum n^2 over all regions and time points.
% Writes ../results/consistent_compare.csv.
here = fileparts(mfilename('fullpath'));
nv = fileparts(here);
nmse = @(s, n) sum((s(:) - n(:)).^2) / sum(n(:).^2);
ld = @(d, m, id) load(fullfile(nv, d, "rec_l2f__" + m + "__ECl__" + id + ".mat"), 'N');
rows = {};
for f = dir(fullfile(nv, 'outputs_consistent', 'rec_l2f__numeric__ECl__*.mat'))'
    id = extractBetween(string(f.name), "ECl__", ".mat");
    no = ld('outputs_recovery', "numeric", id);  nf = ld('outputs_fixedflux', "numeric", id);
    nc = ld('outputs_consistent', "numeric", id);
    mo = ld('outputs_recovery', "mlp", id);  mc = ld('outputs_consistent', "mlp", id);
    rows(end+1, :) = {id, nmse(nf.N, no.N), nmse(nc.N, nf.N), nmse(mo.N, no.N), nmse(mo.N, nf.N), ...
        nmse(mc.N, nc.N), nmse(mc.N, mo.N)}; %#ok<AGROW>
end
R = cell2table(rows, 'VariableNames', {'truth_id', 'fixedflux_vs_original', 'consistent_vs_fixedflux', ...
    'mlp_vs_original', 'mlp_vs_fixedflux', 'mlp_consistent_vs_consistent', 'mlp_consistent_vs_mlp'});
writetable(R, fullfile(nv, 'results', 'consistent_compare.csv'));
disp(R)
end

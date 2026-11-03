function Z = zigzag_fixedflux()
% Oscillation of the numeric runs before and after the flux fix: per region, the
% number of sign changes of the step-to-step change dN between days 18 and 180
% (a smooth curve has at most a few), and the largest |second difference| over
% the region's peak. Also the MLP sNTM for comparison.
% Writes ../results/fixedflux_zigzag.csv.
here = fileparts(mfilename('fullpath'));
nv = fileparts(here);
rows = {};
for f = dir(fullfile(nv, 'outputs_fixedflux', 'rec_l2f__numeric__*.mat'))'
    parts = split(erase(string(f.name), ".mat"), "__");
    seed = parts(3);  id = parts(4);
    src = "outputs_recovery";  if seed ~= "ECl", src = "outputs_seeds"; end
    runs = {"numeric_fixed", fullfile(f.folder, f.name); ...
            "numeric_joint", fullfile(nv, src, "rec_l2f__numeric__" + seed + "__" + id + ".mat"); ...
            "mlp", fullfile(nv, src, "rec_l2f__mlp__" + seed + "__" + id + ".mat")};
    for k = 1:size(runs, 1)
        r = load(runs{k, 2}, 'N', 'trange');
        keep = r.trange * 180 >= 18;
        d = diff(r.N(:, keep), 1, 2);
        flips = sum(diff(sign(d), 1, 2) ~= 0, 2);
        jump = max(abs(diff(r.N(:, keep), 2, 2)), [], 2) ./ max(r.N, [], 2);
        rows(end+1, :) = {seed, id, runs{k, 1}, nnz(flips > 10), max(flips), max(jump)}; %#ok<AGROW>
    end
end
Z = cell2table(rows, 'VariableNames', {'seed', 'truth_id', 'run', 'regions_over_10_flips', 'max_flips', 'max_jump_over_peak'});
writetable(Z, fullfile(nv, 'results', 'fixedflux_zigzag.csv'));
disp(Z)
end

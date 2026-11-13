function cases = fig8_cases(groups)
% The network run sets compared with the numeric NTM (plot_fig8_sets.m,
% plot_trajectory_sets.m). One struct per setting and seed:
%   group, setting (e.g. "g0.005_l0.05_d20_e80"), seed, src_fmt (run files
%   relative to network_validation/, %s = method), num_file (numeric run when it
%   is not sprintf(src_fmt, 'numeric'), else '').
% Groups:
%   "fig8"      the Fig. 8 setting (outputs_fig8/)
%   "interior"  8 settings x seeds ECl, CA1, RSPd (outputs_consistent/rec_l2f__*)
%   "edge"      8 settings on the corners of the training range, ECl
%   "ood"       8 settings outside the training range, ECl
% All numeric runs are in the surrogates' training convention (flux at
% lambda2 = 0.01, W1 at lambda2 = lambda1, per-edge flux). For edge and ood the
% numeric runs are outputs_consistent/<setting>__ECl__numeric_l2f.mat
% (manifests/consistent_edge_ood_numeric.csv) and the surrogate runs are those of
% outputs/ (surrogate runs do not depend on lambda2).
%
%   cases = fig8_cases(["interior", "edge"]);

if nargin < 1, groups = ["fig8", "interior", "edge", "ood"]; end
nv = fileparts(fileparts(mfilename('fullpath')));
cases = struct('group', {}, 'setting', {}, 'seed', {}, 'src_fmt', {}, 'num_file', {});
for g = string(groups)
    switch g
        case "fig8"
            cases(end+1) = struct('group', g, 'setting', "g5e-05_l0.01_d10_e50", 'seed', "ECl", ...
                'src_fmt', 'outputs_fig8/fig8__%s__ECl__g5e-05_l0.01_d10_e50.mat', 'num_file', ''); %#ok<AGROW>
        case "interior"
            d = dir(fullfile(nv, 'outputs_consistent', 'rec_l2f__numeric__*.mat'));
            for k = 1:numel(d)
                t = regexp(d(k).name, 'rec_l2f__numeric__(\w+?)__(g.*)\.mat', 'tokens', 'once');
                cases(end+1) = struct('group', g, 'setting', string(t{2}), 'seed', string(t{1}), ...
                    'src_fmt', char("outputs_consistent/rec_l2f__%s__" + t{1} + "__" + t{2} + ".mat"), ...
                    'num_file', ''); %#ok<AGROW>
            end
        case {"edge", "ood"}
            M = readtable(fullfile(nv, 'manifests', 'consistent_edge_ood_numeric.csv'), 'TextType', 'string');
            M = M(M.group == "consistent_" + g, :);
            for k = 1:height(M)
                sid = extractBefore(M.setting_id(k), "__" + M.seed(k));      % e.g. edge__g5e-05_...
                cases(end+1) = struct('group', g, 'setting', extractAfter(sid, g + "__"), 'seed', M.seed(k), ...
                    'src_fmt', char("outputs/" + sid + "__" + M.seed(k) + "__%s.mat"), ...
                    'num_file', char("outputs_consistent/" + M.case_id(k) + ".mat")); %#ok<AGROW>
            end
    end
end
end

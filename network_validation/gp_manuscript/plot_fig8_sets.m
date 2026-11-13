function S = plot_fig8_sets(groups)
% Fig. 8-style figures (regional tau of each sNTM against the numeric NTM, months
% 1-5, R^2 and NMSE per panel; plot_fig8_gp.m) for every network setting:
%   "fig8"      the Fig. 8 setting (outputs_fig8/)
%   "interior"  8 settings x seeds ECl, CA1, RSPd (outputs_consistent/rec_l2f__*)
%   "edge"      8 settings on the corners of the training range, ECl
%   "ood"       8 settings outside the training range, ECl
% (run files and conventions: fig8_cases.m). Settings whose numeric run does not
% exist yet are skipped.
% Figures: figures/fig8_sets/<group>/<setting>.png/.pdf, and one PDF per group
% (figures/fig8_sets/<group>.pdf). Table: results/fig8_sets.csv (group, setting,
% seed, parameters, method, month, R2, NMSE, slope).
%
%   S = plot_fig8_sets();                      % all groups
%   S = plot_fig8_sets(["edge", "ood"]);

if nargin < 1, groups = ["fig8", "interior", "edge", "ood"]; end
here = fileparts(mfilename('fullpath'));
nv = fileparts(here);
fig_root = fullfile(here, 'figures', 'fig8_sets');
csv_file = fullfile(nv, 'results', 'fig8_sets.csv');

cases = fig8_cases(groups);   % run sets, see fig8_cases.m

S = table();
for k = 1:numel(cases)
    c = cases(k);
    num_path = c.num_file;
    if isempty(num_path), num_path = sprintf(c.src_fmt, 'numeric'); end
    if ~isfile(fullfile(nv, num_path))
        fprintf('%s %s %s: numeric run not available yet, skipped\n', c.group, c.setting, c.seed);
        continue
    end
    v = sscanf(char(c.setting), 'g%f_l%f_d%f_e%f')';
    ttl = sprintf('%s, %s seed: \\gamma = %g, \\lambda = %g, \\delta = %g, \\epsilon = %g', ...
        c.group, c.seed, v(1), v(2), v(3), v(4));
    out_dir = fullfile(fig_root, c.group);
    if ~isfolder(out_dir), mkdir(out_dir); end
    fprintf('\n=== %s %s %s ===\n', c.group, c.setting, c.seed);
    T = plot_fig8_gp(fullfile(out_dir, c.setting + "__" + c.seed), c.src_fmt, ttl, c.num_file);
    close all
    n = height(T);
    T = [table(repmat(c.group, n, 1), repmat(c.setting, n, 1), repmat(c.seed, n, 1), ...
        repmat(v(1), n, 1), repmat(v(2), n, 1), repmat(v(3), n, 1), repmat(v(4), n, 1), ...
        'VariableNames', ["group", "setting", "seed", "gamma", "lambda", "delta", "epsilon"]), T]; %#ok<AGROW>
    S = [S; T]; %#ok<AGROW>
end

% Keep rows of groups not redrawn in this call
if isfile(csv_file) && ~isempty(S)
    old = readtable(csv_file, 'TextType', 'string');
    old = old(~ismember(old.group, unique(S.group)), :);
    S = [old; S];
end
if ~isempty(S), writetable(S, csv_file); fprintf('wrote %s\n', csv_file); end

% One PDF per group, pages in setting order
for g = string(groups)
    f = dir(fullfile(fig_root, g, '*.png'));
    if isempty(f), continue; end
    out_pdf = fullfile(fig_root, g + ".pdf");
    if isfile(out_pdf), delete(out_pdf); end
    for k = 1:numel(f)
        im = imread(fullfile(f(k).folder, f(k).name));
        fig = figure('Visible', 'off', 'Units', 'pixels', 'Position', [0 0 size(im, 2) size(im, 1)] / 3);
        ax = axes(fig, 'Position', [0 0 1 1]);
        image(ax, im);  axis(ax, 'image', 'off');
        exportgraphics(fig, out_pdf, 'Append', true, 'Resolution', 200);
        close(fig);
    end
    fprintf('wrote %s (%d pages)\n', out_pdf, numel(f));
end
end

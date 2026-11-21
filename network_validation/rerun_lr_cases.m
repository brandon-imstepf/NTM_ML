function rerun_lr_cases(max_cases, dirs)
% Rerun every LR network case with the refitted LR surrogate (11/19:
% fit_lr_surrogate.py now selects and fits on z-scored features; the raw-input
% fit had picked gamma*N2 and epsilon*N2 in place of epsilon for phi_tau).
% For each LR run file in the output folders below, the case is rerun from its
% manifest row with the settings stored in the file (cfg.lambda2_w1_is_lambda1,
% cfg.flux_per_edge; files from before these options have neither, which is the
% default), into outputs_lrfix/<folder>/. Nothing in the original folders is
% changed; swap_lr_cases.m moves the old files aside and the new ones in.
%
%   rerun_lr_cases()          all
%   rerun_lr_cases(1)         first case only (timing test)
%   rerun_lr_cases(Inf, "outputs_fig8")

here = fileparts(mfilename('fullpath'));
src = struct('dir', {"outputs_consistent", "outputs_seeds", "outputs", "outputs_fig8"}, ...
    'manifests', {"consistent_l2f_surrogate", "seeds_l2f_surrogate", ["surrogate", "surrogate_interior"], "fig8_gp"});
if nargin < 1, max_cases = Inf; end
if nargin >= 2, src = src(ismember([src.dir], string(dirs))); end
n = 0;
for s = src
    d = dir(fullfile(here, s.dir, '*.mat'));
    for k = 1:numel(d)
        f = fullfile(d(k).folder, d(k).name);
        info = whos('-file', f);
        if ~ismember('meta', {info.name}), continue, end
        r = load(f, 'meta', 'cfg');
        if r.meta.method ~= "lr", continue, end
        % the manifest row of this case
        row = [];  mf = "";
        for m = s.manifests
            M = readtable(fullfile(here, 'manifests', m + ".csv"), 'TextType', 'string', 'Delimiter', ',');
            row = find(M.case_id == r.meta.case_id);
            if ~isempty(row), mf = fullfile(here, 'manifests', m + ".csv"); break, end
        end
        if isempty(row)        % not in the listed manifests: search them all
            for m = string({dir(fullfile(here, 'manifests', '*.csv')).name})
                M = readtable(fullfile(here, 'manifests', m), 'TextType', 'string', 'Delimiter', ',');
                if ~ismember('case_id', M.Properties.VariableNames), continue, end
                row = find(M.case_id == r.meta.case_id);
                if ~isempty(row), mf = fullfile(here, 'manifests', m); break, end
            end
        end
        assert(~isempty(row), 'no manifest row for %s', r.meta.case_id);
        out = fullfile(here, 'outputs_lrfix', s.dir);
        setenv('SNTM_OUT_DIR', out);
        l2 = isfield(r.cfg, 'lambda2_w1_is_lambda1') && r.cfg.lambda2_w1_is_lambda1;
        fe = isfield(r.cfg, 'flux_per_edge') && r.cfg.flux_per_edge;
        if l2, setenv('SNTM_LAMBDA2_W1', 'lambda1'); else, setenv('SNTM_LAMBDA2_W1', ''); end
        if fe, setenv('SNTM_FLUX_PER_EDGE', '1'); else, setenv('SNTM_FLUX_PER_EDGE', ''); end
        fprintf('%s/%s (row %d of %s; lambda2_w1=lambda1: %d)\n', s.dir, r.meta.case_id, row, mf, l2);
        t0 = tic;
        run_sntm_cases(mf, row, row);
        fprintf('  %.1f s\n', toc(t0));
        n = n + 1;
        if n >= max_cases, return, end
    end
end
setenv('SNTM_OUT_DIR', '');  setenv('SNTM_LAMBDA2_W1', '');  setenv('SNTM_FLUX_PER_EDGE', '');
end

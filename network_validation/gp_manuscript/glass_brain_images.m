function glass_brain_images(days, methods, out_dir, src_fmt)
% Glass-brain images for Fig. 7a of GP_manuscript.tex: tau (N + M) at the
% Fig. 8 parameters (gamma1 5e-5, lambda 0.01, delta 10, epsilon 50, seed ECl,
% the Table 3 timing runs, repeat 1) for each method and day, drawn with
% Brainframe (https://github.com/Raj-Lab-UCSF/Brainframe; set BRAINFRAME_ROOT,
% default ~/Documents/MATLAB/Brainframe) as in BrainframePlot.m of the
% simulator: N + M on the 30 regions, mapped to the 426 CCF regions, divided by
% the largest value of the numeric run, values below a quarter of the median
% nonzero value set to 0, square root; red point clouds, dorsal view.
% Writes <out_dir>/bf_<method>_d<day>.png.
%
%   glass_brain_images();                                  % all methods, days 0 7 15 30 60
%   glass_brain_images(30, "mlp", 'scratch');
%   glass_brain_images([], [], out_dir, 'outputs_fig8/fig8__%s__ECl__g5e-05_l0.01_d10_e50.mat')
%       other runs: src_fmt relative to network_validation/, %s = method
%       (make_glass_brain_figure.m)
% Day 0 is drawn for the numeric run only (all methods share the initial state).

here = fileparts(mfilename('fullpath'));
nv = fileparts(here);
if nargin < 1 || isempty(days), days = [0 7 15 30 60]; end
if nargin < 2 || isempty(methods), methods = ["numeric", "mlp", "gp", "sr", "lr"]; end
if nargin < 3, out_dir = fullfile(here, 'figures', 'glass_brain'); end
if ~isfolder(out_dir), mkdir(out_dir); end
addpath(nv);
cfg = sntm_config();
addpath(fullfile(cfg.sim_root, 'Code', 'Transport_Network'));
bf_root = getenv('BRAINFRAME_ROOT');
if isempty(bf_root), bf_root = fullfile(getenv('USERPROFILE'), 'Documents', 'MATLAB', 'Brainframe'); end
addpath(bf_root);
% Brainframe draws with pointCloud/pcshow (Computer Vision Toolbox); without
% that toolbox, use the scatter3 stand-ins in bf_shims/
if isempty(which('pointCloud')), addpath(fullfile(here, 'bf_shims')); end

if nargin < 4 || isempty(src_fmt)
    src_fmt = 'outputs_timing/timing__%s__ECl__g5e-05_l0.01_d10_e50__rep1.mat';
end
file = @(m) fullfile(nv, sprintf(src_fmt, m));
ref = load(file("numeric"), 'N', 'M', 'trange', 'region_names');
X_ref = DataToCCF_Transport(ref.N + ref.M, ref.region_names, 1, cfg.matdir);
scale = max(X_ref(:), [], 'omitnan');

for m = string(methods)
    r = load(file(m), 'N', 'M', 'trange', 'region_names');
    for d = days
        if d == 0 && m ~= "numeric", continue; end
        [~, h] = min(abs(r.trange * 180 - d));
        X = DataToCCF_Transport(r.N(:, h) + r.M(:, h), r.region_names, 1, cfg.matdir);
        X(isnan(X)) = 0;
        X = X / scale;
        X(X < 0.25 * median(nonzeros(X))) = 0;
        X = sqrt(X);
        label = sprintf('bf_%s_d%d', m, d);
        input_struct = brainframe_inputs_mouse(bf_root, 'region_groups', ones(426, 1), ...
            'cmap', [1 0 0], 'voxUreg', 1, 'sphere', 0, 'centered', [0 1], ...
            'data', X, 'norm_method', 'none', 'xfac', 3, 'pointsize', 3, ...
            'bgcolor', 'w', 'img_labels', label, 'img_format', 'png', ...
            'img_directory', char(out_dir), 'img_renderer', 'opengl', ...
            'img_views', [0 0 1], 'savenclose', 1);
        brainframe(input_struct);
        fprintf('%s\n', label);
    end
end
end

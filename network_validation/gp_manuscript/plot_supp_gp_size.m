function plot_supp_gp_size(out_file)
% Supplementary figure (11/10): held-out NMSE (Eq. 7) of the GP against the
% number of training samples, for phi_tau and Delta_tau, from
% results/supp_gp_griddatan.csv (supp_gp_griddatan.m). Circles: the GP
% configuration of the paper trained on n rows of data_bias_e5_1.csv (n = 5000
% is the paper's GP); open square (drawn at n = 1150 to keep it apart): trained on bias_train_e3.csv, the 10^3-sample
% training set of linear and symbolic regression. Error bars: s.d. over 10^3
% bootstrap resamples of 10^3 held-out samples. Lines: the other surrogates
% and the best griddatan variant on the same 10^5 held-out samples
% (results/domain_edge_error.csv, bin "all"; griddatan from supp_gp_griddatan.csv).
%
%   plot_supp_gp_size();    % gp_manuscript/figures/supp_gp_size.png

here = fileparts(mfilename('fullpath'));
nv = fileparts(here);
if nargin < 1, out_file = fullfile(here, 'figures', 'supp_gp_size'); end
T = readtable(fullfile(nv, 'results', 'supp_gp_griddatan.csv'), 'TextType', 'string');
E = readtable(fullfile(nv, 'results', 'domain_edge_error.csv'), 'TextType', 'string');
E = E(E.bin == "all", :);
gp = T(startsWith(T.model, "gp_") & T.model ~= "gp_lrsr_1000", :);
gp = sortrows(gp, 'n_train');
lrsr = T(T.model == "gp_lrsr_1000", :);
gd = T(startsWith(T.model, "griddatan_"), :);
target = ["flux", "emc"];  tlab = ["\phi_\tau", "\Delta_\tau"];  etarget = ["flux", "w1"];
refs = ["mlp", "sr", "lr"];  rlab = ["MLP", "symbolic regression", "linear regression"];
rcol = [0.17 0.40 0.69; 0.85 0.45 0.10; 0.55 0.55 0.55];
fig = figure('Units', 'centimeters', 'Position', [2 2 18 8], 'Color', 'w');
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for t = 1:2
    ax = nexttile(tl);  hold(ax, 'on');  box(ax, 'on');  grid(ax, 'on');
    h = gobjects(0);  lab = {};
    for r = 1:numel(refs)
        v = E.nmse_mean(E.method == refs(r) & E.target == etarget(t));
        h(end+1) = yline(ax, v, '--', 'Color', rcol(r, :), 'LineWidth', 1.2); %#ok<AGROW>
        lab{end+1} = rlab(r); %#ok<AGROW>
    end
    [vg, kg] = min(gd.("heldout_" + target(t) + "_nmse_mean"));
    h(end+1) = yline(ax, vg, ':', 'Color', 'k', 'LineWidth', 1.2);
    lab{end+1} = "griddatan, best (" + strrep(extractAfter(gd.model(kg), "griddatan_"), "_", " ") + ")";
    m = gp.("heldout_" + target(t) + "_nmse_mean");  s = gp.("heldout_" + target(t) + "_nmse_sd");
    h(end+1) = errorbar(ax, gp.n_train, m, s, '-o', 'Color', [0.11 0.62 0.47], 'MarkerFaceColor', [0.11 0.62 0.47], ...
        'LineWidth', 1.4, 'MarkerSize', 6);
    lab{end+1} = 'GP, training rows of the MLP data';
    h(end+1) = errorbar(ax, lrsr.n_train + 150, lrsr.("heldout_" + target(t) + "_nmse_mean"), ...
        lrsr.("heldout_" + target(t) + "_nmse_sd"), 's', 'Color', [0.11 0.62 0.47], 'MarkerFaceColor', 'w', 'MarkerSize', 8, 'LineWidth', 1.4);
    lab{end+1} = 'GP, LR/SR training set (10^3)';
    set(ax, 'YScale', 'log', 'XLim', [0 5500], 'FontSize', 9, 'TickDir', 'out');
    xlabel(ax, 'GP training samples');
    ylabel(ax, "held-out NMSE, " + tlab(t));
    title(ax, tlab(t), 'FontWeight', 'normal');
end
lg = legend(h, lab, 'Orientation', 'horizontal', 'NumColumns', 3, 'FontSize', 8);  lg.Layout.Tile = 'south';
exportgraphics(fig, out_file + ".png", 'Resolution', 300);
exportgraphics(fig, out_file + ".pdf", 'ContentType', 'vector');
fprintf('wrote %s.png/.pdf\n', out_file);
end

function S = plot_supp_error_fig(out_dir)
% Supplementary Fig. 2 (with the GP, 11/10; 1-NN, 11/16): single-edge prediction error of
% each selected surrogate against the true value, on its training and test
% data, as in Table 1 of the paper.
%   LR (feature selected) and SR: trained on bias_train_e3.csv, tested on
%       bias_val_e3.csv (10^3 each)
%   GP: its 5,000 training rows (data_bias_e5_1.csv, rng(42) randperm), tested on
%       data_bias_e5_11.csv (10^5 held-out; the GP was selected on bias_val_e3)
%   1-NN interpolation: data_bias_e5_1.csv (10^5 of its 10^6 lookup-table rows;
%       error zero by construction), tested on data_bias_e5_11.csv (10^5 held-out)
%   MLP: data_bias_e5_1.csv (10^5 of its 10^6 training rows), tested on
%       data_bias_e5_11.csv (10^5 held-out)
% (a) absolute error, predicted - true, against the true value;
% (b) relative error |predicted - true| / |true| against |true| (log axes;
%     true values of 0 are left out).
% Rows: surrogates; columns: phi_tau train, phi_tau test, Delta_tau train,
% Delta_tau test. Each panel title gives the NMSE (Eq. 7) over the whole set.
% Axes are shared down each column (x; also y for the relative error).
% Surrogates from make_surrogate.m (the models of the network simulations).
% Font sizes 25% larger than before (11/16).
% Writes <out_dir>/supp_error_abs.png, supp_error_rel.png (and .pdf), and
% results/supp_error_fig_nmse.csv.
%
%   plot_supp_error_fig();      % gp_manuscript/figures/

here = fileparts(mfilename('fullpath'));
nv = fileparts(here);
addpath(nv);
if nargin < 1, out_dir = fullfile(here, 'figures'); end
if ~isfolder(out_dir), mkdir(out_dir); end
cfg = sntm_config();
root = cfg.ntm_ml_root;
e3 = fullfile(root, 'data', 'data_baseline_models', 'e3_data');
e6 = fullfile(root, 'data', 'data_neural_network_model', 'training_data_e6');
lrsr_train = readmatrix(fullfile(e3, 'bias_train_e3.csv'));
lrsr_test = readmatrix(fullfile(e3, 'bias_val_e3.csv'));
super = readmatrix(fullfile(e6, 'data_bias_e5_1.csv'));
heldout = readmatrix(fullfile(e6, 'data_bias_e5_11.csv'));
G = load(cfg.gp_model_file);
rng(G.training_seed);
gp_train = super(randperm(1e5, G.dataset_size), :);

% rows from the bottom up in the paper order (11/17): LR (bottom), SR, 1-NN, GP, MLP (top)
methods = ["mlp", "gp", "interp", "sr", "lr"];   % top to bottom
labels = ["MLP", "Gaussian process", "1-NN interpolation", "Symbolic regression", "Linear regression"];   % (11/20) "(feature selected)" dropped: too long for the row at 25% larger text; in the caption
sets = struct('lr', {{lrsr_train, lrsr_test}}, 'sr', {{lrsr_train, lrsr_test}}, ...
    'gp', {{gp_train, heldout}}, 'interp', {{super, heldout}}, 'mlp', {{super, heldout}});
P = struct();
rows = {};
for m = methods
    sur = make_surrogate(char(m), cfg);
    for k = 1:2
        D = sets.(m){k};
        p = [chunked(sur.predict_flux, D(:, 1:6)), chunked(sur.predict_w1, D(:, 1:6))];
        P.(m){k} = struct('p', p, 'y', D(:, 7:8));
        rows(end+1, :) = {m, k, size(D, 1), nmse_eq7(p(:, 1), D(:, 7)), nmse_eq7(p(:, 2), D(:, 8))}; %#ok<AGROW>
    end
end
S = cell2table(rows, 'VariableNames', {'method', 'set', 'n', 'flux_nmse', 'emc_nmse'});
S.set = categorical(S.set, [1 2], {'train', 'test'});
writetable(S, fullfile(nv, 'results', 'supp_error_fig_nmse.csv'));
disp(S)

tname = ["\phi_\tau", "\Delta_\tau"];
setname = ["train", "test"];
col = [0.13 0.40 0.67];
fs = 1.25 * 1.25;   % (11/20, user) all text 25% larger than before, then another 25% (x1.5625)
for kind = ["abs", "rel"]
    fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [1 1 30 32], 'Color', 'w');
    tl = tiledlayout(fig, numel(methods), 4, 'TileSpacing', 'compact', 'Padding', 'compact');
    axs = gobjects(numel(methods), 4);
    for r = 1:numel(methods)
        m = methods(r);
        for t = 1:2
            for k = 1:2
                ax = nexttile(tl);  hold(ax, 'on');  box(ax, 'on');
                axs(r, 2 * (t - 1) + k) = ax;
                y = P.(m){k}.y(:, t);  p = P.(m){k}.p(:, t);
                n = nmse_eq7(p, y);
                if kind == "abs"
                    scatter(ax, y, p - y, 3, col, 'filled', 'MarkerFaceAlpha', 0.3, 'MarkerEdgeColor', 'none');
                    yline(ax, 0, '-', 'Color', [0.2 0.2 0.2]);
                    if r == numel(methods), xlabel(ax, "true " + tname(t)); end
                    if t == 1 && k == 1, ylabel(ax, {char(labels(r)), 'predicted - true'}, 'FontSize', 10 * fs); end
                else
                    keep = y ~= 0;
                    scatter(ax, abs(y(keep)), abs(p(keep) - y(keep)) ./ abs(y(keep)), 3, col, 'filled', ...
                        'MarkerFaceAlpha', 0.3, 'MarkerEdgeColor', 'none');
                    set(ax, 'XScale', 'log', 'YScale', 'log', 'YTick', 10.^(-30:10:30));   % (11/20) fixed: at the larger text MATLAB dropped labels
                    yline(ax, 1, ':', 'Color', [0.3 0.3 0.3]);
                    if r == numel(methods), xlabel(ax, "|true " + tname(t) + "|"); end
                    if t == 1 && k == 1, ylabel(ax, {char(labels(r)), 'relative error'}, 'FontSize', 10 * fs); end
                end
                if all(p == y)   % 1-NN on its own lookup table
                    text(ax, 0.5, 0.78, {'zero error', '(lookup table)'}, 'Units', 'normalized', ...
                        'HorizontalAlignment', 'center', 'FontSize', 10 * fs, 'Color', [0.3 0.3 0.3]);
                end
                ttl = sprintf('%s, %s: NMSE %.2g', tname(t), setname(k), n);
                % (11/20) in the absolute-error figure the axis exponent (x10^-n) sits
                % at the top left of each panel: an empty second line lifts the title above it
                if kind == "abs", ttl = {ttl, ''}; end
                title(ax, ttl, 'FontWeight', 'normal', 'FontSize', 10 * fs);
                set(ax, 'FontSize', 8.75 * fs, 'TickDir', 'out');
            end
        end
    end
    % same axes down each column (one target and set): x always, y for the relative error
    for c = 1:4
        if kind == "rel", linkaxes(axs(:, c), 'xy'); else, linkaxes(axs(:, c), 'x'); end
    end
    f = fullfile(out_dir, "supp_error_" + kind);
    exportgraphics(fig, f + ".png", 'Resolution', 250);
    close(fig);
    fprintf('wrote %s.png\n', f);
end
end

function y = chunked(f, X)
y = zeros(size(X, 1), 1);
for a = 1:20000:size(X, 1)
    b = min(a + 19999, size(X, 1));
    y(a:b) = f(X(a:b, :));
end
end

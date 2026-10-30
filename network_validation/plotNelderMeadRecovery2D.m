function plotNelderMeadRecovery2D(nm_file, outdir, label, truth_file)
% Figure for a Nelder-Mead network recovery of two parameters searched together
% (an nm_*.mat file of nelder_mead_network_recovery.m with two parameters):
%   nm2d_<label>.png  the path of the search in the plane of the two parameters,
%                     each parameter and the loss against simulation number
% truth_file (optional) is the truth run the search was fitted to; by default
% the path saved in the result (result.truth_file), used for the title only.
%
%   plotNelderMeadRecovery2D('outputs_nm/nm_l2f__mlp__g0.0005_l0.02_d80_e20__gamma1+delta.mat', 'outputs_nm', 'mlp_truth')

    if ~isfolder(outdir), mkdir(outdir); end
    S = load(nm_file, 'trace', 'result');
    tr = S.trace;  res = S.result;
    assert(numel(res.param) == 2, 'plotNelderMeadRecovery2D is for two parameters');
    if nargin < 4 || isempty(truth_file), truth_file = res.truth_file; end
    truth = load(truth_file, 'meta', 'trange');

    labels = containers.Map(["gamma1", "lambda1", "delta", "epsilon"], ["\gamma_1", "\lambda_1", "\delta", "\epsilon"]);
    p_lab = arrayfun(@(p) string(labels(p)), res.param);
    V = vertcat(tr.value);  L = [tr.nmse];  n = numel(L);
    [~, b] = min(L);
    truth_name = string(truth.meta.method);
    if truth_name == "mlp", truth_name = "MLP"; end
    fig_title = sprintf('Nelder-Mead recovery of %s and %s together, MLP vs %s truth, %g days', ...
        p_lab(1), p_lab(2), truth_name, round(truth.trange(end) * 180, 1));
    blue = [0.15 0.3 0.6];  red = [0.8 0.1 0.1];  path_grey = [0.6 0.6 0.6];

    fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 30 24]);
    tl = tiledlayout(fig, 3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, fig_title, 'FontSize', 13);

    % path in the plane of the two parameters
    ax = nexttile(tl, 1, [3 1]); hold(ax, 'on');
    plot(ax, V(:,1), V(:,2), '-', 'Color', [path_grey 0.7], 'LineWidth', 0.75);
    h = [scatter(ax, V(:,1), V(:,2), 60, 1:n, 'filled', 'MarkerEdgeColor', 'k'), ...
         plot(ax, res.true_value(1), res.true_value(2), 'p', 'MarkerSize', 18, 'Color', 'k', 'MarkerFaceColor', 'y'), ...
         plot(ax, V(b,1), V(b,2), 'o', 'MarkerSize', 15, 'Color', red, 'LineWidth', 2)];
    text(ax, V(1,1), V(1,2), '  start', 'FontSize', 10, 'VerticalAlignment', 'top');
    colormap(ax, sequentialMap());  clim(ax, [1 max(n, 2)]);
    cb = colorbar(ax, 'southoutside');
    xlabel(cb, 'simulation number (order of the search)');
    if res.settings.log(1), set(ax, 'XScale', 'log'); end
    if res.settings.log(2), set(ax, 'YScale', 'log'); end
    set(ax, 'Box', 'off', 'TickDir', 'out', 'FontSize', 11);
    grid(ax, 'on');
    xlabel(ax, p_lab(1));  ylabel(ax, p_lab(2));
    title(ax, sprintf('true (%.4g, %.4g), start (%.4g, %.4g), recovered (%.4g, %.4g)', ...
        res.true_value, res.start_value, V(b,:)), 'FontWeight', 'normal');
    legend(ax, h, {'simulation of the search (colour = order)', 'true values', 'recovered (lowest loss)'}, ...
        'Location', 'best', 'Box', 'off', 'FontSize', 10);
    hold(ax, 'off');

    % each parameter against simulation number
    for k = 1:2
        ax = nexttile(tl, 2 * k); hold(ax, 'on');
        plot(ax, 1:n, V(:,k), '-o', 'Color', blue, 'MarkerFaceColor', blue, 'LineWidth', 1.5);
        yline(ax, res.true_value(k), '--k', 'true', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left');
        if res.settings.log(k), set(ax, 'YScale', 'log'); end
        set(ax, 'Box', 'off', 'TickDir', 'out', 'FontSize', 11);
        grid(ax, 'on');
        ylabel(ax, p_lab(k));
        hold(ax, 'off');
    end

    % loss against simulation number
    ax = nexttile(tl, 6); hold(ax, 'on');
    h = [plot(ax, 1:n, L, 'o', 'Color', blue, 'MarkerFaceColor', blue), ...
         stairs(ax, 1:n, cummin(L), '-k', 'LineWidth', 1.5)];
    set(ax, 'YScale', 'log', 'Box', 'off', 'TickDir', 'out', 'FontSize', 11);
    grid(ax, 'on');
    xlabel(ax, 'simulation number');  ylabel(ax, 'loss: NMSE of node \tau');
    legend(ax, h, {'loss of that simulation', 'lowest loss so far'}, 'Location', 'northeast', 'Box', 'off', 'FontSize', 10);
    hold(ax, 'off');

    exportgraphics(fig, fullfile(outdir, sprintf('nm2d_%s.png', label)), 'Resolution', 150);
    close(fig);
end

function cmap = sequentialMap()
% one hue, light (early simulations) to dark (late), as in plotNelderMeadRecovery.m
    light = [0.78 0.86 0.97];
    dark = [0.05 0.17 0.45];
    t = linspace(0, 1, 256)';
    cmap = light + t .* (dark - light);
end

function recoveryOverviewPlot(rc, s, outdir)
% Two figures for locked set s, every model side by side (columns):
%   by_set/single_<set>.png  one row per parameter: loss vs that parameter, the
%                            others at their true values
%   by_set/joint_<set>.png   one row per pair: loss surface over the pair, the
%                            others at their true values; all panels on one
%                            colour scale so models and pairs can be compared
% Flat (not identifiable) recoveries have no recovered marker. A pair flat along
% both parameters is left blank.
% rc is the struct saved by synthetic_parameter_recovery_all_models.m in
% recovery_curves.mat. To redraw without rerunning the recovery:
%     load(fullfile('recovery_all_models', 'flux', 'recovery_curves.mat'));
%     for s = 1:numel(rc.set_ids), recoveryOverviewPlot(rc, s, fullfile('recovery_all_models', 'flux')); end
% rc can also be a cell of such structs, one per output (e.g. {flux rc, w1 rc}):
% each model then gets one column per output, next to each other
% (recoverySideBySidePlot.m).

    if ~iscell(rc), rc = {rc}; end
    n_out = numel(rc);
    rc1 = rc{1};   % models, parameters, sets and true values are the same in every output
    n_models = numel(rc1.model_labels);
    n_params = numel(rc1.params);
    n_pairs = size(rc1.pairs, 1);
    n_cols = n_models * n_out;
    theta_true = rc1.true_params(s,:);   % [gamma1, lambda1, delta, epsilon]
    set_title = sprintf('%s:  \\gamma_1 = %.3g,  \\lambda_1 = %.3g,  \\delta = %.3g,  \\epsilon = %.3g', ...
        strrep(rc1.set_ids(s), '_', '\_'), theta_true);
    if n_out == 1 && isfield(rc1, 'output_label')
        set_title = sprintf('%s  [loss on %s]', set_title, rc1.output_label);
    end
    set_outdir = fullfile(outdir, 'by_set');
    if ~exist(set_outdir, 'dir'), mkdir(set_outdir); end

    %% single-parameter figure
    fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 10*n_cols 8.5*n_params]);
    tl = tiledlayout(fig, n_params, n_cols, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, set_title + "  (single-parameter)", 'FontSize', 14);
    for i = 1:n_params
        p = rc1.params(i);
        for m = 1:n_models
            for o = 1:n_out
                r = rc{o};
                ax = nexttile(tl, (i-1)*n_cols + (m-1)*n_out + o); hold(ax, 'on');
                plot(ax, p.curve_grid, r.single_curves(:,i,m,s), 'Color', hex2rgb(r.model_colors(m)), 'LineWidth', 2);
                xline(ax, theta_true(i), '--k', 'LineWidth', 1);
                if r.single_flat(i,m,s)
                    title(ax, sprintf('%s:  %s flat (not identifiable)', panelName(r, m, n_out), p.label), 'FontWeight', 'normal');
                else
                    plot(ax, r.single_hat(i,m,s), r.single_loss(i,m,s), 'o', 'MarkerFaceColor', 'r', ...
                        'MarkerEdgeColor', 'k', 'MarkerSize', 7);
                    title(ax, sprintf('%s:  %s = %.3g', panelName(r, m, n_out), p.label, r.single_hat(i,m,s)), 'FontWeight', 'normal');
                end
                if p.log, set(ax, 'XScale', 'log'); end
                set(ax, 'YScale', 'log', 'Box', 'off', 'TickDir', 'out');
                xlim(ax, p.range);
                grid(ax, 'on');
                xlabel(ax, p.label);
                if m == 1 && o == 1, ylabel(ax, 'loss'); end
                hold(ax, 'off');
            end
        end
    end
    ax1 = nexttile(tl, 1); hold(ax1, 'on');
    h = [plot(ax1, NaN, NaN, '--k'), ...
         plot(ax1, NaN, NaN, 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k')];
    lg = legend(h, {'true value', 'recovered (minimum loss)'}, 'Orientation', 'horizontal', 'Box', 'off', 'FontSize', 11);
    lg.Layout.Tile = 'south';
    hold(ax1, 'off');
    exportgraphics(fig, fullfile(set_outdir, sprintf('single_%s.png', rc1.set_ids(s))), 'Resolution', 130);
    close(fig);

    %% joint figure
    surfaces = cellfun(@(r) log10(r.joint_surfaces(:,:,:,:,s)), rc, 'UniformOutput', false);
    all_values = cell2mat(cellfun(@(x) x(:), surfaces(:), 'UniformOutput', false));
    clim_all = [min(all_values), max(all_values)];
    if clim_all(1) == clim_all(2), clim_all = clim_all + [-1 1]; end
    % contour levels, the same in every panel so a shade means the same loss everywhere
    n_levels = 20;
    levels = linspace(clim_all(1), clim_all(2), n_levels + 1);

    fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 10*n_cols 8.5*n_pairs]);
    tl = tiledlayout(fig, n_pairs, n_cols, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, set_title + "  (joint)", 'FontSize', 14);
    for k = 1:n_pairs
        a = rc1.pairs(k,1);  pa = rc1.params(a);   % y axis
        b = rc1.pairs(k,2);  pb = rc1.params(b);   % x axis
        [to_x, x_label] = paramAxis(pb);
        [to_y, y_label] = paramAxis(pa);
        for m = 1:n_models
            for o = 1:n_out
                r = rc{o};
                ax = nexttile(tl, (k-1)*n_cols + (m-1)*n_out + o); hold(ax, 'on');
                flat = r.joint_flat(k,:,m,s);
                if all(flat)
                    % flat along both: nothing to show
                    axis(ax, 'off');
                    title(ax, sprintf('%s:  %s, %s flat (not identifiable)', panelName(r, m, n_out), pa.label, pb.label), ...
                        'FontWeight', 'normal');
                    hold(ax, 'off');
                    continue
                end
                % filled contours of log10 loss (shading), level lines drawn thin
                [X, Y] = meshgrid(to_x(pb.joint_grid), to_y(pa.joint_grid));
                contourf(ax, X, Y, surfaces{o}(:,:,k,m), levels, 'LineColor', [0.25 0.25 0.25], 'LineWidth', 0.3);
                img_ax = ax; % any surface panel carries the shared colour scale
                set(ax, 'Box', 'off', 'TickDir', 'out', 'CLim', clim_all, 'Layer', 'top');
                axis(ax, 'tight');
                plot(ax, to_x(theta_true(b)), to_y(theta_true(a)), 'x', 'Color', 'm', 'MarkerSize', 12, 'LineWidth', 2.5);
                hat = r.joint_hat(k,:,m,s);
                if any(flat)
                    % one of the two is not identifiable: report the other only
                    ok = find(~flat);
                    title(ax, sprintf('%s:  %s = %.3g,  %s flat', panelName(r, m, n_out), ...
                        rc1.params(rc1.pairs(k,ok)).label, hat(ok), rc1.params(rc1.pairs(k,3-ok)).label), 'FontWeight', 'normal');
                else
                    plot(ax, to_x(hat(2)), to_y(hat(1)), 'o', 'MarkerFaceColor', 'r', ...
                        'MarkerEdgeColor', 'k', 'MarkerSize', 9, 'LineWidth', 1.5);
                    title(ax, sprintf('%s:  %s = %.3g,  %s = %.3g', panelName(r, m, n_out), ...
                        pa.label, hat(1), pb.label, hat(2)), 'FontWeight', 'normal');
                end
                xlabel(ax, x_label);
                if m == 1 && o == 1, ylabel(ax, y_label); end
                hold(ax, 'off');
            end
        end
    end
    cb = colorbar(img_ax);
    cb.Layout.Tile = 'east';
    ylabel(cb, 'log_{10} loss');
    ax1 = nexttile(tl, 1); hold(ax1, 'on');
    h = [plot(ax1, NaN, NaN, 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k'), ...
         plot(ax1, NaN, NaN, 'x', 'Color', 'm', 'LineWidth', 2)];
    lg = legend(h, {'recovered (minimum loss)', 'true'}, 'Orientation', 'horizontal', 'Box', 'off', 'FontSize', 11);
    lg.Layout.Tile = 'south';
    hold(ax1, 'off');
    exportgraphics(fig, fullfile(set_outdir, sprintf('joint_%s.png', rc1.set_ids(s))), 'Resolution', 130);
    close(fig);
end

function name = panelName(r, m, n_out)
% model name, plus the output the loss is on when several outputs are shown
    name = char(r.model_labels(m));
    if n_out > 1, name = sprintf('%s, %s', name, r.output_label); end
end

function [to_axis, label] = paramAxis(p)
% plot axis for a parameter: log10 p for log-searched parameters, p otherwise
    if p.log
        to_axis = @(x) log10(x);
        label = ['log_{10} ' p.label];
    else
        to_axis = @(x) x;
        label = p.label;
    end
end

function rgb = hex2rgb(hex)
    h = char(hex);
    rgb = hex2dec({h(2:3), h(4:5), h(6:7)})' / 255;
end

function plotNetworkValidation(out_dir)
% Figures comparing every surrogate network simulation with its numeric run
% (same setting_id), for the interior, edge and out-of-range groups.
%
% Error measures (t = one simulated time point; t = 0 is skipped, the runs
% are identical there by construction):
%   NMSE(t)          sum_i (s_i - n_i)^2 / sum_i n_i^2 over regions i (Eq. 7 at one time point)
%   node error e_i(t)  (s_i - n_i)^2 / mean_j n_j^2, so mean_i e_i(t) = NMSE(t)
%   edge residual r_e(t)  (s_e - n_e) / max_t sqrt(mean_e' n_e'(t)^2), signed, for
%                    phi_tau (F) and Delta_tau (W1); the scale is one number per
%                    simulation, the largest RMS of the numeric values over the
%                    edge set (seed edges, or all edges of the connectome) and time
% s = surrogate, n = numeric. Seed edges are the edges into and out of the
% seed region (the region with the largest N at t = 0).
% Group curves are means over the group's simulations; a simulation whose
% surrogate N is not finite is left out (counted in the titles).
%
% Writes to <out_dir>/summary/:
%   fig_nmse_time.png            NMSE(t) per group: median over simulations, IQR band
%   fig_scatter_<group>.png      surrogate vs numeric N at months 1-5, every
%                                region of every simulation (as Fig. 8 of the paper)
%   fig_errors_<group>.png       models (columns) x node error per region and its
%                                mean, seed-edge and all-edge residuals of phi_tau and Delta_tau
%
%   plotNetworkValidation()               % sntm_config out_dir (outputs/)
%   plotNetworkValidation('outputs')

if nargin < 1, cfg = sntm_config(); out_dir = cfg.out_dir; end
sum_dir = fullfile(out_dir, 'summary');
if ~isfolder(sum_dir), mkdir(sum_dir); end

methods = ["mlp", "sr", "gp", "interp", "lr"];
mlabels = ["MLP", "SR", "GP", "NN-int.", "LR"];
colors = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"];
groups = ["interior", "edge", "ood"];
glabels = ["Interior", "Edge", "Out-of-range"];
months = 1:5;   % scatter time points, 30 days apart

%% 1. Pair every surrogate run with its numeric run
files = dir(fullfile(out_dir, '*.mat'));
meta = cell(numel(files), 1);
for k = 1:numel(files)
    S = load(fullfile(files(k).folder, files(k).name), 'meta');
    meta{k} = S.meta;
end
meta = [meta{:}];
case_method = string({meta.method});
case_group = string({meta.group});
case_setting = string({meta.setting_id});

% per group and model: the simulations' error curves
R = struct();
for g = 1:numel(groups)
    num_k = find(case_method == "numeric" & case_group == groups(g));
    for m = 1:numel(methods)
        % seed edges differ between seeds, so they are kept per simulation (cells)
        r = struct('nmse', [], 'node', [], 'F', [], 'W1', [], 'seedF', {{}}, 'seedW1', {{}}, ...
            'seed_out', {{}}, 'seed', [], 'scatter', [], 'n_sims', 0, 'n_nonfinite', 0);
        for k = num_k
            j = find(case_method == methods(m) & case_setting == case_setting(k), 1);
            if isempty(j), continue, end
            n = load(fullfile(files(k).folder, files(k).name), 'N', 'F', 'W1', 'C', 'trange');
            s = load(fullfile(files(j).folder, files(j).name), 'N', 'F', 'W1');
            r.n_sims = r.n_sims + 1;
            if ~all(isfinite(s.N(:)))
                r.n_nonfinite = r.n_nonfinite + 1;
                continue
            end
            days = n.trange * 180;
            keep = 2:numel(days);   % skip t = 0
            [~, seed] = max(n.N(:,1));
            edges = find(n.C > 0 & ~eye(size(n.C)));   % C keeps 24 self-connections; they carry no flux between regions
            [ei, ej] = ind2sub(size(n.C), edges);

            sq = (s.N(:,keep) - n.N(:,keep)).^2;
            r.nmse(:, end+1) = sum(sq, 1)' ./ sum(n.N(:,keep).^2, 1)';
            r.node(:,:, end+1) = sq ./ mean(n.N(:,keep).^2, 1);
            r.F(:,:, end+1) = edgeResidual(s.F, n.F, edges, keep);
            r.W1(:,:, end+1) = edgeResidual(s.W1, n.W1, edges, keep);
            at_seed = ei == seed | ej == seed;
            r.seedF{end+1} = edgeResidual(s.F, n.F, edges(at_seed), keep);
            r.seedW1{end+1} = edgeResidual(s.W1, n.W1, edges(at_seed), keep);
            r.seed_out{end+1} = ei(at_seed) == seed;   % true: out of the seed, false: into it
            r.seed(end+1) = seed;
            [~, idx] = min(abs(days(:) - 30*months), [], 1);
            r.scatter(:,:, end+1, :) = permute(cat(3, n.N(:,idx), s.N(:,idx)), [1 2 4 3]);   % regions x months x sim x [numeric, surrogate]
        end
        if ~isempty(r.nmse)
            % first slice of each 3-D stack is the zero the stack was started from
            r.node = r.node(:,:, 2:end);  r.F = r.F(:,:, 2:end);  r.W1 = r.W1(:,:, 2:end);
            r.scatter = r.scatter(:,:, 2:end, :);
            r.days = days(keep);
        end
        R.(groups(g)).(methods(m)) = r;
    end
end

%% 2. NMSE(t), all groups
fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 36 10]);
tl = tiledlayout(fig, 1, numel(groups), 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, 'Network NMSE of N over simulated time (median over simulations, IQR band)', 'FontSize', 14);
for g = 1:numel(groups)
    ax = nexttile(tl); hold(ax, 'on');
    h = gobjects(0);  names = strings(0);
    for m = 1:numel(methods)
        r = R.(groups(g)).(methods(m));
        if isempty(r.nmse), continue, end
        rgb = hex2rgb(colors(m));
        q = prctile(r.nmse, [25 50 75], 2);
        fill(ax, [r.days(:); flipud(r.days(:))], [q(:,1); flipud(q(:,3))], rgb, ...
            'FaceAlpha', 0.15, 'EdgeColor', 'none');
        h(end+1) = plot(ax, r.days, q(:,2), 'Color', rgb, 'LineWidth', 2); %#ok<AGROW>
        names(end+1) = sprintf('%s (n = %d%s)', mlabels(m), size(r.nmse, 2), nonfiniteNote(r)); %#ok<AGROW>
    end
    set(ax, 'YScale', 'log', 'Box', 'off', 'TickDir', 'out', 'XLim', [0 180], 'XTick', 0:30:180);
    grid(ax, 'on');
    xlabel(ax, 'days');
    if g == 1, ylabel(ax, 'NMSE(t)'); end
    title(ax, glabels(g), 'FontWeight', 'normal');
    if isempty(h)
        text(ax, 0.5, 0.5, 'no surrogate runs found', 'Units', 'normalized', 'HorizontalAlignment', 'center');
    end
    if ~isempty(h), legend(ax, h, names, 'Location', 'best', 'Box', 'off'); end
    hold(ax, 'off');
end
exportgraphics(fig, fullfile(sum_dir, 'fig_nmse_time.png'), 'Resolution', 200);
close(fig);

%% 3. Per group: scatter at months 1-5, and the error panels
for g = 1:numel(groups)
    present = arrayfun(@(m) ~isempty(R.(groups(g)).(m).nmse), methods);
    if ~any(present), continue, end
    gm = methods(present);  gl = mlabels(present);  gc = colors(present);

    % surrogate vs numeric N (Fig. 8 of the paper), all regions of all simulations
    fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 8*numel(months) + 4, 7.5*numel(gm) + 2]);
    tl = tiledlayout(fig, numel(gm), numel(months), 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, sprintf('%s: regional N, surrogate vs numeric', glabels(g)), 'FontSize', 14);
    xlabel(tl, 'N, numeric');  ylabel(tl, 'N, surrogate');
    for m = 1:numel(gm)
        r = R.(groups(g)).(gm(m));
        for t = 1:numel(months)
            ax = nexttile(tl); hold(ax, 'on');
            x = reshape(r.scatter(:,t,:,1), [], 1);  y = reshape(r.scatter(:,t,:,2), [], 1);
            lim = [0, max([x; y])];
            plot(ax, lim, lim, '-k', 'LineWidth', 0.75);
            plot(ax, x, y, '.', 'Color', hex2rgb(gc(m)), 'MarkerSize', 8);
            set(ax, 'Box', 'on', 'TickDir', 'out', 'XLim', lim, 'YLim', lim);
            axis(ax, 'square');
            text(ax, 0.05, 0.92, sprintf('NMSE = %.3g', sum((y - x).^2) / sum(x.^2)), ...
                'Units', 'normalized', 'FontSize', 10);
            if m == 1, title(ax, sprintf('Month %d', months(t)), 'FontWeight', 'normal'); end
            if t == 1, ylabel(ax, gl(m), 'FontWeight', 'bold'); end
            hold(ax, 'off');
        end
    end
    exportgraphics(fig, fullfile(sum_dir, sprintf('fig_scatter_%s.png', groups(g))), 'Resolution', 150);
    close(fig);

    % error panels: rows = node error, seed edges (phi_tau, Delta_tau), all edges (phi_tau, Delta_tau)
    rows = {'node', 'node error e_i(t) (mean over simulations)'; ...
            'seedF', 'seed edges: \phi_\tau residual'; ...
            'seedW1', 'seed edges: \Delta_\tau residual'; ...
            'allF', 'all edges: \phi_\tau residual'; ...
            'allW1', 'all edges: \Delta_\tau residual'};
    fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 9*numel(gm) + 4, 7*size(rows, 1) + 2]);
    tl = tiledlayout(fig, size(rows, 1), numel(gm), 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, sprintf('%s: surrogate minus numeric over simulated time', glabels(g)), 'FontSize', 14);
    for p = 1:size(rows, 1)
        for m = 1:numel(gm)
            r = R.(groups(g)).(gm(m));
            rgb = hex2rgb(gc(m));
            ax = nexttile(tl); hold(ax, 'on');
            switch rows{p,1}
                case 'node'
                    E = mean(r.node, 3);   % regions x time
                    plot(ax, r.days, E', 'Color', [rgb 0.35], 'LineWidth', 0.75);
                    plot(ax, r.days, mean(E, 1), '-k', 'LineWidth', 2);
                    set(ax, 'YScale', 'log');
                case {'seedF', 'seedW1'}
                    % mean over the simulations of each seed, then every edge of every seed
                    for seed = unique(r.seed)
                        same_seed = find(r.seed == seed);
                        Q = mean(cat(3, r.(rows{p,1}){same_seed}), 3);   % seed edges x time
                        out = r.seed_out{same_seed(1)};
                        plot(ax, r.days, Q(out,:)', '-', 'Color', [rgb 0.6], 'LineWidth', 0.75);
                        plot(ax, r.days, Q(~out,:)', '--', 'Color', [rgb 0.6], 'LineWidth', 0.75);
                    end
                    yline(ax, 0, '-k');
                case {'allF', 'allW1'}
                    Q = mean(r.(rows{p,1}(4:end)), 3);
                    q = prctile(Q, [10 50 90], 1);
                    fill(ax, [r.days(:); flipud(r.days(:))], [q(1,:)'; flipud(q(3,:)')], rgb, ...
                        'FaceAlpha', 0.2, 'EdgeColor', 'none');
                    plot(ax, r.days, q(2,:), 'Color', rgb, 'LineWidth', 1.5);
                    plot(ax, r.days, mean(Q, 1), '-k', 'LineWidth', 1.5);
                    yline(ax, 0, ':k');
            end
            set(ax, 'Box', 'off', 'TickDir', 'out', 'XLim', [0 180], 'XTick', 0:30:180);
            grid(ax, 'on');
            if p == 1
                title(ax, sprintf('%s (n = %d%s)', gl(m), size(r.nmse, 2), nonfiniteNote(r)), 'FontWeight', 'normal');
            end
            if m == 1, ylabel(ax, rows{p,2}); end
            if p == size(rows, 1), xlabel(ax, 'days'); end
            hold(ax, 'off');
        end
    end
    lg_ax = nexttile(tl, 1); hold(lg_ax, 'on');
    h = [plot(lg_ax, NaN, NaN, '-', 'Color', [0.5 0.5 0.5]), plot(lg_ax, NaN, NaN, '-k', 'LineWidth', 2), ...
         plot(lg_ax, NaN, NaN, '--', 'Color', [0.5 0.5 0.5])];
    lg = legend(h, {'one region / out of the seed (seed edges) / median over edges (all edges, 10-90% band)', ...
        'mean (over regions = NMSE(t), or over all edges)', 'into the seed (seed edges)'}, ...
        'Orientation', 'horizontal', 'Box', 'off', 'FontSize', 10);
    lg.Layout.Tile = 'south';
    hold(lg_ax, 'off');
    exportgraphics(fig, fullfile(sum_dir, sprintf('fig_errors_%s.png', groups(g))), 'Resolution', 150);
    close(fig);
end
fprintf('wrote figures to %s\n', sum_dir);
end

function Q = edgeResidual(S, N, edges, keep)
% (s - n) on the given edges, divided by the largest RMS of n over those edges
% at any time point (one scale for the whole simulation); edges x time
S = reshape(S, [], size(S, 3));  N = reshape(N, [], size(N, 3));
S = S(edges, keep);  N = N(edges, keep);
Q = (S - N) / max(sqrt(mean(N.^2, 1)));
% the simulator leaves W1 at 0 at the last time point: not a residual of 0
Q(:, all(N == 0, 1) & all(S == 0, 1)) = NaN;
end

function txt = nonfiniteNote(r)
txt = '';
if r.n_nonfinite > 0, txt = sprintf(', %d non-finite left out', r.n_nonfinite); end
end

function rgb = hex2rgb(hex)
hex = char(hex);
rgb = hex2dec({hex(2:3), hex(4:5), hex(6:7)})' / 255;
end

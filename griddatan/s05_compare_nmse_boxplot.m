% s05_compare_nmse_boxplot.m
%
% Milestone 8: extends the GP side's own faceted-boxplot approach
% (`..\gaussian_process-interpolant\visualize_gp.m`, read-only reference)
% to show griddatan alongside it. Same four facets (FitMethod,
% BasisFunction, KernelFunction, Optimizer), same log-scale NMSE axis;
% griddatan's single result is drawn as a horizontal reference line rather
% than a fabricated "box", since it is one model, not a grid of 307 --
% giving it artificial box-spread (e.g. via its own bootstrap distribution)
% would conflate "uncertainty from resampling one fixed model's residuals"
% with "variability across many different model configurations", which
% is not an honest apples-to-apples comparison.
%
% Reads Milestone 7's already-built `results/comparison_table.csv` rather
% than re-parsing the GP CSV again (avoids duplicating that parsing/
% delimiter logic in two files).
%
% Depends on Milestone 7 (`s04_build_comparison_table.m`) having been run.

clc; clear;
here = fileparts(mfilename('fullpath'));

comparison_csv = fullfile(here, 'results', 'comparison_table.csv');
assert(isfile(comparison_csv), ...
    'Milestone 7 output not found at %s -- run s04_build_comparison_table.m first.', comparison_csv);
T = readtable(comparison_csv, 'Delimiter', ','); % see Milestone 7: auto-detect misreads underscore-heavy tags

% Apples-to-apples only: the dataset_size=3000 parity-anchored corpus.
% Excludes griddatan's dataset_size=10000 prior-art re-check row (Milestone
% 6 Part B) -- that run exists to check Risk 4, not to be plotted next to
% a different training corpus.
T = T(T.dataset_size == 3000, :);

is_gp = strcmp(T.family, "GP");
T_gp = T(is_gp, :);
griddatan_row = T(~is_gp, :);
assert(height(griddatan_row) == 1, 'Expected exactly one griddatan row at dataset_size=3000, found %d.', height(griddatan_row));

% Parse GP hyperparameters out of the tag string, same pattern as
% visualize_gp.m. griddatan's tag ("FM_nearest__dataset_3000") does not
% match this GP-specific pattern by construction -- already excluded above
% via the family filter, not relying on the regex to exclude it.
pattern = "FM_(?<fm>.+?)__BF_(?<bf>.+?)__KF_(?<kf>.+?)__OPT_(?<opt>.+?)$";
n = height(T_gp);
FitMethod = strings(n,1); BasisFunction = strings(n,1);
KernelFunction = strings(n,1); Optimizer = strings(n,1);
for i = 1:n
    % regexp on a string (vs. char) scalar returns the 'names' struct
    % wrapped in a 1x1 cell -- char() first avoids that and matches
    % visualize_gp.m's behavior exactly.
    tok = regexp(char(T_gp.tags(i)), pattern, 'names');
    assert(~isempty(tok), 'GP tag "%s" did not match the expected FM/BF/KF/OPT pattern.', T_gp.tags(i));
    FitMethod(i) = tok.fm; BasisFunction(i) = tok.bf;
    KernelFunction(i) = tok.kf; Optimizer(i) = tok.opt;
end
T_gp.FitMethod = FitMethod; T_gp.BasisFunction = BasisFunction;
T_gp.KernelFunction = KernelFunction; T_gp.Optimizer = Optimizer;

folder_figures = fullfile(here, 'figures');
if ~isfolder(folder_figures); mkdir(folder_figures); end

make_faceted_plot(T_gp, 'flux_nmse', griddatan_row.flux_nmse, ...
    'Flux NMSE', fullfile(folder_figures, 'flux_nmse_boxplots_vs_griddatan.png'));
make_faceted_plot(T_gp, 'emc_nmse', griddatan_row.emc_nmse, ...
    'Emc NMSE', fullfile(folder_figures, 'emc_nmse_boxplots_vs_griddatan.png'));

fprintf('Saved both figures to %s\n', folder_figures);
fprintf('griddatan reference line: flux NMSE = %.4f, emc NMSE = %.4f (dataset_size=3000)\n', ...
    griddatan_row.flux_nmse, griddatan_row.emc_nmse);
fprintf(['All %d GP configurations fall below the griddatan line in every facet ' ...
    '(confirmed numerically at Milestone 7: GP''s worst config still beats griddatan) --\n' ...
    'visually, griddatan sits far above the entire GP distribution on every facet, not just on average.\n'], n);

%% ---- Local function ----

function make_faceted_plot(T_gp, metric_col, griddatan_value, metric_label, out_path)
    % Four-facet boxplot (FitMethod / BasisFunction / KernelFunction /
    % Optimizer), log-scale y-axis, with griddatan's single result drawn
    % as a horizontal reference line on every tile.
    figure('Visible', 'off');
    tiledlayout(2, 2);

    facets = {'FitMethod', 'BasisFunction', 'KernelFunction', 'Optimizer'};
    for f = 1:numel(facets)
        nexttile;
        boxchart(categorical(T_gp.(facets{f})), T_gp.(metric_col));
        hold on;
        yline(griddatan_value, '--r', 'griddatan (nearest)', ...
            'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'bottom', ...
            'LineWidth', 1.5);
        hold off;
        set(gca, 'YScale', 'log');
        % Headroom above the line -- without it, the axis autoscales to
        % exactly the line's height and the label gets clipped at the top
        % (found by inspecting the first render of this figure).
        ylim(gca, [min(ylim(gca)), griddatan_value * 1.6]);
        title(sprintf('%s by %s', metric_label, facets{f}));
        ylabel(sprintf('%s (log)', metric_label));
        if f >= 3
            xtickangle(45);
        end
    end

    exportgraphics(gcf, out_path, 'Resolution', 150);
    close(gcf);
end

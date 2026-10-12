% Faceted box plot to show GP model combo NMSE
% Brandon, 10/12/2026

clc;
clear;

fprintf("Running visualize_gp.m... \n")

%dataset_size = 5000;   % must match the folder the grid-search run produced
folder_name = sprintf("saved_models");

files = dir(fullfile(folder_name, "FM_*__BF_*__KF_*__OPT_*.mat"));
n = numel(files);

FitMethod      = strings(n,1);
BasisFunction  = strings(n,1);
KernelFunction = strings(n,1);
Optimizer      = strings(n,1);
flux_nmse      = zeros(n,1);
emc_nmse       = zeros(n,1);

pattern = "FM_(?<fm>.+?)__BF_(?<bf>.+?)__KF_(?<kf>.+?)__OPT_(?<opt>.+?)\.mat";

fprintf("Now looking through folder... \n")
for i = 1:n
    tok = regexp(files(i).name, pattern, 'names');
    FitMethod(i)      = tok.fm;
    BasisFunction(i)  = tok.bf;
    KernelFunction(i) = tok.kf;
    Optimizer(i)      = tok.opt;

    d = load(fullfile(files(i).folder, files(i).name), 'flux_nmse_gpr', 'emc_nmse_gpr');
    flux_nmse(i) = d.flux_nmse_gpr;
    emc_nmse(i)  = d.emc_nmse_gpr;

    if mod(i,200) == 0

    end
end

T = table(FitMethod, BasisFunction, KernelFunction, Optimizer, flux_nmse, emc_nmse);

figure();
tiledlayout(2,2);

nexttile;
boxchart(categorical(T.FitMethod), T.flux_nmse);
set(gca, 'YScale', 'log');
title('Flux NMSE by Fit Method'); ylabel('Flux NMSE (log)');

nexttile;
boxchart(categorical(T.BasisFunction), T.flux_nmse);
set(gca, 'YScale', 'log');
title('Flux NMSE by Basis Function'); ylabel('Flux NMSE (log)');

nexttile;
boxchart(categorical(T.KernelFunction), T.flux_nmse);
set(gca, 'YScale', 'log');
title('Flux NMSE by Kernel Function'); ylabel('Flux NMSE (log)');
xtickangle(45);

nexttile;
boxchart(categorical(T.Optimizer), T.flux_nmse);
set(gca, 'YScale', 'log');
title('Flux NMSE by Optimizer'); ylabel('Flux NMSE (log)');

exportgraphics(gcf, fullfile(folder_name, 'flux_nmse_boxplots.png'), 'Resolution', 150);
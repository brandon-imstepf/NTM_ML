% Gaussian Process & Gridded Interpolant
clc;
clear;

% options:
% grid-search - do a grid search over all possible argument combinations
% best-model - run only on the best performing models found
grid_search = 0;
best_model = 1;
assert(grid_search + best_model < 2 && grid_search + best_model > 0, "Error with Option input.")

% load datasets 
tr_filepath = fullfile('..','data/data_neural_network_model/training_data_e6/','data_bias_e5_1.csv'); %100k rows
training_superset = readmatrix(tr_filepath);
% subset data based off seed
training_seed = 42;
rng(training_seed)
dataset_size = 3000;
training_idx = randperm(100000,dataset_size);
training_dataset = training_superset(training_idx,:);

% load validation
va_filepath = fullfile('..','data/data_baseline_models/e3_data/','bias_val_e3.csv');
validation_dataset = readmatrix(va_filepath);

% load training
features = training_dataset(:,1:6);
flux = training_dataset(:,7);
emc = training_dataset(:,8);


% need to standardize the responses -- they are too small!!
flux_mean = mean(flux); flux_std = std(flux);
flux_z = (flux - flux_mean) / flux_std; % z-scored flux

emc_mean = mean(emc); emc_std = std(emc);
emc_z = (emc - emc_mean) / emc_std; % z-scored emc


% all arguments possible for fitrgp
FitMethodList = ["none","exact","sd","sr","fic"];
BasisFunctionList = ["none","constant","linear","pureQuadratic"];
KernelFunctionList = ["Exponential","SquaredExponential","Matern32",...
    "Matern52","RationalQuadratic","ARDExponential","ARDSquaredExponential",...
    "ARDMatern32","ARDMatern52","ARDRationalQuadratic"];
OptimizerList = ["quasinewton","lbfgs","fminsearch","fminunc","fmincon"];

% Best performing model list -
best_models = [
    "exact","linear","ARDMatern32","lbfgs";
    "sd",   "linear","ARDMatern32","lbfgs";
    "exact","linear","ARDMatern32","fmincon";
    "sd",   "linear","ARDMatern32","fmincon";
    "exact","linear","ARDMatern32","fminunc";
    "sd",   "linear","ARDMatern32","fminunc";
    "exact","linear","ARDMatern32","quasinewton";
    "sd",   "linear","ARDMatern32","quasinewton";
    "fic",  "linear","ARDMatern32","lbfgs";
    "fic",  "linear","ARDMatern32","fmincon"
];


% load validation
val_features = validation_dataset(:,1:6);
val_flux = validation_dataset(:,7);
val_emc = validation_dataset(:,8);

folder_name = sprintf("saved_models_superset_%d", dataset_size);
mkdir(folder_name);
Results = table();
row = 1;

% bootstrapping parameters
n_boot = 1000;
boot_size = 1000;

totalcount = numel(FitMethodList)*numel(BasisFunctionList) ...
    *numel(KernelFunctionList)*numel(OptimizerList);

fprintf("Subset size = %d \n",dataset_size);


if grid_search == 1
    fprintf("Performing Grid-Search")

    for i = 1:length(FitMethodList)
    
        for j = 1:length(BasisFunctionList)
    
            for k = 1:length(KernelFunctionList)
    
                for l = 1:length(OptimizerList)
    
                    tag = sprintf("FM_%s__BF_%s__KF_%s__OPT_%s", ...
                        FitMethodList{i}, BasisFunctionList{j}, ...
                        KernelFunctionList{k}, OptimizerList{l});
    
                    tic
    
                    fprintf("gpr model %d/%d using %s, %s, %s, %s... \n", ...
                        row,totalcount,FitMethodList(i),BasisFunctionList(j),KernelFunctionList(k), ...
                        OptimizerList(l));
                    stzd_flux_gprModel = fitrgp(features,flux_z, 'Standardize',true, ...
                        'FitMethod',FitMethodList(i),'BasisFunction',BasisFunctionList(j), ...
                        'KernelFunction',KernelFunctionList(k),'Optimizer',OptimizerList(l));
                    fprintf("gpr model for emc using %s, %s, %s, %s... \n", ...
                        FitMethodList(i),BasisFunctionList(j),KernelFunctionList(k), ...
                        OptimizerList(l));
                    stzd_emc_gprModel = fitrgp(features,emc_z,'Standardize',true, ...
                        'FitMethod',FitMethodList(i),'BasisFunction',BasisFunctionList(j), ...
                        'KernelFunction',KernelFunctionList(k),'Optimizer',OptimizerList(l));
                    % % use matlabs predict()
                    [gpr_flux_ypred_z, gpr_flux_ysd] = predict(stzd_flux_gprModel, val_features);
                    [gpr_emc_ypred_z, gpr_emc_ysd] = predict(stzd_emc_gprModel, val_features);
    
                    gpr_flux_ypred = gpr_flux_ypred_z * flux_std + flux_mean;
                    gpr_emc_ypred  = gpr_emc_ypred_z  * emc_std  + emc_mean;
    
                    % training-set predictors (for bootstrapping)
                    [gpr_flux_xpred_z, gpr_flux_xsd] = predict(stzd_flux_gprModel, features);
                    [gpr_emc_xpred_z, gpr_emc_xsd] = predict(stzd_emc_gprModel, features);
    
                    gpr_flux_xpred = gpr_flux_xpred_z * flux_std + flux_mean;
                    gpr_emc_xpred  = gpr_emc_xpred_z  * emc_std  + emc_mean;
    
                    % Squared error into NMSE
                    squared_error_flux_y = (val_flux - gpr_flux_ypred).^2;
                    squared_error_emc_y = (val_emc - gpr_emc_ypred).^2;
                    squared_error_flux_x = (flux - gpr_flux_xpred).^2;
                    squared_error_emc_x = (emc - gpr_emc_xpred).^2;
    
                    flux_nmse_gpr = mean(squared_error_flux_y) / var(val_flux);
                    emc_nmse_gpr = mean(squared_error_emc_y) / var(val_emc);
    
                    wall_time = toc;
    
                    boot_flux_nmse_test = zeros(n_boot,1);
                    for b = 1:n_boot
                        boot_idx = randi(numel(squared_error_flux_y),boot_size,1);
                        boot_flux_nmse_test(b) = mean(squared_error_flux_y(boot_idx)) / var(val_flux);
                    end
                    boot_flux_nmse_train = zeros(n_boot,1);
                    for b = 1:n_boot
                        boot_idx = randi(numel(squared_error_flux_x),boot_size,1);
                        boot_flux_nmse_train(b) = mean(squared_error_flux_x(boot_idx)) / var(flux);
                    end
                    boot_emc_nmse_test = zeros(n_boot,1);
                    for b = 1:n_boot
                        boot_idx = randi(numel(squared_error_emc_y),boot_size,1);
                        boot_emc_nmse_test(b) = mean(squared_error_emc_y(boot_idx)) / var(val_emc);
                    end
                    boot_emc_nmse_train = zeros(n_boot,1);
                    for b = 1:n_boot
                        boot_idx = randi(numel(squared_error_emc_x),boot_size,1);
                        boot_emc_nmse_train(b) = mean(squared_error_emc_x(boot_idx)) / var(emc);
                    end
    
    
                    % record everything in the table for later use
                    Results(row,:) = {
                        FitMethodList(i), ...
                        BasisFunctionList(j),...
                        KernelFunctionList(k), ...
                        OptimizerList(l),...
                        flux_nmse_gpr, ...
                        emc_nmse_gpr};
    
                    row = row + 1;
    
                    % save everything
                    filename = fullfile(folder_name, tag + ".mat");
                    save(filename, ...
                        "stzd_flux_gprModel", "stzd_emc_gprModel", ...
                        "gpr_flux_ypred", "gpr_emc_ypred", ...
                        "gpr_flux_ysd", "gpr_emc_ysd", ...
                        "gpr_flux_xsd", "gpr_emc_xsd", ...
                        "flux_nmse_gpr", "emc_nmse_gpr","wall_time","dataset_size", ...
                        "boot_flux_nmse_train","boot_flux_nmse_test",...
                        "boot_emc_nmse_train","boot_emc_nmse_test", ...
                        "flux_mean", "flux_std","emc_mean","emc_std");
                end
            end
        end
    end
    
    files = dir(fullfile(folder_name,'*.mat'));
    n = numel(files);
    tags = strings(n,1);
    flux_nmse = zeros(n,1);
    emc_nmse  = zeros(n,1);
    for i = 1:n
        d = load(fullfile(files(i).folder, files(i).name), 'flux_nmse_gpr', 'emc_nmse_gpr');
        tags(i) = erase(files(i).name, '.mat');
        flux_nmse(i) = d.flux_nmse_gpr;
        emc_nmse(i)  = d.emc_nmse_gpr;
    end
    T = table(tags, flux_nmse, emc_nmse);
    T = sortrows(T, 'flux_nmse');
    writetable(T, sprintf('gp_comparison_results_superset_%d.csv',dataset_size))
    head(T, 10)   % best 10 configs by flux NMSE

end 

if best_model > 0

    fprintf("Performing best_model...")
    totalcount_bestmodels = size(best_models,1);
    
    for m = 1:totalcount_bestmodels
    
        tag = sprintf("FM_%s__BF_%s__KF_%s__OPT_%s", ...
            best_models(m,1), best_models(m,2), best_models(m,3), best_models(m,4));
    
        tic
    
        fprintf("gpr model %d/%d using %s, %s, %s, %s... \n", ...
            row,totalcount_bestmodels,best_models(m,1), best_models(m,2), best_models(m,3), best_models(m,4));
        stzd_flux_gprModel = fitrgp(features,flux_z, 'Standardize',true, ...
            'FitMethod',best_models(m,1),'BasisFunction',best_models(m,2), ...
            'KernelFunction',best_models(m,3),'Optimizer',best_models(m,4));
        fprintf("gpr model for emc using %s, %s, %s, %s... \n", ...
            best_models(m,1), best_models(m,2), best_models(m,3), best_models(m,4));
        stzd_emc_gprModel = fitrgp(features,emc_z,'Standardize',true, ...
            'FitMethod',best_models(m,1),'BasisFunction',best_models(m,2), ...
            'KernelFunction',best_models(m,3),'Optimizer',best_models(m,4));
        % % use matlabs predict()
        [gpr_flux_ypred_z, gpr_flux_ysd] = predict(stzd_flux_gprModel, val_features);
        [gpr_emc_ypred_z, gpr_emc_ysd] = predict(stzd_emc_gprModel, val_features);
    
        gpr_flux_ypred = gpr_flux_ypred_z * flux_std + flux_mean;
        gpr_emc_ypred  = gpr_emc_ypred_z  * emc_std  + emc_mean;
    
        % training-set predictors (for bootstrapping)
        [gpr_flux_xpred_z, gpr_flux_xsd] = predict(stzd_flux_gprModel, features);
        [gpr_emc_xpred_z, gpr_emc_xsd] = predict(stzd_emc_gprModel, features);
    
        gpr_flux_xpred = gpr_flux_xpred_z * flux_std + flux_mean;
        gpr_emc_xpred  = gpr_emc_xpred_z  * emc_std  + emc_mean;
    
        % Squared error into NMSE
        squared_error_flux_y = (val_flux - gpr_flux_ypred).^2;
        squared_error_emc_y = (val_emc - gpr_emc_ypred).^2;
        squared_error_flux_x = (flux - gpr_flux_xpred).^2;
        squared_error_emc_x = (emc - gpr_emc_xpred).^2;
    
        flux_nmse_gpr = mean(squared_error_flux_y) / var(val_flux);
        emc_nmse_gpr = mean(squared_error_emc_y) / var(val_emc);
    
        wall_time = toc;
    
        boot_flux_nmse_test = zeros(n_boot,1);
        for b = 1:n_boot
            boot_idx = randi(numel(squared_error_flux_y),boot_size,1);
            boot_flux_nmse_test(b) = mean(squared_error_flux_y(boot_idx)) / var(val_flux);
        end
        boot_flux_nmse_train = zeros(n_boot,1);
        for b = 1:n_boot
            boot_idx = randi(numel(squared_error_flux_x),boot_size,1);
            boot_flux_nmse_train(b) = mean(squared_error_flux_x(boot_idx)) / var(flux);
        end
        boot_emc_nmse_test = zeros(n_boot,1);
        for b = 1:n_boot
            boot_idx = randi(numel(squared_error_emc_y),boot_size,1);
            boot_emc_nmse_test(b) = mean(squared_error_emc_y(boot_idx)) / var(val_emc);
        end
        boot_emc_nmse_train = zeros(n_boot,1);
        for b = 1:n_boot
            boot_idx = randi(numel(squared_error_emc_x),boot_size,1);
            boot_emc_nmse_train(b) = mean(squared_error_emc_x(boot_idx)) / var(emc);
        end
    
    
        % record everything in the table for later use
        Results(row,:) = {
            best_models(m,1), ...
            best_models(m,2), ...
            best_models(m,3), ...
            best_models(m,4),...
            flux_nmse_gpr, ...
            emc_nmse_gpr};
    
        row = row + 1;
    
        % save everything
        filename = fullfile(folder_name, tag + ".mat");
        save(filename, ...
            "stzd_flux_gprModel", "stzd_emc_gprModel", ...
            "gpr_flux_ypred", "gpr_emc_ypred", ...
            "gpr_flux_ysd", "gpr_emc_ysd", ...
            "gpr_flux_xsd", "gpr_emc_xsd", ...
            "flux_nmse_gpr", "emc_nmse_gpr","wall_time","dataset_size", ...
            "boot_flux_nmse_train","boot_flux_nmse_test",...
            "boot_emc_nmse_train","boot_emc_nmse_test", ...
            "flux_mean", "flux_std","emc_mean","emc_std");

    end
    
    files = dir(fullfile(folder_name,'*.mat'));
    n = numel(files);
    tags = strings(n,1);
    flux_nmse = zeros(n,1);
    emc_nmse  = zeros(n,1);
    for i = 1:n
        d = load(fullfile(files(i).folder, files(i).name), 'flux_nmse_gpr', 'emc_nmse_gpr');
        tags(i) = erase(files(i).name, '.mat');
        flux_nmse(i) = d.flux_nmse_gpr;
        emc_nmse(i)  = d.emc_nmse_gpr;
    end
    T = table(tags, flux_nmse, emc_nmse);
    T = sortrows(T, 'flux_nmse');
    writetable(T, sprintf('gp_comparison_results_superset_%d.csv',dataset_size))
    head(T, 10)   % best 10 configs by flux NMSE
end 

% 
% % run gpr model
% fprintf("gpr model for Flux... \n");
% stzd_flux_gprModel = fitrgp(features,flux_z, 'Standardize',true);
% fprintf("gpr model for EMC... \n");
% stzd_emc_gprModel = fitrgp(features,emc_z,'Standardize',true);
% 
% % load validation
% val_features = validation_dataset(:,1:6);
% val_flux = validation_dataset(:,7);
% val_emc = validation_dataset(:,8);
% 
% % use matlabs predict()
% [gpr_flux_ypred_z, gpr_flux_ysd] = predict(stzd_flux_gprModel, val_features);
% [gpr_emc_ypred_z, gpr_emc_ysd] = predict(stzd_emc_gprModel, val_features);
% 
% gpr_flux_ypred = gpr_flux_ypred_z * flux_std + flux_mean;
% gpr_emc_ypred  = gpr_emc_ypred_z  * emc_std  + emc_mean;
% 
% % get NMSE
% flux_nmse_gpr = mean((val_flux - gpr_flux_ypred).^2) / var(val_flux);
% emc_nmse_gpr = mean((val_emc - gpr_emc_ypred).^2) / var(val_emc);

% plot to see residuals



% Grid search
% kernel, function


% k nearest neighbors? Might as well!

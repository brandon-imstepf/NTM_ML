% Synthetic Parameter Recovery - Gaussian Process%
% 10/06/26 %
clc;
clear;

% Load GP Model(s) 
gp_model_file_1 = fullfile('..','gaussian_process-interpolant','saved_models_superset_3000', ...
    'FM_exact__BF_linear__KF_ARDMatern32__OPT_fminunc.mat');
gp_1 = load(gp_model_file_1, 'stzd_flux_gprModel', 'stzd_emc_gprModel', 'flux_mean','flux_std','emc_mean','emc_std');

gp_model_file_2 = fullfile('..','gaussian_process-interpolant','saved_models_superset_3000', ...
    'FM_exact__BF_linear__KF_ARDMatern32__OPT_lbfgs.mat');
gp_2 = load(gp_model_file_2, 'stzd_flux_gprModel', 'stzd_emc_gprModel', 'flux_mean','flux_std','emc_mean','emc_std');

gp_model_file_3 = fullfile('..','gaussian_process-interpolant','saved_models_superset_3000', ...
    'FM_exact__BF_linear__KF_ARDMatern32__OPT_fmincon.mat');
gp_3 = load(gp_model_file_3, 'stzd_flux_gprModel', 'stzd_emc_gprModel', 'flux_mean','flux_std','emc_mean','emc_std');

gp_model_file_4 = fullfile('..','gaussian_process-interpolant','saved_models_superset_3000', ...
    'FM_exact__BF_linear__KF_ARDMatern32__OPT_quasinewton.mat');
gp_4 = load(gp_model_file_4, 'stzd_flux_gprModel', 'stzd_emc_gprModel', 'flux_mean','flux_std','emc_mean','emc_std');
% ^ Loaded these 4 as they tied to the 1e-6 on Flux NMSE (0.0073629)
% and EMC NMSE (0.037501). Testing to see which are best.

gp_models = {gp_1, gp_2, gp_3, gp_4};
gp_tags   = ["fminunc","lbfgs","fmincon","quasinewton"];  % all models the same except for OPT



% load holdout validation dataset

filepath = fullfile('..','/data/data_baseline_models/e3_data/','bias_val_e3.csv');
dataset = readmatrix(filepath);

% load matching z-score mean used to fit

flux_scaling = load(fullfile("flux_scaling_params.mat"));
emc_scaling = load(fullfile("emc_scaling_params.mat"));

zf = @(raw,idx) (raw - flux_scaling.feature_means(idx))./ flux_scaling.feature_stds(idx);
ze = @(raw,idx) (raw - emc_scaling.feature_means(idx)) ./ emc_scaling.feature_stds(idx);
% [raw,idx] is [gamma1,1], [lambda1,2],[delta,3],[epsilon,4],[N1,5],[N2,6]

    numiterations = 1000;
    gamma1_range = [5e-5,5e-2]; 
    lambda_range = [0.01, 0.1];
    delta_range = [10, 100];
    epsilon_range = [10, 100];
    logdist_fun = @(x,n) 10.^(log10(x(1)) + (log10(x(2)/x(1)) * rand(n,1)));
    lindist_fun = @(x,n) x(1) + (x(2) - x(1)) * rand(n,1);
    lambdaspace = logdist_fun(lambda_range,numiterations);
    deltaspace = lindist_fun(delta_range,numiterations);
    epsilonspace = lindist_fun(epsilon_range,numiterations);
    gammaspace = logdist_fun(gamma1_range,numiterations);
    gridN = 50; % Important, changed from 200 because super duper long calculations otherwise.
    gamma_grid = sort(logdist_fun(gamma1_range, gridN));
    lambda_grid = sort(logdist_fun(lambda_range, gridN));
    [GG,LL] = ndgrid(gamma_grid, lambda_grid);

function runGPRecovery(gp_model, tag, dataset, flux_scaling, emc_scaling, ...
    gammaspace, lambdaspace, gamma1_range, lambda_range, gamma_grid, lambda_grid, GG, LL)

    outdir = "GP_flux_" + tag;
    numiterations = 1000;

    gp_flux_scaled_1 = @(gamma1,lambda1,delta,epsilon,N1,N2) ...
        ( ( reshape(predict(gp_model.stzd_flux_gprModel, [gamma1(:), lambda1(:), ...
            repmat(delta,numel(gamma1),1), repmat(epsilon,numel(gamma1),1), ...
            repmat(N1,numel(gamma1),1), repmat(N2,numel(gamma1),1)]), size(gamma1)) ...
          * gp_model.flux_std + gp_model.flux_mean ) - flux_scaling.target_mean ) / flux_scaling.target_std;
    
    gp_emc_scaled_1 = @(gamma1,lambda1,delta,epsilon,N1,N2) ...
        ( ( reshape(predict(gp_model.stzd_emc_gprModel, [gamma1(:), lambda1(:), ...
            repmat(delta,numel(gamma1),1), repmat(epsilon,numel(gamma1),1), ...
            repmat(N1,numel(gamma1),1), repmat(N2,numel(gamma1),1)]), size(gamma1)) ...
          * gp_model.emc_std + gp_model.emc_mean ) - emc_scaling.target_mean ) / emc_scaling.target_std;
    
    
    
    

    
    
    % Phase 3 -- loop over the whole thing.
    
    [refined_gamma1, refined_lambda, refined_delta, refined_loss_gamma, ...
        refined_loss_lambda, refined_loss_delta] = deal(zeros(size(dataset,1),1));
    
    
    for row = 1:size(dataset,1)
        g_row = dataset(row,1);
        l_row = dataset(row,2);
        d_row = dataset(row,3);
        e_row = dataset(row,4);
        n1_row = dataset(row,5);
        n2_row = dataset(row,6);
    
        % scaled
        flux_scaled = ( dataset(row,7) - flux_scaling.target_mean ) / flux_scaling.target_std;
        emc_scaled = ( dataset(row,8) - emc_scaling.target_mean ) / emc_scaling.target_std;
        % gamma
        gamma_loss_fun = @(gamma1_candidate) ...
            (gp_flux_scaled_1(gamma1_candidate,l_row,d_row,e_row,n1_row,n2_row) - flux_scaled)^2;
        [refined_gamma1(row), refined_loss_gamma(row)] = fminbnd(gamma_loss_fun,5e-5,5e-2);
        % lambda
        lambda_loss_fun = @(lambda_candidate) ...
            (gp_flux_scaled_1(g_row,lambda_candidate,d_row,e_row,n1_row,n2_row) - flux_scaled)^2;
        [refined_lambda(row), refined_loss_lambda(row)] = fminbnd(lambda_loss_fun,0.01, 0.1);
        % delta
        delta_loss_fun = @(delta_candidate) ...
            (gp_flux_scaled_1(g_row,l_row,delta_candidate,e_row,n1_row,n2_row) - flux_scaled)^2;
        [refined_delta(row), refined_loss_delta(row)] = fminbnd(delta_loss_fun,10, 100);
    
    end
    
    
    % sort and plot gamma data
    [sorted_gamma, gamma_idx] = sort(gammaspace);
    [gamma_minloss,gamma_count] = min(refined_loss_gamma);
    recovered_gamma1 = gammaspace(gamma_count);
    
    
    % gamma_scatter_figure = figure();
    % scatter(dataset(:,1), refined_gamma1)
    % set(gca, 'XScale', 'log', 'YScale', 'log')
    % xlabel('true gamma1')
    % ylabel('recovered gamma1')
    % saveas(gamma_scatter_figure,'scatter_gamma_figure.png')
    % 
    % find rows with the best recovery aka lowest error
    log_error_gamma = abs(log10(refined_gamma1) - log10(dataset(:,1)));
    lowest_count = 10;
    [lowest_count_errors, best_gamma_recovery_rows] = mink(log_error_gamma,lowest_count);
    
    % find row with lowest minimum (just to see, maybe it's not the best
    % recovered??)
    [~, min_min_gamma_row] = min(refined_loss_gamma);
    
    function gammaRecoveryBuilder(dataset,row,gammaspace,gp_flux_scaled,flux_scaling,outdir)
    
        l = dataset(row,2);
        d = dataset(row,3);
        e = dataset(row,4);
        n1 = dataset(row,5);
        n2 = dataset(row,6);
    
        numiterations = 1000;
        % prefilled guesses
        squared_error = zeros(numiterations,1);
        
        for i = 1:numiterations
            % flux part
            gamma_guess_flux = gp_flux_scaled(gammaspace(i),l,d,e,n1,n2);
            flux_scaled = ( dataset(row,7) - flux_scaling.target_mean ) / flux_scaling.target_std;
            gamma_flux_residual = gamma_guess_flux - flux_scaled;
     
        
            % total error
            squared_error(i) = gamma_flux_residual^2; 
        
        end
    
        gamma_loss_fun = @(gamma1_candidate) ...
            (gp_flux_scaled(gamma1_candidate,l,d,e,n1,n2) - flux_scaled)^2;
        [refined_gamma1, refined_loss_gamma] = fminbnd(gamma_loss_fun,5e-5,5e-2);
        
        [sorted_gamma, idx] = sort(gammaspace);
        
        %[minloss,count] = min(squared_error);
        %recovered_gamma1 = gammaspace(count);
        
        figure()
        plot(sorted_gamma, squared_error(idx),'b','LineWidth',3);
        xlabel('Gamma Value','FontSize',15)
        ylabel('Loss','FontSize',15)
        hold on
        loss_thresh = refined_loss_gamma + std(squared_error);
        in_band = sorted_gamma(squared_error(idx) <= loss_thresh);
        band_lo = min(in_band);
        band_hi = max(in_band);
        yl = ylim;
        patch([band_lo band_hi band_hi band_lo], [yl(1) yl(1) yl(2) yl(2)], ...
            [0.5 0.5 0.5], 'FaceAlpha', 0.15, 'EdgeColor','none', ...
            'DisplayName', 'within 1\sigma of min loss');
        plot(refined_gamma1, refined_loss_gamma, 'ro', 'MarkerFaceColor','r')
        xline(dataset(row,1), '--r', sprintf('True Gamma, %.4g', dataset(row,1)));
        %text(refined_gamma1,refined_loss_gamma+0.05,sprintf('   min at %4g', refined_gamma1))
        %text(0.95, 0.95, sprintf('Minimum loss = %4g', refined_loss_gamma), ...
        %'Units', 'normalized', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
        %'BackgroundColor', 'white', 'EdgeColor', 'black', 'Margin', 2, 'FontSize',20)
        h_dot = plot(refined_gamma1, refined_loss_gamma, 'ro', 'MarkerFaceColor', 'r');
        legend(h_dot, sprintf('Minimum loss = %.4g', refined_loss_gamma), 'Location', 'northeast', 'FontSize', 24)
        hold off
        if ~exist(outdir, 'dir'), mkdir(outdir); end
        filename = fullfile(outdir, sprintf('flux_recovery_row%04d.png', row));
        exportgraphics(gcf, filename, 'Resolution', 150);
    
    end
    
    %gammaRecoveryBuilder(dataset,best_gamma_recovery_row,gammaspace,gp_flux_scaled,flux_scaling)
    %gammaRecoveryBuilder(dataset,min_min_gamma_row,gammaspace,gp_flux_scaled_1,flux_scaling)
    
    for row = best_gamma_recovery_rows'
        gammaRecoveryBuilder(dataset,row,gammaspace,gp_flux_scaled_1,flux_scaling,outdir)
    end
    
    % lambda function
    
    function lambdaRecoveryBuilder(dataset,row,lambdaspace,gp_flux_scaled,flux_scaling)
    
        g = dataset(row,1);
        d = dataset(row,3);
        e = dataset(row,4);
        n1 = dataset(row,5);
        n2 = dataset(row,6);
        
        numiterations = 1000;
        % prefilled guesses
        squared_error = zeros(numiterations,1);
        
        for i = 1:numiterations
            % flux part
            lambda_guess_flux = gp_flux_scaled(g,lambdaspace(i),d,e,n1,n2);
            flux_scaled = ( dataset(row,7) - flux_scaling.target_mean ) / flux_scaling.target_std;
            lambda_flux_residual = lambda_guess_flux - flux_scaled;
        
            % total error
            squared_error(i) = lambda_flux_residual^2; 
        
        end
        
        lambda_loss_fun = @(lambda_candidate) ...
            (gp_flux_scaled(g,lambda_candidate,d,e,n1,n2) - flux_scaled)^2;
        [refined_lambda, refined_loss_lambda] = fminbnd(lambda_loss_fun,0.01, 0.1);
        
        [sorted_lambda, idx] = sort(lambdaspace);
        
        %[minloss,count] = min(squared_error);
        %recovered_gamma1 = gammaspace(count);
        
        figure()
        plot(sorted_lambda, squared_error(idx),'b','LineWidth',3);
        xlabel('Lambda Value','FontSize',15)
        ylabel('Loss','FontSize',15)
        hold on
        plot(refined_lambda, refined_loss_lambda, 'ro', 'MarkerFaceColor','r')
        xline(dataset(row,2), '--r', sprintf('True Lambda, %.4g', dataset(row,2)));
        %text(refined_gamma1,refined_loss_gamma+0.05,sprintf('   min at %4g', refined_gamma1))
        %text(0.95, 0.95, sprintf('Minimum loss = %4g', refined_loss_gamma), ...
        %'Units', 'normalized', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
        %'BackgroundColor', 'white', 'EdgeColor', 'black', 'Margin', 2, 'FontSize',20)
        h_dot = plot(refined_lambda, refined_loss_lambda, 'ro', 'MarkerFaceColor', 'r');
        legend(h_dot, sprintf('Minimum loss = %.4g', refined_loss_lambda), 'Location', 'northeast', 'FontSize', 24)
        hold off
        %saveas(gamma_recovery_fig, 'gamma_recovery.png')
    
    end
    
    % % find rows with the best recovery aka lowest error
    % log_error_lambda = abs(log10(refined_lambda) - log10(dataset(:,2)));
    % lowest_count_lambda = 10;
    % [lowest_count_errors_lambda, best_lambda_recovery_rows] = mink(log_error_lambda,lowest_count_lambda);
    % 
    % for row = best_lambda_recovery_rows'
    %     lambdaRecoveryBuilder(dataset,row,lambdaspace,gp_flux_scaled,flux_scaling)
    % end
    
    % for i = 40:60
    %     lambdaRecoveryBuilder(dataset,i,lambdaspace,gp_flux_scaled,tau_emc_scaled,flux_scaling,emc_scaling)
    % end
    

    
    function [gamma_hat, lambda_hat, loss_hat] = jointRecoveryFit( ...
            gamma_grid,lambda_grid,GG,LL,gamma1_range,lambda_range, ...
            gp_flux_scaled,flux_scaled,d_row,e_row,n1_row,n2_row)
    
        loss_surface = (gp_flux_scaled(GG,LL,d_row,e_row,n1_row,n2_row) - flux_scaled).^2;
        [~,lin_idx] = min(loss_surface(:));
        [gi,li] = ind2sub(size(loss_surface),lin_idx);
    
        polish_fun = @(logp) (gp_flux_scaled(10^logp(1),10^logp(2),d_row,e_row,n1_row,n2_row) - flux_scaled)^2;
        [logp_hat, loss_hat] = fminsearch(polish_fun, [log10(gamma_grid(gi)), log10(lambda_grid(li))]);
        logp_hat(1) = min(max(logp_hat(1), log10(gamma1_range(1))), log10(gamma1_range(2)));
        logp_hat(2) = min(max(logp_hat(2), log10(lambda_range(1))), log10(lambda_range(2)));
        gamma_hat  = 10^logp_hat(1);
        lambda_hat = 10^logp_hat(2);
    end
    
    fprintf("Running jointRecoveryFit...\n")
    
    [refined_gamma1_joint, refined_lambda_joint, refined_loss_joint] = deal(zeros(size(dataset,1),1));
    for row = 1:size(dataset,1)/10
        d_row = dataset(row,3); e_row = dataset(row,4);
        n1_row = dataset(row,5); n2_row = dataset(row,6);
        flux_scaled = ( dataset(row,7) - flux_scaling.target_mean ) / flux_scaling.target_std;
        [refined_gamma1_joint(row), refined_lambda_joint(row), refined_loss_joint(row)] = ...
            jointRecoveryFit(gamma_grid,lambda_grid,GG,LL,gamma1_range,lambda_range,gp_flux_scaled_1,flux_scaled,d_row,e_row,n1_row,n2_row);
    end
    
    fprintf("Finding joint error...\n")
    
    log_error_joint = sqrt((log10(refined_gamma1_joint) - log10(dataset(:,1))).^2 + ...
        (log10(refined_lambda_joint) - log10(dataset(:,2))).^2);
    [lowest_joint_errors, best_joint_rows] = mink(log_error_joint, 10);
    
    function jointRecoveryPlot(gamma_grid,lambda_grid,GG,LL,gp_flux_scaled, ...
        flux_scaled,d_row,e_row,n1_row,n2_row,row,g_row,l_row,gamma_hat,lambda_hat,outdir)

        % In case the truth is outside explored range, don't plot at all.
        % if g_row < gamma_grid(1) || g_row > gamma_grid(end) || ...
        %         l_row < lambda_grid(1) || l_row > lambda_grid(end)
        %     fprintf("Skipping row %d: true (gamma,lambda) outside explored grid.\n", row);
        %     return
        % end
        % 
        loss_surface = (gp_flux_scaled(GG,LL,d_row,e_row,n1_row,n2_row) - flux_scaled).^2;
        figure()
        imagesc(log10(lambda_grid), log10(gamma_grid), log10(loss_surface));
        set(gca, 'YDir', 'normal');
        xlabel('log_{10} \lambda_1','FontSize',15); ylabel('log_{10} \gamma_1','FontSize',15);
        colorbar; title(sprintf('Joint loss surface, row %d', row));
        ylabel(colorbar, 'log_{10} squared error','FontSize',15)
        log_loss = log10(loss_surface(:));
        mean_log_loss = mean(log_loss);
        std_log_loss = std(log_loss);
        % move plot to the range
        xlim([min([log10(lambda_grid(1)), log10(l_row), log10(lambda_hat)]), ...
            max([log10(lambda_grid(end)), log10(l_row), log10(lambda_hat)])]);
        ylim([min([log10(gamma_grid(1)), log10(g_row), log10(gamma_hat)]), ...
            max([log10(gamma_grid(end)), log10(g_row), log10(gamma_hat)])]);
        hold on
        contour(log10(lambda_grid), log10(gamma_grid), log10(loss_surface), ...
            [mean_log_loss, mean_log_loss], 'Color','k', 'LineWidth',2, 'LineStyle','--', ...
            'DisplayName','mean log-loss');
        contour(log10(lambda_grid), log10(gamma_grid), log10(loss_surface), ...
            [mean_log_loss-std_log_loss, mean_log_loss-std_log_loss], 'Color','k', 'LineWidth',1, 'LineStyle',':', ...
            'DisplayName','mean \pm 1\sigma log-loss');
        contour(log10(lambda_grid), log10(gamma_grid), log10(loss_surface), ...
            [mean_log_loss+std_log_loss, mean_log_loss+std_log_loss], 'Color','w', 'LineWidth',1, 'LineStyle',':', ...
            'HandleVisibility','off');
        
        plot(log10(lambda_hat), log10(gamma_hat), 'o', 'MarkerFaceColor','r', ...
            'MarkerEdgeColor','k', 'MarkerSize', 18, 'LineWidth', 2, 'DisplayName', 'recovered');
        plot(log10(l_row), log10(g_row), 'x', 'Color','m', ...
            'MarkerSize', 18, 'LineWidth', 3, 'DisplayName', 'true');
        legend('Location', 'northeast','FontSize',15);
        if ~exist(outdir, 'dir'), mkdir(outdir); end
        filename = fullfile(outdir, sprintf('joint_recovery_row%04d.png', row));
        exportgraphics(gcf, filename, 'Resolution', 150);
    end
    
    fprintf("Plotting most accurate...\n")
    
    for row = best_joint_rows'
        g_row = dataset(row,1); l_row = dataset(row,2);
        d_row = dataset(row,3); e_row = dataset(row,4);
        n1_row = dataset(row,5); n2_row = dataset(row,6);
        flux_scaled = ( dataset(row,7) - flux_scaling.target_mean ) / flux_scaling.target_std;
        jointRecoveryPlot(gamma_grid,lambda_grid,GG,LL,gp_flux_scaled_1,flux_scaled,d_row,e_row,n1_row,n2_row,row,g_row,l_row,refined_gamma1_joint(row),refined_lambda_joint(row),outdir);
    end

end


for m = 1:numel(gp_models)
    fprintf('=== Running recovery for %s ===\n', gp_tags(m));
    runGPRecovery(gp_models{m}, gp_tags(m), dataset, flux_scaling, emc_scaling, ...
        gammaspace, lambdaspace, gamma1_range, lambda_range, gamma_grid, lambda_grid, GG, LL);
end
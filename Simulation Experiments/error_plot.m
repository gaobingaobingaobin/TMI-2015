function [ error_opt, error_const, error_RF_compensated, ...
    error_T1_effective, error_SNR ] ...
    = error_plot(index, true_parameter_value, fig_string, ...
        noise_values, berkeley_colors, ...
        parameters_of_interest_est_opt, ...
        parameters_of_interest_est_const, ...
        parameters_of_interest_est_RF_compensated, ...
        parameters_of_interest_est_T1_effective, ...
        parameters_of_interest_est_SNR )

    for param_count = 1:length(noise_values)
        error_opt(param_count)               = sqrt(mean(abs(parameters_of_interest_est_opt(:, index, param_count)            - true_parameter_value).^2)); 
        error_const(param_count)             = sqrt(mean(abs(parameters_of_interest_est_const(:, index, param_count)          - true_parameter_value).^2)); 
        error_RF_compensated(param_count)    = sqrt(mean(abs(parameters_of_interest_est_RF_compensated(:, index, param_count) - true_parameter_value).^2)); 
        error_T1_effective(param_count)      = sqrt(mean(abs(parameters_of_interest_est_T1_effective(:, index, param_count)   - true_parameter_value).^2)); 
        error_SNR(param_count)               = sqrt(mean(abs(parameters_of_interest_est_SNR(:, index, param_count)            - true_parameter_value).^2)); 
    end
    
    figure
    set(gca,'ColorOrder', berkeley_colors([2 4 1 6 3], :), 'NextPlot', 'replacechildren')
    loglog(noise_values, error_T1_effective, 'o-', 'LineWidth', 2)
    hold on
    loglog(noise_values, error_RF_compensated, 'o-', 'LineWidth', 2)
    loglog(noise_values, error_const, 'o-', 'LineWidth', 2)
    loglog(noise_values, error_SNR, 'o-', 'LineWidth', 2)
    loglog(noise_values, error_opt, 'o-', 'LineWidth', 2)
    hold off
    leg = legend('T1 effective', 'RF compensated', 'constant', 'total SNR', 'Fisher information'); 
    xlabel('\sigma^2')
    ylabel('RMS error (average of 25 trials)')
    axis([1e03 1e06 1e-04 1e-01]) 
    set(leg,'FontSize',20);
    set(gca,'FontSize',20);
    tightfig(gcf)    
    print(gcf, '-dpdf', [fig_string '_error.pdf']);


end


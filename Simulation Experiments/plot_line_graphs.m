function [  ] = plot_line_graphs(parameter_values, error_opt_array, error_const_array, error_RF_compensated_array, error_T1_effective_array, error_SNR_array, berkeley_colors, xaxis_labels, axis_limits, axis_type)

for k=1:size(parameter_values, 1) 
    
    figure
    set(gca,'ColorOrder', berkeley_colors([2 4 1 6 3], :), 'NextPlot', 'replacechildren')
    if strcmp(axis_type(k), 'linear')
        plot(parameter_values(k, :), error_T1_effective_array(k, :), 'o-', 'LineWidth', 2)
        hold on
        plot(parameter_values(k, :), error_RF_compensated_array(k, :), 'o-', 'LineWidth', 2)
        plot(parameter_values(k, :), error_const_array(k, :), 'o-', 'LineWidth', 2)
        plot(parameter_values(k, :), error_SNR_array(k, :), 'o-', 'LineWidth', 2)
        plot(parameter_values(k, :), error_opt_array(k, :), 'o-', 'LineWidth', 2)
        hold off
    elseif strcmp(axis_type(k), 'ylog')
        semilogy(parameter_values(k, :), error_T1_effective_array(k, :), 'o-', 'LineWidth', 2)
        hold on
        semilogy(parameter_values(k, :), error_RF_compensated_array(k, :), 'o-', 'LineWidth', 2)
        semilogy(parameter_values(k, :), error_const_array(k, :), 'o-', 'LineWidth', 2)
        semilogy(parameter_values(k, :), error_SNR_array(k, :), 'o-', 'LineWidth', 2)
        semilogy(parameter_values(k, :), error_opt_array(k, :), 'o-', 'LineWidth', 2)
        hold off    
    else
        error('Failed to set axis type')
    end
    leg = legend('T1 effective', 'RF compensated', 'constant', 'total SNR', 'Fisher information'); 
    xlabel(xaxis_labels(k))
    ylabel('RMS error (average of 25 trials)')
    set(leg,'FontSize',20);
    set(gca,'FontSize',20);
    set(gca,'XTick', parameter_values(k, :)); 
    axis(axis_limits(k, :))
    tightfig(gcf); 
    print(gcf, '-dpdf', ['robustness_' num2str(k) '_line.pdf']);
end

end


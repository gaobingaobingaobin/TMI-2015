function [  ] = plot_bar_graphs(parameter_values, error_opt_array, error_const_array, error_RF_compensated_array, error_T1_effective_array, berkeley_colors, xaxis_labels, axis_limits)

bc = berkeley_colors([2 4 1 3], :); 
for k=1:size(parameter_values, 1) 
    
    figure  
    set(gca,'ColorOrder', bc, 'NextPlot', 'replacechildren')
    h_bar_graph = bar(parameter_values(k, :),  ...
        [error_T1_effective_array(k, :)./error_opt_array(k, :); 
        error_RF_compensated_array(k, :)./error_opt_array(k, :); 
        error_const_array(k, :)./error_opt_array(k, :); 
        ones(size(error_opt_array(k, :)))]'); 
    for i = 1:length(h_bar_graph)
        set(h_bar_graph(i), 'FaceColor', bc(i,:)) 
    end
    ylabel('normalized RMS error') 
    xlabel(xaxis_labels(k))
    box on
    leg = legend('T1 effective', 'RF compensated', 'constant', 'Fisher information');
    set(leg,'FontSize',20);
    set(gca,'FontSize',20);
    set(gca,'XTick', parameter_values(k, :)); 
    axis(axis_limits(k, :))
    tightfig(gcf); 
    print(gcf, '-dpdf', ['robustness_' num2str(k) '_bar.pdf']);

end

end


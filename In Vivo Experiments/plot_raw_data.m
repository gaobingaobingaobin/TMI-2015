function [fig1, fig2] = plot_raw_data( pyruvate, lactate, pyruvate_normalization, lactate_normalization )

% normalize pyruvate 
pa_norm_2nd_ex = pyruvate*pyruvate_normalization; 

% plot pyruvate 
fig1 = figure; 
for i=1:size(pa_norm_2nd_ex, 3)
    subplot(5, 6, i)
    image(pa_norm_2nd_ex(:, :, i, 1))
    axis off
    title(['t=' num2str(2*i)])
end

% normalize lactate 
lac_norm_2nd_ex = lactate*lactate_normalization;  

% plot lactate
fig2 = figure;  
for i=1:size(lac_norm_2nd_ex, 3)
    subplot(5, 6, i)
    image(lac_norm_2nd_ex(:, :, i, 1))
    axis off
    title(['t=' num2str(2*i)])
end

end


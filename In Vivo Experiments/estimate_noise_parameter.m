for exp=1:3
    load(['exp' num2str(exp) '_raw_data.mat'])
    pa_sample = pa(1, :, :, 1); 
    lac_sample = lac(1, :, :, 1);
    y = pa_sample(:); 
    N = length(y); 
    sigma_squared_hat_pa = 1/(2*N)*sum(y.^2) 
    
    y = lac_sample(:); 
    N = length(y); 
    sigma_squared_hat_lac = 1/(2*N)*sum(y.^2) 
    
end
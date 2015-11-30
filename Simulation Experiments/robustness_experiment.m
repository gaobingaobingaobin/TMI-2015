function [ error_opt_array, error_const_array, error_RF_compensated_array, error_T1_effective_array ] = robustness_experiment( model, parameters_to_vary, parameter_value_array, num_trials, thetas_opt, thetas_const, thetas_RF_compensated, thetas_T1_effective, kTRANS_val, kPL_val, R1P_val, R1L_val, u_est, sigma_squared )
% ROBUSTNESS_EXPERIMENT is a wrapper function for test_robustness that 
%    loops through all parameters in the vector parameters_to_vary 
%    and stores the results of test_robustness in an array 


error_opt_array = zeros(size(parameter_value_array)); 
error_const_array = zeros(size(parameter_value_array)); 
error_RF_compensated_array = zeros(size(parameter_value_array)); 
error_T1_effective_array = zeros(size(parameter_value_array)); 
for k=1:length(parameters_to_vary)
    parameter_to_vary = parameters_to_vary(k); 
    parameter_values = parameter_value_array(k, :); 
    [error_opt, error_const, ...
        error_RF_compensated, error_T1_effective ] = test_robustness(model, ...
        parameter_to_vary, parameter_values , num_trials, thetas_opt, thetas_const, ...
        thetas_RF_compensated, thetas_T1_effective, kTRANS_val, kPL_val, R1P_val, ...
        R1L_val, u_est, sigma_squared); 
    error_opt_array(k, :)            = error_opt; 
    error_const_array(k, :)          = error_const; 
    error_RF_compensated_array(k, :) = error_RF_compensated; 
    error_T1_effective_array(k, :)   = error_T1_effective; 
end

end


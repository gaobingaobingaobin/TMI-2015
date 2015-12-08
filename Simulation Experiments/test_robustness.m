function [ kPL_error_opt, kPL_error_const, kPL_error_RF_compensated, kPL_error_T1_effective, kPL_error_SNR] = test_robustness(model, parameter_to_vary, parameter_values, num_trials, thetas_opt, thetas_const, thetas_RF_compensated, thetas_T1_effective, thetas_SNR, kTRANS_val, kPL_val, R1P_val, R1L_val, u_est, sigma_squared )
% TEST_ROBUSTNESS loops through values of parameter_to_vary, 
%   computing the estimation error for each of three flip angle sequences
%   assuming that parameter_values(i) is the true value of parameter


    [ parameters_of_interest_est_opt, ...
    parameters_of_interest_est_const, ...
    parameters_of_interest_est_RF_compensated, ...
    parameters_of_interest_est_T1_effective, ...
    parameters_of_interest_est_SNR] ...
         = test_accuracy(model, parameter_to_vary, parameter_values, ...
         num_trials, thetas_opt, thetas_const, thetas_RF_compensated, ...
         thetas_T1_effective, thetas_SNR, ...
         kTRANS_val, kPL_val, R1P_val, R1L_val, u_est, sigma_squared ); 
     
    syms kPL 
    for param_count = 1:length(parameter_values)
        if isAlways(parameter_to_vary == kPL)
            kPL_val = parameter_values(param_count); 
        end
        kPL_error_opt(param_count)               = sqrt(mean(abs(parameters_of_interest_est_opt(:, 2, param_count)             - kPL_val).^2)); 
        kPL_error_const(param_count)             = sqrt(mean(abs(parameters_of_interest_est_const(:, 2, param_count)           - kPL_val).^2)); 
        kPL_error_RF_compensated(param_count)    = sqrt(mean(abs(parameters_of_interest_est_RF_compensated(:, 2, param_count)  - kPL_val).^2)); 
        kPL_error_T1_effective(param_count)      = sqrt(mean(abs(parameters_of_interest_est_T1_effective(:, 2, param_count)    - kPL_val).^2)); 
        kPL_error_SNR(param_count)               = sqrt(mean(abs(parameters_of_interest_est_SNR(:, 2, param_count)             - kPL_val).^2)); 
    end

end


function [ parameters_of_interest_est_opt, ...
    parameters_of_interest_est_const, ...
    parameters_of_interest_est_RF_compensated, ...
    parameters_of_interest_est_T1_effective, ...
    parameters_of_interest_est_SNR] ...
         = test_accuracy(model, parameter_to_vary, parameter_values, ...
         num_trials, thetas_opt, thetas_const, thetas_RF_compensated, ...
         thetas_T1_effective, thetas_SNR, ...
         kTRANS_val, kPL_val, R1P_val, R1L_val, u_est, sigma_squared )
% TEST_ROBUSTNESS loops through values of parameter_to_vary, 
%   computing the estimation error for each of three flip angle sequences
%   assuming that parameter_values(i) is the true value of parameter
    syms R1P R1L kPL kTRANS input_scale sigma_2 t0 B1

    %% Generate simulated data sets 

    y_thetas_opt            = zeros(model.n, model.N, num_trials, length(parameter_values)); 
    y_thetas_const          = zeros(model.n, model.N, num_trials, length(parameter_values)); 
    y_thetas_RF_compensated = zeros(model.n, model.N, num_trials, length(parameter_values)); 
    y_thetas_T1_effective   = zeros(model.n, model.N, num_trials, length(parameter_values)); 
    y_thetas_SNR            = zeros(model.n, model.N, num_trials, length(parameter_values)); 
    for param_count = 1:length(parameter_values)
        parameter_to_vary == sigma_2
        isAlways(parameter_to_vary == sigma_2)
        warning('off', 'symbolic:sym:isAlways:TruthUnknown')
        if isAlways(parameter_to_vary == R1P) || isAlways(parameter_to_vary == R1L)
            k = find(isAlways(eq(model.nuisance_parameters, parameter_to_vary)));  
            model.nuisance_parameters_nominal_values(k) = parameter_values(param_count); 
            model = discretize(model);  
        elseif isAlways(parameter_to_vary == kTRANS)
            kTRANS_val = parameter_values(param_count); 
            model.known_parameter_values = parameter_values(param_count); 
            model = discretize(model);  
        elseif isAlways(parameter_to_vary == kPL)
            kPL_val = parameter_values(param_count); 
            model.parameters_of_interest_nominal_values = parameter_values(param_count);
            model = discretize(model); 
        elseif isAlways(parameter_to_vary == input_scale) 
            input_params = [2.1430    3.4658/parameter_values(param_count)   10.4105*(parameter_values(param_count))^2.1430    3.2596];  % [gamma, beta, A0/1000, t0] 
            u_est = gamma_variate_input(input_params, 90/180*pi*ones(25, 1));  % function gamma_variate_input is meant to generate observed AIFs % calling it with a flip angle sequence of 90 degrees gives the true AIF
            if length(u_est) < model.N
                u_est = [u_est; zeros(model.N-length(u_est), 1)]; % pad end of u_est with zeros if model.N is greater than 25
            end
            model.nuisance_parameters_nominal_values = [R1P_val R1L_val u_est(1:end-1)' ];
            model = discretize(model); 
        elseif isAlways(parameter_to_vary == t0) 
            input_params = [2.1430    3.4658   10.4105  parameter_values(param_count)];  % [gamma, beta, A0/1000, t0] 
            u_est = gamma_variate_input(input_params, 90/180*pi*ones(25, 1));  % function gamma_variate_input is meant to generate observed AIFs % calling it with a flip angle sequence of 90 degrees gives the true AIF
            if length(u_est) < model.N
                u_est = [u_est; zeros(model.N-length(u_est), 1)]; % pad end of u_est with zeros if model.N is greater than 25
            end
            model.nuisance_parameters_nominal_values = [R1P_val R1L_val u_est(1:end-1)' ];
            model = discretize(model);     
        elseif isAlways(parameter_to_vary == sigma_2)
            model.noise_parameters = parameter_values(param_count)*[1 1];
            sigma_squared = parameter_values(param_count); 
            model = discretize(model); 
        elseif isAlways(parameter_to_vary == B1)
            model.flip_angle_input_matrix = parameter_values(param_count)*eye(2); 
            B1_val = parameter_values(param_count); 
            model = discretize(model); 
        else
            error('Parameter to vary not chosen successfully') 
        end
        warning('on', 'symbolic:sym:isAlways:TruthUnknown')
        var_initial = [ kPL_val, R1P_val, R1L_val, u_est', u_est', u_est', u_est', u_est']; 
        known_parameters = [kTRANS_val]; 
        for trial=1:num_trials
            y_thetas_opt(:, :, trial, param_count)   = generate_data(model, model.flip_angle_input_matrix*thetas_opt); 
            y_thetas_const(:, :, trial, param_count) = generate_data(model, model.flip_angle_input_matrix*thetas_const); 
            y_thetas_RF_compensated(:, :, trial, param_count) = generate_data(model, model.flip_angle_input_matrix*thetas_RF_compensated); 
            y_thetas_T1_effective(:, :, trial, param_count) = generate_data(model, model.flip_angle_input_matrix*thetas_T1_effective); 
            y_thetas_SNR(:, :, trial, param_count) = generate_data(model, model.flip_angle_input_matrix*thetas_SNR); 
        end
    % end

    
    %% Estimate nuisance parameters jointly from collected data sets  

    fact = 1
    options = optimset( 'Display', 'iter', 'MaxFunEvals', 5000); 
    % for param_count = 1:length(parameter_values) 
        
        % model.nuisance_parameters_nominal_values(k) = parameter_values(param_count); 
        var(1:3)
        % obj = @(var) joint_least_squares(var(1:3), known_parameters, sigma_squared, thetas_opt, thetas_const,  thetas_RF_compensated, mean(y_thetas_opt(:, :, :, param_count), 3), mean(y_thetas_const(:, :, :, param_count), 3), mean(y_thetas_RF_compensated(:, :, :, param_count), 3), var(4:3+model.N)/fact, var(4+model.N:3+2*model.N)/fact, var(4+model.N+model.N:3+model.N+model.N+model.N)/fact); 
        obj = @(var) joint_least_squares_5flipangles(var(1:3), known_parameters, sigma_squared, ...
            thetas_opt, thetas_const,  thetas_RF_compensated, thetas_T1_effective, thetas_SNR, ...
            mean(y_thetas_opt(:, :, :, param_count), 3), ...
            mean(y_thetas_const(:, :, :, param_count), 3), ...
            mean(y_thetas_RF_compensated(:, :, :, param_count), 3), ...
            mean(y_thetas_T1_effective(:, :, :, param_count), 3), ...
            mean(y_thetas_SNR(:, :, :, param_count), 3), ...
            var(4:3+model.N)/fact, ...
            var(4+model.N:3+2*model.N)/fact, ...
            var(4+model.N+model.N:3+model.N+model.N+model.N)/fact, ...
            var(4+model.N+model.N+model.N:3+model.N+model.N+model.N+model.N)/fact, ... 
            var(4+model.N+model.N+model.N+model.N:3+model.N+model.N+model.N+model.N+model.N)/fact); 
        joint_param_est(:, param_count) = fmincon(obj, var_initial, [], [], [], [], zeros(size(var_initial)), [1 1 1], [], options);
        param_of_interest(:, param_count)    = joint_param_est(1:3, param_count)
        u_est_opt(:, param_count)            = joint_param_est(4:3+model.N, param_count)
        u_est_const(:, param_count)          = joint_param_est(4+model.N:3+model.N+model.N, param_count)
        u_est_RF_compensated(:, param_count) = joint_param_est(4+model.N+model.N:3+model.N+model.N+model.N, param_count)
        u_est_T1_effective(:, param_count)   = joint_param_est(4+model.N+model.N+model.N:3+model.N+model.N+model.N+model.N, param_count)
        u_est_SNR(:, param_count)            = joint_param_est(4+model.N+model.N+model.N+model.N:3+model.N+model.N+model.N+model.N+model.N, param_count)
    % end


    %% Fit parameters of interest voxel-wise 

    %known_parameters = []; 
    % for param_count = 1:length(parameter_values)   
        
        
        
        var_initial = [kTRANS_val, kPL_val, R1P_val, R1L_val]; 
        for trial=1:num_trials

            [param_count, trial]

            obj = @(var) negative_log_likelihood_rician(var, [], sigma_squared, thetas_opt, y_thetas_opt(:, :, trial, param_count), abs(u_est_opt(:, param_count))); 
            options = optimset( 'Display', 'iter', 'MaxFunEvals', 1000); 
            parameters_of_interest_est_opt(trial, :, param_count) = fminunc(obj, var_initial, options);

            obj = @(var) negative_log_likelihood_rician(var, [], sigma_squared, thetas_const, y_thetas_const(:, :, trial, param_count), abs(u_est_const(:, param_count))); 
            options = optimset( 'Display', 'iter', 'MaxFunEvals', 1000); 
            parameters_of_interest_est_const(trial, :, param_count) = fminunc(obj, var_initial, options);

            obj = @(var) negative_log_likelihood_rician(var, [], sigma_squared, thetas_RF_compensated, y_thetas_RF_compensated(:, :, trial, param_count), abs(u_est_RF_compensated(:, param_count))); 
            options = optimset( 'Display', 'iter', 'MaxFunEvals', 1000); 
            parameters_of_interest_est_RF_compensated(trial, :, param_count) = fminunc(obj, var_initial, options);
            
            obj = @(var) negative_log_likelihood_rician(var, [], sigma_squared, thetas_T1_effective, y_thetas_T1_effective(:, :, trial, param_count), abs(u_est_T1_effective(:, param_count))); 
            options = optimset( 'Display', 'iter', 'MaxFunEvals', 1000); 
            parameters_of_interest_est_T1_effective(trial, :, param_count) = fminunc(obj, var_initial, options);
            
            obj = @(var) negative_log_likelihood_rician(var, [], sigma_squared, thetas_SNR, y_thetas_SNR(:, :, trial, param_count), abs(u_est_SNR(:, param_count))); 
            options = optimset( 'Display', 'iter', 'MaxFunEvals', 1000); 
            parameters_of_interest_est_SNR(trial, :, param_count) = fminunc(obj, var_initial, options);
        end
    %end
    
    %for param_count = 1:length(parameter_values)  
        kPL_error_opt(param_count)               = sqrt(mean(abs(parameters_of_interest_est_opt(:, 2, param_count)             - kPL_val).^2)); 
        kPL_error_const(param_count)             = sqrt(mean(abs(parameters_of_interest_est_const(:, 2, param_count)           - kPL_val).^2)); 
        kPL_error_RF_compensated(param_count)    = sqrt(mean(abs(parameters_of_interest_est_RF_compensated(:, 2, param_count)  - kPL_val).^2)); 
        kPL_error_T1_effective(param_count)      = sqrt(mean(abs(parameters_of_interest_est_T1_effective(:, 2, param_count)    - kPL_val).^2)); 
        kPL_error_SNR(param_count)               = sqrt(mean(abs(parameters_of_interest_est_SNR(:, 2, param_count)             - kPL_val).^2)); 
    end


end


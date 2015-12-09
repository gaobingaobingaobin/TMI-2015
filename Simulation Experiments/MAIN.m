% Script generates figures from paper "Optimizing flip angle sequences 
% for physiological parameter estimation in hyperpolarized carbon-13
% magnetic resonance imaging experiments"
%
% John Maidens
% July 2015 

clear all
close all
clc

% set seed for random number generation 
rng(42); 

% verify that required toolboxes are installed 
check_system_requirements(); 

% set colors 
berkeley_colors = ...
 1/256*[  0,   0,   0;
         45,  99, 127; 
        224, 158,  25; 
        194, 185, 167;
        217, 102, 31;
        185, 211, 182]; 
     

%% Specify system model 

% initialize model object 
model = linear_exchange_model; 

% define number of acquisitions 
model.N = 30; 

% define model parameters
syms R1P R1L kPL kTRANS 

% define input parameters 
U = sym('u', [1, model.N-1]); 

% parameter values
kTRANS_val = 0.0550; 
kPL_val = 0.0700; 
R1P_val = 1/20; 
R1L_val = 1/20; 

% parameters of interest 
% (those for which we wish to compute an estimate with minimal variance) 
model.parameters_of_interest = [ kPL ]; 
model.parameters_of_interest_nominal_values = [ kPL_val ]; 

% nuisance parameters
% (those parameters that are unknown but whose estimates we only care about
% insofar as they allow us to estamate the parameters of interest) 
model.nuisance_parameters = [R1P  R1L U];

input_params = [2.1430    3.4658   10.4105    3.2596];  % [gamma, beta, A0/1000, t0] 
u_est = gamma_variate_input(input_params, 90/180*pi*ones(25, 1));  % function gamma_variate_input is meant to generate observed AIFs % calling it with a flip angle sequence of 90 degrees gives the true AIF
u_est = [u_est; zeros(model.N-length(u_est), 1)]; % pad end of u_est with zeros if model.N is greater than 25
model.nuisance_parameters_nominal_values = [R1P_val R1L_val u_est(1:end-1)' ]; 

% known parameters
% (those whose values are assumed to be known constants) 
model.known_parameters       = [ kTRANS]; 
model.known_parameter_values = [ kTRANS_val ];  

% define system matrices for differential eq. 
%   dx/dt = A*x(t) + B*u(t)
%    y(t) = C*x(t) + D*u(t) 

% two-site exchange model with input feedthrough 
model.A = [ -kPL-R1P   0  ;
               kPL   -R1L];  
         
model.B = [kTRANS; 0]; 

model.C = [1 0; 
           0 1]; 
       
model.D = [0; 
           0]; 

% define input function shape  
model.u = [U 0]; 

% define initial condition 
model.x0 = [0; 0]; 

% define repetition time
model.TR = 2; 

% choose noise type
model.noise_type = 'Rician';
% model.noise_type = 'None';

% choose noise magnitude  
sigma_2_star = 2.3608e+04; 
model.noise_parameters = sigma_2_star*[1 1]; % sigma^2 values for the noise 

% choose flip angle input matrix 
%   This allows you to set linear equality constraints on the flip angles
%   for example setting: 
%
%      model.flip_angle_input_matrix = [1 0; 
%                                       0 1; 
%                                       1 0]; 
%
%   fixes the first and third flip angles to be equal one another. 
%   Consider defaulting to
% 
%      model.flip_angle_input_matrix = eye(model.n) 
% 
%   if you wish to compute all flip angles separately. 
model.flip_angle_input_matrix = eye(2); 
                             
% model.flip_angle_input_matrix = eye(model.m + model.n)                              

% choose design criterion 
design_criterion = 'D-optimal'; 
% design_criterion = 'E-optimal'; 
% design_criterion = 'A-optimal'; 
% design_criterion = 'T-optimal'; 
% design_criterion = 'totalSNR'; 

% discretize model (doing this in advance makes things run faster) 
model = discretize(model);  

% compute sensitivities (doing this in advance makes things run faster)
model = sensitivities(model);  


%% Plot simulated trajectories with constant flip angles 

thetas_const = 15*pi/180*ones(2, model.N); 

% generate simulated trajecories
[y, ~, x_true] = generate_data(model, thetas_const); 

% plot simulated input trajectories 
figure
set(gca,'ColorOrder', berkeley_colors(1:end, :), 'NextPlot', 'replacechildren')
plot(model.TR*(0:model.N-1), u_est, 'o-', 'LineWidth', 2)
title('Simulated input trajectory', 'FontSize', 20) 
xlabel('time (s)', 'FontSize', 20)
ylabel('u_t (au)', 'FontSize', 20)
set(gca,'FontSize',20)
tightfig(gcf);
print(gcf, '-dpdf', 'sim_input.pdf');

% plot simulated state trajectories 
figure
set(gca,'ColorOrder', berkeley_colors(2:end, :), 'NextPlot', 'replacechildren')
plot(model.TR*(0:model.N-1), x_true', 'o-', 'LineWidth', 2)
title('Simulated state trajectories', 'FontSize', 20) 
xlabel('time (s)', 'FontSize', 20)
ylabel('x_t (au)', 'FontSize', 20)
leg = legend('pyruvate (x_{1t} )', 'lactate (x_{2t} )'); 
set(leg,'FontSize',20);
set(gca,'FontSize',20);
tightfig(gcf);
print(gcf, '-dpdf', 'sim_state.pdf');


% plot simulated output trajectories 
figure
set(gca,'ColorOrder', berkeley_colors(2:end, :), 'NextPlot', 'replacechildren')
plot(model.TR*(0:model.N-1), y(1:2, :)', 'o-', 'LineWidth', 2)
title('Simulated measurement trajectories', 'FontSize', 20) 
xlabel('time (s)', 'FontSize', 20)
ylabel('Y_t (au)', 'FontSize', 20)
leg = legend('pyruvate (Y_{1t} )', 'lactate (Y_{2t} )'); 
set(leg,'FontSize',20);
set(gca,'FontSize',20);
tightfig(gcf);
print(gcf, '-dpdf', 'sim_measurement.pdf');


%% Design optimal flip angles

% specify optimization start point and options for MATLAB optimization toolbox 
initial_q_value = 5*pi/180*ones(size(model.flip_angle_input_matrix, 2), model.N);
options = optimset('MaxFunEvals', 3000, 'MaxIter', 500, 'Display', 'iter');

% perform optimization 
[thetas_opt, ~, q_opt] = optimal_flip_angle_design_regularized(model, design_criterion, ...
    initial_q_value, 0.1, options); 


%% Plot optimal flip angles 

figure 
set(gca,'ColorOrder', berkeley_colors(2:end, :), 'NextPlot', 'replacechildren')
plot(q_opt'.*180./pi, 'x-', 'LineWidth', 2) 
title('Optimized flip angle sequence') 
xlabel('acquisition number')
ylabel('flip angle (degrees)')
legend('pyruvate', 'lactate')


%% Visualize state and output trajectories for optimized sequence

[y, y_true, x_true] = generate_data(model, thetas_opt); 

% plot simulated state trajectories 
figure
set(gca,'ColorOrder', berkeley_colors(2:end, :), 'NextPlot', 'replacechildren')
plot(model.TR*(0:model.N-1), x_true', 'o-', 'LineWidth', 2)
title('Simulated state trajectories', 'FontSize', 20) 
xlabel('time (s)', 'FontSize', 20)
ylabel('x_t (au)', 'FontSize', 20)
leg = legend('pyruvate (x_{1t})', 'lactate (x_{2t})'); 
set(leg,'FontSize',20);
set(gca,'FontSize',20);
tightfig(gcf);
print(gcf, '-dpdf', 'sim_optimized_state.pdf');


% plot simulated output trajectories 
figure
set(gca,'ColorOrder', berkeley_colors(2:end, :), 'NextPlot', 'replacechildren')
plot(model.TR*(0:model.N-1), y(1:2, :)', 'o-', 'LineWidth', 2)
title('Simulated measurement trajectories', 'FontSize', 20) 
xlabel('time (s)', 'FontSize', 20)
ylabel('Y_t (au)', 'FontSize', 20)
leg = legend('pyruvate (Y_{1t})', 'lactate (Y_{2t})'); 
set(leg,'FontSize',20);
set(gca,'FontSize',20);
tightfig(gcf);
print(gcf, '-dpdf', 'sim_optimized_measurement.pdf');


%% Perform accuracy experiment 

noise_vals = sort([logspace(3, 6, 5), sigma_2_star ]);
num_trials = 25;
syms sigma_2

% load other time-varying flip angle sequences 
load('thetas_RF_compensated.mat') 
load('thetas_T1_effective.mat') 
load('thetas_SNR.mat') 

[ parameters_of_interest_est_opt, ...
    parameters_of_interest_est_const, ...
    parameters_of_interest_est_RF_compensated, ...
    parameters_of_interest_est_T1_effective, ...
    parameters_of_interest_est_SNR] ...
         = test_accuracy(model, sigma_2, noise_vals, ...
         num_trials, thetas_opt, thetas_const, thetas_RF_compensated, ...
         thetas_T1_effective, thetas_SNR, ...
         kTRANS_val, kPL_val, R1P_val, R1L_val, u_est, sigma_2_star ); 


 %% Generate scatterplot of parameter estimates  

sigma_index = 3 % index where true noise value lies 

figure
set(gca,'ColorOrder', berkeley_colors([2 4 3], :), 'NextPlot', 'replacechildren')
plot(parameters_of_interest_est_const(:, 1, sigma_index), parameters_of_interest_est_const(:, 2, sigma_index), 'o', 'MarkerFaceColor', berkeley_colors(2, :) )
hold on
plot(parameters_of_interest_est_RF_compensated(:, 1, sigma_index), parameters_of_interest_est_RF_compensated(:, 2, sigma_index), 'o', 'MarkerFaceColor', berkeley_colors(4, :))
plot(parameters_of_interest_est_opt(:, 1, sigma_index), parameters_of_interest_est_opt(:, 2, sigma_index), 'o', 'MarkerFaceColor', berkeley_colors(3, :))
plot(kTRANS_val, kPL_val, 'kx', 'MarkerSize', 20, 'LineWidth', 4)
hold off
leg = legend('constant', 'RF compensated', 'Fisher information', 'ground truth'); 
set(leg,'FontSize',20);
set(gca,'FontSize',20);
xlabel('kTRANS')
ylabel('kPL')
tightfig(gcf)
print(gcf, '-dpdf', 'kPL_kTRANS_numerical_est.pdf');

figure
set(gca,'ColorOrder', berkeley_colors([2 4 3], :), 'NextPlot', 'replacechildren')
plot(parameters_of_interest_est_const(:, 3, sigma_index), parameters_of_interest_est_const(:, 4, sigma_index),  'o', 'MarkerFaceColor', berkeley_colors(2, :) )
hold on
plot(parameters_of_interest_est_RF_compensated(:, 3, sigma_index), parameters_of_interest_est_RF_compensated(:, 4, sigma_index),  'o', 'MarkerFaceColor', berkeley_colors(4, :) )
plot(parameters_of_interest_est_opt(:, 3, sigma_index), parameters_of_interest_est_opt(:, 4, sigma_index),  'o', 'MarkerFaceColor', berkeley_colors(3, :) )
plot(R1P_val, R1L_val, 'kx', 'MarkerSize', 20, 'LineWidth', 4)
hold off
leg = legend('constant', 'RF compensated', 'Fisher information', 'ground truth');
set(leg,'FontSize',20);
set(gca,'FontSize',20);
xlabel('R1P')
ylabel('R1L')
tightfig(gcf)
print(gcf, '-dpdf', 'R1P_R1L_numerical_est.pdf');
     

%% Generate error plots

error_plot(1, kTRANS_val, 'kTRANS', ...
        noise_vals([1:2 4:end]), berkeley_colors, ...
        parameters_of_interest_est_opt(:, :, [1:2 4:end]), ...
        parameters_of_interest_est_const(:, :, [1:2 4:end]), ...
        parameters_of_interest_est_RF_compensated(:, :, [1:2 4:end]), ...
        parameters_of_interest_est_T1_effective(:, :, [1:2 4:end]), ...
        parameters_of_interest_est_SNR(:, :, [1:2 4:end]))
    
error_plot(2, kPL_val, 'kPL', ...
        noise_vals([1:2 4:end]), berkeley_colors, ...
        parameters_of_interest_est_opt(:, :, [1:2 4:end]), ...
        parameters_of_interest_est_const(:, :, [1:2 4:end]), ...
        parameters_of_interest_est_RF_compensated(:, :, [1:2 4:end]), ...
        parameters_of_interest_est_T1_effective(:, :, [1:2 4:end]), ...
        parameters_of_interest_est_SNR(:, :, [1:2 4:end]))
    
error_plot(3, R1P_val, 'R1P', ...
        noise_vals([1:2 4:end]), berkeley_colors, ...
        parameters_of_interest_est_opt(:, :, [1:2 4:end]), ...
        parameters_of_interest_est_const(:, :, [1:2 4:end]), ...
        parameters_of_interest_est_RF_compensated(:, :, [1:2 4:end]), ...
        parameters_of_interest_est_T1_effective(:, :, [1:2 4:end]), ...
        parameters_of_interest_est_SNR(:, :, [1:2 4:end]))
    
error_plot(4, R1L_val, 'R1L', ...
        noise_vals([1:2 4:end]), berkeley_colors, ...
        parameters_of_interest_est_opt(:, :, [1:2 4:end]), ...
        parameters_of_interest_est_const(:, :, [1:2 4:end]), ...
        parameters_of_interest_est_RF_compensated(:, :, [1:2 4:end]), ...
        parameters_of_interest_est_T1_effective(:, :, [1:2 4:end]), ...
        parameters_of_interest_est_SNR(:, :, [1:2 4:end]))
    % 
%     
% %% Estimate model parameters from simulated data 
% 
% goodness_of_fit_criterion = 'maximum-likelihood'; 
% noise_vals = sort([logspace(3, 6, 5), sigma_2_star ])
% num_trials = 25
% 
% var_initial = [ kPL_val, R1P_val, R1L_val, u_est', u_est', u_est']; 
% known_parameters = [kTRANS_val]; 
% 
% % load other time-varying flip angle sequences 
% load('thetas_RF_compensated.mat') 
% load('thetas_T1_effective.mat') 
% load('thetas_SNR.mat') 
% 
% % generate simulated data sets 
% for noise = 1:length(noise_vals)
%     noise
%     model.noise_parameters = noise_vals(noise)*[1 1]; 
%     for trial=1:num_trials
%         y_thetas_opt(:, :, trial, noise)   = generate_data(model, model.flip_angle_input_matrix*thetas_opt); 
%         y_thetas_const(:, :, trial, noise) = generate_data(model, model.flip_angle_input_matrix*thetas_const); 
%         y_thetas_RF_compensated(:, :, trial, noise) = generate_data(model, model.flip_angle_input_matrix*thetas_RF_compensated); 
%     end
% end
% 
% 
% %% Estimate nuisance parameters jointly from collected data sets  
% 
% fact = 1
% options = optimset( 'Display', 'iter', 'MaxFunEvals', 10000); 
% for noise = 1:length(noise_vals)
%     noise
%     obj = @(var) joint_least_squares(var(1:3), known_parameters, noise_vals(noise), thetas_opt, thetas_const,  thetas_RF_compensated, mean(y_thetas_opt(:, :, :, noise), 3), mean(y_thetas_const(:, :, :, noise), 3), mean(y_thetas_RF_compensated(:, :, :, noise), 3), var(4:3+model.N)/fact, var(4+model.N:3+2*model.N)/fact, var(4+model.N+model.N:3+model.N+model.N+model.N)/fact); 
%     joint_param_est(:, noise) = fmincon(obj, var_initial, [], [], [], [], [0 0 0], [1 1 1], [], options);
%     param_of_interest(:, noise)    = joint_param_est(1:3, noise)
%     u_est_opt(:, noise)            = joint_param_est(4:3+model.N, noise)
%     u_est_const(:, noise)          = joint_param_est(4+model.N:3+model.N+model.N, noise)
%     u_est_RF_compensated(:, noise) = joint_param_est(4+model.N+model.N:3+model.N+model.N+model.N, noise)
% end
% 
% 
% %% Fit parameters of interest voxel-wise 
% 
% var_initial = [kTRANS_val, kPL_val, R1P_val, R1L_val]; 
% known_parameters = []; 
% for noise = 1:length(noise_vals)
%     for trial=1:num_trials
%         
%         [noise, trial]
%   
%         obj = @(var) negative_log_likelihood_rician(var, known_parameters, noise_vals(noise), thetas_opt, y_thetas_opt(:, :, trial, noise), u_est_opt(:, noise)); 
%         options = optimset( 'Display', 'iter', 'MaxFunEvals', 1000); 
%         parameters_of_interest_est_opt(trial, :, noise) = fminunc(obj, var_initial, options);
%         
%         obj = @(var) negative_log_likelihood_rician(var, known_parameters, noise_vals(noise), thetas_const, y_thetas_const(:, :, trial, noise), u_est_const(:, noise)); 
%         options = optimset( 'Display', 'iter', 'MaxFunEvals', 1000); 
%         parameters_of_interest_est_const(trial, :, noise) = fminunc(obj, var_initial, options);
%                
%         obj = @(var) negative_log_likelihood_rician(var, known_parameters, noise_vals(noise), thetas_RF_compensated, y_thetas_RF_compensated(:, :, trial, noise), u_est_const(:, noise)); 
%         options = optimset( 'Display', 'iter', 'MaxFunEvals', 1000); 
%         parameters_of_interest_est_RF_compensated(trial, :, noise) = fminunc(obj, var_initial, options);
%     end
% end
% 
% 
% %% Generate plots of parameter error 
% 
% for noise = 1:length(noise_vals)
%         kPL_error_opt(noise)               = sqrt(mean(abs(parameters_of_interest_est_opt(:, 2, noise)             - kPL_val).^2)); 
%         kPL_error_const(noise)             = sqrt(mean(abs(parameters_of_interest_est_const(:, 2, noise)           - kPL_val).^2)); 
%         kPL_error_RF_compensated(noise)    = sqrt(mean(abs(parameters_of_interest_est_RF_compensated(:, 2, noise)  - kPL_val).^2)); 
%         kTRANS_error_opt(noise)            = sqrt(mean(abs(parameters_of_interest_est_opt(:, 1, noise)             - kTRANS_val).^2)); 
%         kTRANS_error_const(noise)          = sqrt(mean(abs(parameters_of_interest_est_const(:, 1, noise)           - kTRANS_val).^2)); 
%         kTRANS_error_RF_compensated(noise) = sqrt(mean(abs(parameters_of_interest_est_RF_compensated(:, 1, noise)  - kTRANS_val).^2)); 
%         R1P_error_opt(noise)               = sqrt(mean(abs(parameters_of_interest_est_opt(:, 3, noise)             - R1P_val).^2)); 
%         R1P_error_const(noise)             = sqrt(mean(abs(parameters_of_interest_est_const(:, 3, noise)           - R1P_val).^2)); 
%         R1P_error_RF_compensated(noise)    = sqrt(mean(abs(parameters_of_interest_est_RF_compensated(:, 3, noise)  - R1P_val).^2)); 
%         R1L_error_opt(noise)               = sqrt(mean(abs(parameters_of_interest_est_opt(:, 4, noise)             - R1L_val).^2)); 
%         R1L_error_const(noise)             = sqrt(mean(abs(parameters_of_interest_est_const(:, 4, noise)           - R1L_val).^2)); 
%         R1L_error_RF_compensated(noise)    = sqrt(mean(abs(parameters_of_interest_est_RF_compensated(:, 4, noise)  - R1L_val).^2)); 
% end
% 
% kPL_error_opt               = kPL_error_opt([1:2 4:end]);
% kPL_error_const             = kPL_error_const([1:2 4:end]);
% kPL_error_RF_compensated    = kPL_error_RF_compensated([1:2 4:end]);
% kTRANS_error_opt            = kTRANS_error_opt([1:2 4:end]);
% kTRANS_error_const          = kTRANS_error_const([1:2 4:end]);
% kTRANS_error_RF_compensated = kTRANS_error_RF_compensated([1:2 4:end]);
% R1P_error_opt               = R1P_error_opt([1:2 4:end]);
% R1P_error_const             = R1P_error_const([1:2 4:end]);
% R1P_error_RF_compensated    = R1P_error_RF_compensated([1:2 4:end]);
% R1L_error_opt               = R1L_error_opt([1:2 4:end]);
% R1L_error_const             = R1L_error_const([1:2 4:end]);
% R1L_error_RF_compensated    = R1L_error_RF_compensated([1:2 4:end]);
% noise_vals_small = noise_vals([1:2 4:end]);
% 
% figure
% set(gca,'ColorOrder', berkeley_colors([2 4 3], :), 'NextPlot', 'replacechildren')
% loglog(noise_vals_small, kPL_error_const, 'o-', 'LineWidth', 2)
% hold on
% loglog(noise_vals_small, kPL_error_RF_compensated, 'o-', 'LineWidth', 2)
% loglog(noise_vals_small, kPL_error_opt, 'o-', 'LineWidth', 2)
% hold off
% legend('constant', 'RF compensated', 'Fisher information')
% xlabel('\sigma^2')
% ylabel('RMS error (average of 25 trials)')
% % title('kPL estimation error comparison')
% set(leg,'FontSize',20);
% set(gca,'FontSize',20);
% tightfig(gcf)
% print(gcf, '-dpdf', 'kPL_error.pdf');
% 
% figure
% set(gca,'ColorOrder', berkeley_colors([2 4 3], :), 'NextPlot', 'replacechildren')
% loglog(noise_vals_small, kTRANS_error_const, 'o-', 'LineWidth', 2)
% hold on
% loglog(noise_vals_small, kTRANS_error_RF_compensated, 'o-', 'LineWidth', 2)
% loglog(noise_vals_small, kTRANS_error_opt, 'o-', 'LineWidth', 2)
% hold off
% legend('constant', 'RF compensated', 'Fisher information')
% xlabel('\sigma^2')
% ylabel('RMS error (average of 25 trials)')
% title('kTRANS estimation error comparison')
% set(leg,'FontSize',20);
% set(gca,'FontSize',20);
% tightfig(gcf)
% print(gcf, '-dpdf', 'kTRANS_error.pdf');
% 
% figure
% set(gca,'ColorOrder', berkeley_colors([2 4 3], :), 'NextPlot', 'replacechildren')
% loglog(noise_vals_small, R1P_error_const, 'o-', 'LineWidth', 2)
% hold on
% loglog(noise_vals_small, R1P_error_RF_compensated, 'o-', 'LineWidth', 2)
% loglog(noise_vals_small, R1P_error_opt, 'o-', 'LineWidth', 2)
% hold off
% legend('constant', 'RF compensated', 'Fisher information')
% xlabel('\sigma^2')
% ylabel('RMS error (average of 25 trials)')
% title('R1P estimation error comparison')
% set(leg,'FontSize',20);
% set(gca,'FontSize',20);
% tightfig(gcf)
% print(gcf, '-dpdf', 'R1P_error.pdf');
% 
% figure
% set(gca,'ColorOrder', berkeley_colors([2 4 3], :), 'NextPlot', 'replacechildren')
% loglog(noise_vals_small, R1L_error_const, 'o-', 'LineWidth', 2)
% hold on
% loglog(noise_vals_small, R1L_error_RF_compensated, 'o-', 'LineWidth', 2)
% loglog(noise_vals_small, R1L_error_opt, 'o-', 'LineWidth', 2)
% hold off
% legend('constant', 'RF compensated', 'Fisher information')
% xlabel('\sigma^2')
% ylabel('RMS error (average of 25 trials)')
% title('R1L estimation error comparison')
% set(leg,'FontSize',20);
% set(gca,'FontSize',20);
% tightfig(gcf)
% print(gcf, '-dpdf', 'R1L_error.pdf');
% 
% 
% %% Generate scatterplot of data 
% 
% sigma_index = 3 % index where true noise value lies 
% 
% figure
% set(gca,'ColorOrder', berkeley_colors([2 4 3], :), 'NextPlot', 'replacechildren')
% plot(parameters_of_interest_est_const(:, 1, sigma_index), parameters_of_interest_est_const(:, 2, sigma_index), 'o', 'MarkerFaceColor', berkeley_colors(2, :) )
% hold on
% plot(parameters_of_interest_est_RF_compensated(:, 1, sigma_index), parameters_of_interest_est_RF_compensated(:, 2, sigma_index), 'o', 'MarkerFaceColor', berkeley_colors(4, :))
% plot(parameters_of_interest_est_opt(:, 1, sigma_index), parameters_of_interest_est_opt(:, 2, sigma_index), 'o', 'MarkerFaceColor', berkeley_colors(3, :))
% plot(kTRANS_val, kPL_val, 'kx', 'MarkerSize', 20, 'LineWidth', 4)
% hold off
% leg = legend('constant', 'RF compensated', 'Fisher information', 'ground truth'); 
% set(leg,'FontSize',20);
% set(gca,'FontSize',20);
% xlabel('kTRANS')
% ylabel('kPL')
% tightfig(gcf)
% print(gcf, '-dpdf', 'kPL_kTRANS_numerical_est.pdf');
% 
% figure
% set(gca,'ColorOrder', berkeley_colors([2 4 3], :), 'NextPlot', 'replacechildren')
% plot(parameters_of_interest_est_const(:, 3, sigma_index), parameters_of_interest_est_const(:, 4, sigma_index),  'o', 'MarkerFaceColor', berkeley_colors(2, :) )
% hold on
% plot(parameters_of_interest_est_RF_compensated(:, 3, sigma_index), parameters_of_interest_est_RF_compensated(:, 4, sigma_index),  'o', 'MarkerFaceColor', berkeley_colors(4, :) )
% plot(parameters_of_interest_est_opt(:, 3, sigma_index), parameters_of_interest_est_opt(:, 4, sigma_index),  'o', 'MarkerFaceColor', berkeley_colors(3, :) )
% plot(R1P_val, R1L_val, 'kx', 'MarkerSize', 20, 'LineWidth', 4)
% hold off
% leg = legend('constant', 'RF compensated', 'Fisher information', 'ground truth');
% set(leg,'FontSize',20);
% set(gca,'FontSize',20);
% xlabel('R1P')
% ylabel('R1L')
% tightfig(gcf)
% print(gcf, '-dpdf', 'R1P_R1L_numerical_est.pdf');
% 
% 
% %% Plot bar graph of estimated error 
% 
% figure(13) 
% bc = berkeley_colors([2 4 3], :); 
% set(gca,'ColorOrder', bc, 'NextPlot', 'replacechildren')
% h_bar_graph = bar(log10(noise_vals_small),  ...
%     [kPL_error_const./kPL_error_opt; 
%     kPL_error_RF_compensated./kPL_error_opt; 
%     ones(size(kPL_error_opt))]'); 
% for i = 1:length(h_bar_graph)
%     set(h_bar_graph(i), 'FaceColor', bc(i, :)) 
% end
% ylabel('normalized RMS error') 
% xlabel('log_{10} \sigma^2') 
% legend('constant', 'RF compensated', 'Fisher information')
% box on
% set(leg,'FontSize',20);
% set(gca,'FontSize',20);
% axis([2.5 6.5 0 5])
% set(gca,'XTick', log10(noise_vals_small))
% tightfig(gcf)
% print(gcf, '-dpdf', 'kPL_error_bar.pdf');


%% Perform robustness experiment 

syms input_scale t0 B1
parameters_to_vary = [kTRANS, kPL, R1P, R1L, input_scale, B1];  
parameter_values = [linspace(0.01, 0.09, 5); 
                    linspace(0.03, 0.11, 5);
                    linspace(0.02, 0.08, 5);
                    linspace(0.02, 0.08, 5);
                    linspace(0.6, 1.4 , 5);
                    linspace(0.8, 1.2, 5)]; 


[ error_opt_array, error_const_array, error_RF_compensated_array, ...
    error_T1_effective_array, error_SNR_array ] = ...
    robustness_experiment(model, ...
    parameters_to_vary, parameter_values , 25, thetas_opt, thetas_const, ...
    thetas_RF_compensated, thetas_T1_effective, thetas_SNR, ...
    kTRANS_val, kPL_val, R1P_val, R1L_val, u_est, sigma_2_star); 


%% Plot the results of robustness experiment 

xaxis_labels = {'k_{TRANS}', 'k_{PL}', 'R_{1P}', 'R_{1L}', '\kappa', 'relative B_1'}; 
axis_limits_line  = [0.00 0.10 1e-03 1e02; 
                     0.02 0.12 0 0.07;
                     0.01 0.09 0 0.07;
                     0.01 0.09 0 0.07;
                     0.5  1.5  0 0.07;
                     0.75 1.25 0 0.07]; 
axis_limits_bar  = [0.00 0.10 0 5; 
                0.02 0.12 0 5;
                0.01 0.09 0 5;
                0.01 0.09 0 5;
                0.5  1.5  0 5;
                0.75 1.25 0 5]; 
axis_type = {'ylog', 'linear', 'linear', 'linear', 'linear', 'linear'}; 
plot_line_graphs(parameter_values, error_opt_array, error_const_array, ...
    error_RF_compensated_array, error_T1_effective_array, error_SNR_array,  ...
    berkeley_colors, xaxis_labels, axis_limits_line, axis_type)
plot_bar_graphs(parameter_values, error_opt_array, error_const_array, ...
    error_RF_compensated_array, error_T1_effective_array, error_SNR_array, ...
    berkeley_colors, xaxis_labels, axis_limits_bar)


%% Compute numerical value of kPL improvement

index = 2; 
for param_count = 1:length(noise_vals)
    error_opt(param_count)               = sqrt(mean(abs(parameters_of_interest_est_opt(:, index, param_count)            - kPL_val).^2)); 
    error_const(param_count)             = sqrt(mean(abs(parameters_of_interest_est_const(:, index, param_count)          - kPL_val).^2)); 
    error_RF_compensated(param_count)    = sqrt(mean(abs(parameters_of_interest_est_RF_compensated(:, index, param_count) - kPL_val).^2)); 
    error_T1_effective(param_count)      = sqrt(mean(abs(parameters_of_interest_est_T1_effective(:, index, param_count)   - kPL_val).^2)); 
    error_SNR(param_count)               = sqrt(mean(abs(parameters_of_interest_est_SNR(:, index, param_count)            - kPL_val).^2)); 
end
 
percent_improvement_over_constant_flip_angle_sequence       = mean(100*(error_const./error_opt - 1))
percent_improvement_over_RF_compensated_flip_angle_sequence = mean(100*(error_RF_compensated./error_opt - 1))
percent_improvement_over_T1_effective_flip_angle_sequence   = mean(100*(error_T1_effective./error_opt - 1))
percent_improvement_over_max_SNR_flip_angle_sequence        = mean(100*(error_SNR./error_opt - 1))



% Script generates figures from paper "Optimizing flip angle sequences 
% for physiological parameter estimation in hyperpolarized carbon-13
% magnetic resonance imaging experiments"
%
% John Maidens
% July 2015 

clear all 
close all
clc 

berkeley_colors = ...
 1/256*[  0,   0,   0;
         45,  99, 127; 
        224, 158,  25; 
        194, 185, 167;
        217, 102, 31;
        185, 211, 182]; 

    
%% Visualize proton data 

slice = 5; 
load('fsems_anatomy.mat', 'im_h')
figure
image(im_h(:, :, slice)/300) 
colormap(gray)
axis off
daspect([1 1 1])
tightfig(gcf); 
print(gcf, '-dpdf', 'proton_image.pdf');


%% Extract time series data from particular voxels from tumour region 

% choose location from which to extract time series 
x_tissue = 9; % coordinates for center of left kindney region 
y_tissue = 8; % coordinates for center of left kindney region
r_tissue = 2;  % radius for left kidney region
x_blood = 8;  % coordinates for center of blood vessel region
y_blood = 10;  % coordinates for center of blood vessel region
r_blood = 1.1;   % radius for blood vessel region
location_params = [x_tissue, y_tissue, r_tissue, x_blood, y_blood, r_blood];

% extract time series from experiment 1 data 
%   -- Fisher information optimized flip angles 
load('exp1_raw_data.mat') 
[ mean_AIF_exp1, mean_pyruvate_exp1, mean_lactate_exp1,  ...
      all_AIF_exp1,     all_pyruvate_exp1,     all_lactate_exp1 ] ...
    = extract_time_series( pa(:, :, :, 1), lac(:, :, :, 1), im_h, location_params);
alpha_exp1 = pi/180*flips; % flip angles used 

% extract time series from experiment 1 data 
%   -- Fisher information optimized flip angles 
load('exp2_raw_data.mat') 
[ mean_AIF_exp2, mean_pyruvate_exp2, mean_lactate_exp2,  ...
      all_AIF_exp2,     all_pyruvate_exp2,     all_lactate_exp2 ] ...
    = extract_time_series( pa(:, :, :, 1), lac(:, :, :, 1), im_h, location_params);
alpha_exp2 = pi/180*flips; % flip angles used 

% extract time series from experiment 1 data 
%   -- Fisher information optimized flip angles 
load('exp3_raw_data.mat') 
[ mean_AIF_exp3, mean_pyruvate_exp3, mean_lactate_exp3,  ...
      all_AIF_exp3,     all_pyruvate_exp3,     all_lactate_exp3 ] ...
    = extract_time_series( pa(:, :, :, 1), lac(:, :, :, 1), im_h, location_params);
alpha_exp3 = flips; % flip angles used  


%% Plot flip angles used in experiments 

figure 
set(gca,'ColorOrder', berkeley_colors(2:3, :), 'NextPlot', 'replacechildren')
plot(180/pi*alpha_exp1', 'x-', 'LineWidth', 2) 
title('Optimized flip angle sequence') 
xlabel('acquisition number')
ylabel('flip angle (degrees)')
leg = legend('pyruvate', 'lactate'); 
axis([0 30 0 100])
set(leg,'FontSize',14);
set(gca,'FontSize',14);
tightfig(gcf);
print(gcf, '-dpdf', 'flip_angles_fisher.pdf');

figure 
set(gca,'ColorOrder', berkeley_colors(2:3, :), 'NextPlot', 'replacechildren')
plot(180/pi*alpha_exp2', 'x-', 'LineWidth', 2) 
title('T1 effective flip angle sequence') 
xlabel('acquisition number')
ylabel('flip angle (degrees)')
leg = legend('pyruvate', 'lactate'); 
axis([0 30 0 100])
set(leg,'FontSize',14);
set(gca,'FontSize',14);
tightfig(gcf);
print(gcf, '-dpdf', 'flip_angles_T1_effective.pdf');

figure 
set(gca,'ColorOrder', berkeley_colors(2:3, :), 'NextPlot', 'replacechildren')
plot(180/pi*alpha_exp3', 'x-', 'LineWidth', 2) 
title('RF compensated flip angle sequence') 
xlabel('acquisition number')
ylabel('flip angle (degrees)')
leg = legend('pyruvate', 'lactate'); 
axis([0 30 0 100])
set(leg,'FontSize',14);
set(gca,'FontSize',14);
tightfig(gcf);
print(gcf, '-dpdf', 'flip_angles_RF_compensated.pdf');


%% Visulaize raw data 

load('exp1_raw_data.mat') 
[fig1, fig2] = plot_raw_data( pa, lac, 1/100, 1/30 );
fig1 = tightfig(fig1); 
fig2 = tightfig(fig2); 
print(fig1, '-dpdf', 'raw_data_pyruvate.pdf');
print(fig2, '-dpdf', 'raw_data_lactate.pdf');


%% Plot slice profiles 

% load actual slice profile 
load('correction_factors_profiles.mat')
lo = find(zloc > -30, 1);
hi = find(zloc > 30, 1);
Pi = flip_profile(lo:hi); 
z = zloc(lo:hi); 

% generate ideal slice profile 
Pi_ideal = (z <= 7.5).*(z >= -7.5); 

% number of discretization points 
M = length(Pi); 

figure
set(gca,'ColorOrder', berkeley_colors(2:3, :), 'NextPlot', 'replacechildren')
plot(z*8/15, Pi_ideal, z*8/15, Pi, 'Linewidth', 2)
xlabel('distance from center of slice (mm)')
ylabel('fraction of flip angle applied') 
title('Flip angle profile') 
leg = legend('ideal profile', 'realistic profile') ; 
axis([-15 15 0 1.05])
set(leg,'FontSize',14);
set(gca,'FontSize',14);
tightfig(gcf); 
print(gcf, '-dpdf', 'profile.pdf');


%% Visualize locations where time series are extracted 

figure
image(im_h(:, :, 5)/300) 
colormap(gray)
axis off
hold on 
gold_im_h = cat(3, berkeley_colors(3, 1)*ones(size(im_h(:, :, 5))), berkeley_colors(3, 2)*ones(size(im_h(:, :, 5))), berkeley_colors(3, 3)*ones(size(im_h(:, :, 5))));
h = imshow(gold_im_h);
[X, Y] = meshgrid(1:size(pa, 1), 1:size(pa, 1)); 
mask_left = ((X-x_tissue).^2 + (Y-y_tissue).^2 < r_tissue^2); 
mask_left_im_h = kron(mask_left, [ones(15, 15) zeros(15, 1); zeros(1, 15) 0]); 
set(h, 'AlphaData', mask_left_im_h); 
hold off
daspect([1 1 1])
tightfig(gcf); 
print(gcf, '-dpdf', 'voxels_extracted.pdf');


%% Perform first step of parameter estimation process: fit a single AIF to all the voxels simultaneously 

TR = 2; 
N = 30; 

% nominal parameter values 
kTRANS = 0.05; 
kPL    = 0.05; 
R1P    = 1/20; 
R1L    = 1/20; 
noise_level = 4.2844e04;   % estimate of sigma squared for this data set

parameters_initial = [kTRANS kPL R1L]; 
known_parameters = [ R1P  ]; 

% compute mean state and input trajectories 
u_exp1    = mean_AIF_exp1./sum(sin(Pi*alpha_exp1(1, :)), 1); 
pa_mean_exp1   = mean(all_pyruvate_exp1, 1); 
lac_mean_exp1  = mean(all_lactate_exp1, 1); 
data_exp1 = [pa_mean_exp1; lac_mean_exp1]; 

u_exp2    = mean_AIF_exp2./sum(sin(Pi*alpha_exp2(1, :)), 1); 
pa_mean_exp2   = mean(all_pyruvate_exp2, 1); 
lac_mean_exp2  = mean(all_lactate_exp2, 1); 
data_exp2 = [pa_mean_exp2; lac_mean_exp2]; 

u_exp3    = mean_AIF_exp3./sum(sin(Pi*alpha_exp3(1, :)), 1);
pa_mean_exp3   = mean(all_pyruvate_exp3, 1); 
lac_mean_exp3  = mean(all_lactate_exp3, 1); 
data_exp3 = [pa_mean_exp3; lac_mean_exp3]; 

% initialization of the optimization decision variables 
t = 0:TR:(N-1)*TR;
alph = 2;
A0 = 5000;
beta = 5; 
fact = 1e-06 % normalizing factor for AIF so that decision variables are of similar scale
u_initial = fact*A0*t.^alph.*exp(-t/beta); 
var_initial = [parameters_initial, u_initial, u_initial, u_initial]; 

% perform optimization for joint parameter fit 
obj = @(var) joint_least_squares(var(1:3), known_parameters, noise_level, Pi, alpha_exp1, alpha_exp2, alpha_exp3, data_exp1, data_exp2, data_exp3, var(4:3+N)/fact, var(4+N:3+2*N)/fact, var(4+N+N:3+3*N)/fact); 
options = optimset( 'Display', 'iter', 'MaxFunEvals', 10000); 
param_est = fmincon(obj, var_initial, [], [], [], [], [0 0 0 ], [1 1 1], [], options) 

% extract resulting estimates of the arterial input function 
u_est_exp1 = param_est(4:3+N)/fact; 
u_est_exp2 = param_est(4+N:3+N+N)/fact; 
u_est_exp3 = param_est(4+N+N:3+N+N+N)/fact; 


%% Perform the second step of parameter estimation: fit kPL, kTRANS separately for each voxel 

% fixed parameter values  
R1P    = 1/20; 
R1L = param_est(3); 
known_parameters = [R1P R1L]; 

% inital parameter values for search 
kTRANS = param_est(1); 
kPL = param_est(2); 
parameters_initial = [kTRANS kPL ]; 

% Normalize magnitude of the input functions 
scaling_factor = sum([u_exp1  u_exp2 u_exp3])/sum(param_est(3:end))
% scaling_factor = 1/fact

% Do ML estimation 
var_initial = parameters_initial; 
for i=1:size(all_pyruvate_exp1, 1)
    i
    
    data_exp1 = abs([all_pyruvate_exp1(i, :); all_lactate_exp1(i, :)]); 
    obj = @(var) negative_log_likelihood_rician(var, known_parameters, noise_level, Pi, alpha_exp1, data_exp1, scaling_factor*abs(param_est(3:2+N))); 
 %   obj = @(var) negative_log_likelihood_rician(var, known_parameters, noise_level, Pi, alpha_exp1, data_exp1, u_est_exp1); 
    options = optimset( 'Display', 'iter', 'MaxFunEvals', 1000); 
    param_est_exp1(i, :) = fminunc(obj, var_initial, options);

    data_exp2 = abs([all_pyruvate_exp2(i, :); all_lactate_exp2(i, :)]); 
    obj = @(var) negative_log_likelihood_rician(var, known_parameters, noise_level, Pi, alpha_exp2, data_exp2, scaling_factor*abs(param_est(3+N:2+N+N))); 
%    obj = @(var) negative_log_likelihood_rician(var, known_parameters, noise_level, Pi, alpha_exp2, data_exp2, u_est_exp2); 
    options = optimset( 'Display', 'iter', 'MaxFunEvals', 1000); 
    param_est_exp2(i, :) = fminunc(obj, var_initial, options);

    data_exp3 = abs([all_pyruvate_exp3(i, :); all_lactate_exp3(i, :)]); 
    obj = @(var) negative_log_likelihood_rician(var, known_parameters, noise_level, Pi, alpha_exp3, data_exp3, scaling_factor*abs(param_est(3+N+N:2+N+N+N))); 
%    obj = @(var) negative_log_likelihood_rician(var, known_parameters, noise_level, Pi, alpha_exp3, data_exp3, u_est_exp3); 
    options = optimset( 'Display', 'iter', 'MaxFunEvals', 1000); 
    param_est_exp3(i, :) = fminunc(obj, var_initial, options);
end


%% plot estimated input functions 

figure 
set(gca,'ColorOrder', berkeley_colors([2 4 3], :), 'NextPlot', 'replacechildren')
plot(0:TR:(N-1)*TR,  scaling_factor*param_est(3+N:2+N+N),  0:TR:(N-1)*TR,  scaling_factor*param_est(3+N+N:2+N+N+N), 0:TR:(N-1)*TR,  scaling_factor*param_est(3:2+N), 'LineWidth', 2)
xlabel('t (seconds)')
ylabel('Estimated u(t)')
leg = legend('T1 effective', 'RF compensated', 'Fisher information');  
title('Inputs estimated from state trajectories') 
set(leg,'FontSize',14);
set(gca,'FontSize',14);
tightfig(gcf);
print(gcf, '-dpdf', 'estimated_inputs.pdf');


%% Visulalize spread of kTRANS and kPL estimates 

figure 
set(gca,'ColorOrder', berkeley_colors([2 4 3], :), 'NextPlot', 'replacechildren')
hold on
plot(param_est_exp2(:, 1), param_est_exp2(:, 2), 'o', 'MarkerFaceColor', berkeley_colors(2, :))
plot(param_est_exp3(:, 1), param_est_exp3(:, 2), 'o', 'MarkerFaceColor', berkeley_colors(4, :)) 
plot(param_est_exp1(:, 1), param_est_exp1(:, 2), 'o', 'MarkerFaceColor', berkeley_colors(3, :))
hold off
xlabel('k_{TRANS}')
ylabel('k_{PL}')
leg = legend('T1 effective', 'RF compensated', 'Fisher information');  
title('Spread of parameter estimates for multiple voxels')
set(leg,'FontSize',14);
set(gca,'FontSize',14);
xlabel('kTRANS')
ylabel('kPL')
tightfig(gcf);
print(gcf, '-dpdf', 'kPL_kTRANS_in_vivo_est.pdf');


%% Plot fit of model to data 

i = 6 % choose voxel 

figure 
set(gca,'ColorOrder', berkeley_colors([2 3 1 4], :), 'NextPlot', 'replacechildren')
data_exp1 = abs([all_pyruvate_exp1(i, :); all_lactate_exp1(i, :)]); 
y_exp1 = trajectories_2param(param_est_exp1(i, :), known_parameters, Pi, alpha_exp1, scaling_factor*param_est(3:2+N));
plot(0:TR:TR*(N-1), data_exp1, '-o', 0:TR:TR*(N-1), y_exp1, '-', 'LineWidth', 2) 
xlabel('time (s)')
ylabel('measured signal (au)') 
title('Output fit -- Fisher information-optimized flip angles')
leg = legend('pyruvate (data)', 'lactate (data)', 'pyruvate (model)', 'lactate (model)');  
axis([0 60 0 10000])
set(leg,'FontSize',14);
set(gca,'FontSize',14);
tightfig(gcf);
print(gcf, '-dpdf', 'fit_Fisher_information.pdf');

figure
set(gca,'ColorOrder', berkeley_colors([2 3 1 4], :), 'NextPlot', 'replacechildren')
data_exp2 = abs([all_pyruvate_exp2(i, :); all_lactate_exp2(i, :)]); 
y_exp2 = trajectories_2param(param_est_exp2(i, :), known_parameters, Pi, alpha_exp2, scaling_factor*param_est(3+N:2+N+N));
plot(0:TR:TR*(N-1), data_exp2, '-o', 0:TR:TR*(N-1), y_exp2, '-', 'LineWidth', 2) 
xlabel('time (s)')
ylabel('measured signal (au)') 
title('Output fit -- T1 effective flip angles')
leg = legend('pyruvate (data)', 'lactate (data)', 'pyruvate (model)', 'lactate (model)');  
axis([0 60 0 10000])
set(leg,'FontSize',14);
set(gca,'FontSize',14);
tightfig(gcf);
print(gcf, '-dpdf', 'fit_T1_effective.pdf');

figure
set(gca,'ColorOrder', berkeley_colors([2 3 1 4], :), 'NextPlot', 'replacechildren')
data_exp3 = abs([all_pyruvate_exp3(i, :); all_lactate_exp3(i, :)]); 
y_exp3 = trajectories_2param(param_est_exp3(i, :), known_parameters, Pi, alpha_exp3, scaling_factor*param_est(3+N+N:2+N+N+N));
plot(0:TR:TR*(N-1), data_exp3, '-o', 0:TR:TR*(N-1), y_exp3, '-', 'LineWidth', 2) 
xlabel('time (s)')
ylabel('measured signal (au)') 
title('Output fit -- RF compensated flip angles')
leg = legend('pyruvate (data)', 'lactate (data)', 'pyruvate (model)', 'lactate (model)');  
axis([0 60 0 10000])
set(leg,'FontSize',14);
set(gca,'FontSize',14);
tightfig(gcf);
print(gcf, '-dpdf', 'fit_RF_compensated.pdf');


%% Do ML estimation for Fisher information data

load('exp1_raw_data.mat')

var_initial = parameters_initial; 
for i=1:size(pa, 1)
    for j=1:size(pa, 2)

        [i, j] 

        data = abs([pa(i, j, :, 1); lac(i, j, :, 1)]); 
        obj = @(var) negative_log_likelihood_rician(var, known_parameters, noise_level, Pi, alpha_exp1, data, scaling_factor*abs(param_est(3:2+N))); 
        options = optimset( 'Display', 'iter', 'MaxFunEvals', 1000); 
        parameter_map_Fisher_information(i, j, :) = fminunc(obj, var_initial, options);

    end
end


%% Plot parameter maps for Fisher information data

% plot kTRANS map 
figure
imagesc(parameter_map_Fisher_information(:, :, 1))
title('Map of parameter kTRANS')
axis off

% plot kPL map superimposed over proton image 
proton_im(:, :, 1) = im_h(:, :, 5)/20000; 
proton_im(:, :, 2) = im_h(:, :, 5)/20000; 
proton_im(:, :, 3) = im_h(:, :, 5)/20000; 

figure 
image(proton_im) 
hold on
% transparency is set using kTRANS value 
image(imresize(parameter_map_Fisher_information(:, :, 2), [256, 256]), 'AlphaData', imresize(parameter_map_Fisher_information(:, :, 1), [256, 256])*10, 'CDataMapping','scaled')
caxis([0.09 0.13])
colorbar
hold off
axis off
daspect([1 1 1])
print(gcf, '-dpdf', 'parameter_map_Fisher_information.pdf');


%% Do ML estimation for T1 effective data 

load('exp2_raw_data.mat')

var_initial = parameters_initial; 
for i=1:size(pa, 1)
    for j=1:size(pa, 2)

        [i, j] 

        data = abs([pa(i, j, :, 1); lac(i, j, :, 1)]); 
        obj = @(var) negative_log_likelihood_rician(var, known_parameters, noise_level, Pi, alpha_exp2, data, scaling_factor*abs(param_est(3+N:2+N+N))); 
        options = optimset( 'Display', 'iter', 'MaxFunEvals', 1000); 
        parameter_map_T1_effective(i, j, :) = fminunc(obj, var_initial, options);

    end
end


%% Plot parameter maps for T1 effective data

% plot kTRANS map 
figure
imagesc(parameter_map_T1_effective(:, :, 1))
title('Map of parameter kTRANS')
axis off

% plot kPL map superimposed over proton image 
figure 
image(proton_im) 
hold on
% transparency is set using kTRANS value 
image(imresize(parameter_map_T1_effective(:, :, 2), [256, 256]), 'AlphaData', imresize(parameter_map_T1_effective(:, :, 1), [256, 256])*10, 'CDataMapping','scaled')
caxis([0.09 0.13])
colorbar
hold off
axis off
daspect([1 1 1])
print(gcf, '-dpdf', 'parameter_map_T1_effective.pdf');


%% Do ML estimation for RF compensated data 

load('exp3_raw_data.mat')

var_initial = parameters_initial; 
for i=1:size(pa, 1)
    for j=1:size(pa, 2)

        [i, j] 

        data = abs([pa(i, j, :, 1); lac(i, j, :, 1)]); 
        obj = @(var) negative_log_likelihood_rician(var, known_parameters, noise_level, Pi, alpha_exp3, data, scaling_factor*abs(param_est(3+N+N:2+N+N+N))); 
        options = optimset( 'Display', 'iter', 'MaxFunEvals', 1000); 
        parameter_map_RF_compensated(i, j, :) = fminunc(obj, var_initial, options);

    end
end


%% Plot parameter maps for RF compensated data

% plot kTRANS map 
figure
imagesc(parameter_map_RF_compensated(:, :, 1))
title('Map of parameter kTRANS')
axis off

% plot kPL map superimposed over proton image 
figure 
image(proton_im) 
hold on
% transparency is set using kTRANS value 
image(imresize(parameter_map_RF_compensated(:, :, 2), [256, 256]), 'AlphaData', imresize(parameter_map_RF_compensated(:, :, 1), [256, 256])*10, 'CDataMapping','scaled')
caxis([0.09 0.13])
colorbar
hold off
axis off
daspect([1 1 1])
print(gcf, '-dpdf', 'parameter_map_RF_compensated.pdf');




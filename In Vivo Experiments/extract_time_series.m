function [  mean_AIF_15, mean_pyruvate_15, mean_lactate_15, all_AIF, all_pyruvate, all_lactate, left_kidney_15,  blood_15] = extract_time_series( pa, lac, proton, location_params)

x_left_15  = location_params(1); 
y_left_15  = location_params(2); 
r_left_15  = location_params(3); 
x_blood_15 = location_params(4); 
y_blood_15 = location_params(5); 
r_blood_15 = location_params(6); 

N = 30; 

[X, Y] = meshgrid(1:size(pa, 1), 1:size(pa, 1)); 
mask_left = ((X-x_left_15).^2 + (Y-y_left_15).^2 < r_left_15^2); 
mask_blood = ((X-x_blood_15).^2 + (Y-y_blood_15).^2 < r_blood_15^2); 


% generate array of left kidney and blood points 
left_kidney_15 = []; 
blood_15 = []; 
for i=1:size(pa, 1)
    for j=1:size(pa, 1)
        if mask_left(i, j)
            left_kidney_15 = [left_kidney_15; X(i,j), Y(i,j)]; 
        end
        if mask_blood(i, j)
            blood_15 = [blood_15; X(i, j), Y(i,j)]; 
        end
    end
end

 mean_pyruvate_15 = zeros(1, N); 
 mean_lactate_15 = zeros(1, N);
 mean_AIF_15 = zeros(1, N); 
 
 all_pyruvate = []; 
 all_lactate = []; 
 all_AIF = []; 
 
 % plot time series for each voxel in left kidney region 
 for i=1:min(length(left_kidney_15), length(blood_15))
     x = left_kidney_15(i, 1); 
     y = left_kidney_15(i, 2); 
     % get time series for that point 
    pyruvate = reshape(pa(y, x, :, 1), 1, length(pa(y, x, :, 1))); 
    lactate = reshape(lac(y, x, :, 1), 1, length(lac(y, x, :, 1))); 
    mean_pyruvate_15 = mean_pyruvate_15 + pyruvate; 
    mean_lactate_15 = mean_lactate_15 + lactate; 
    all_pyruvate = [all_pyruvate; pyruvate]; 
    all_lactate  = [all_lactate;  lactate ]; 
    

     x = blood_15(i, 1); 
     y = blood_15(i, 2); 
     % get time series for that point 
    AIF = reshape(pa(y, x, :, 1), 1, length(pa(y, x, :, 1))); 
    mean_AIF_15 = mean_AIF_15 + AIF; 
    all_AIF = [all_AIF; AIF]; 
 end
 
 % plot time series for each voxel in blood region 
 for i=min(length(left_kidney_15), length(blood_15))+1:length(left_kidney_15)
     x = left_kidney_15(i, 1); 
     y = left_kidney_15(i, 2); 
     % get time series for that point 
    pyruvate = reshape(pa(y, x, :, 1), 1, length(pa(y, x, :, 1))); 
    lactate = reshape(lac(y, x, :, 1), 1, length(lac(y, x, :, 1))); 
    mean_pyruvate_15 = mean_pyruvate_15 + pyruvate; 
    mean_lactate_15 = mean_lactate_15 + lactate; 
    all_pyruvate = [all_pyruvate; pyruvate]; 
    all_lactate  = [all_lactate;  lactate ]; 
    
 end

 mean_lactate_15 = mean_lactate_15/length(left_kidney_15); 
 mean_pyruvate_15 = mean_pyruvate_15/length(left_kidney_15); 
 mean_AIF_15 = mean_AIF_15/length(blood_15); 
 
end


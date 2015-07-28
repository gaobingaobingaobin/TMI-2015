function [ l ] = joint_least_squares(parameters, known_parameters, noise_level, Pi, alpha1, alpha2, alpha3, z1, z2, z3, u1, u2, u3)
    
    % compute trajectory of the model with parameter values 
    y1 = trajectories(parameters, known_parameters, Pi, alpha1, u1)';
    y2 = trajectories(parameters, known_parameters, Pi, alpha2, u2)';
    y3 = trajectories(parameters, known_parameters, Pi, alpha3, u3)';
    
    err = [y1-z1; y2-z2; y3-z3]; 
    l = norm(err(:)); 

end

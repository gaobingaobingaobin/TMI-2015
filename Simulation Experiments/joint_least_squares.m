function [ l ] = joint_least_squares(parameters, known_parameters, noise_level, alpha1, alpha2, alpha3, z1, z2, z3, u1, u2, u3)

    % compute trajectory of the model with parameter values 
    y1 = trajectories_joint(parameters, known_parameters, alpha1, u1);
    y2 = trajectories_joint(parameters, known_parameters, alpha2, u2);
    y3 = trajectories_joint(parameters, known_parameters, alpha3, u3);
    
    err = [y1-z1; y2-z2; y3-z3]; 
 
    l = norm(err(:)); 

end

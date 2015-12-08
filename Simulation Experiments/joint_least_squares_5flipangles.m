function [ l ] = joint_least_squares_5flipangles(parameters, known_parameters, noise_level, alpha1, alpha2, alpha3, alpha4, alpha5, z1, z2, z3, z4, z5, u1, u2, u3, u4, u5)

    % compute trajectory of the model with parameter values 
    y1 = trajectories_joint(parameters, known_parameters, alpha1, u1);
    y2 = trajectories_joint(parameters, known_parameters, alpha2, u2);
    y3 = trajectories_joint(parameters, known_parameters, alpha3, u3);
    y4 = trajectories_joint(parameters, known_parameters, alpha4, u4);
    y5 = trajectories_joint(parameters, known_parameters, alpha5, u5);
    
    err = [y1-z1; y2-z2; y3-z3; y4-z4; y5-z5]; 
 
    l = norm(err(:)); 

end

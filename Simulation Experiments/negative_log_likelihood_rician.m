function [ l1 ] = negative_log_likelihood_rician(parameters, known_parameters, noise_level, alpha1, z1, u1)
    %FUNCTION NEGATIVE_LOG_LIKELIHOOD_RICIAN Computes log likelihood for 
    %    compartmental model with Rician noise 
    
    N = length(alpha1); 

    % compute trajectory of the model with parameter values 
    y1 = trajectories_4param(parameters, known_parameters, alpha1, u1);

    % compute negative log likelihood 
    l1 = 0;
    for t = 1:N
        for k = 1:2
            l1 = l1 - (...
                log(z1(k, t)) - log(noise_level) ...
                - (z1(k, t)^2 + y1(k, t)^2)/(2*noise_level) ...
                + z1(k, t)*y1(k, t)/noise_level ...
                + log(besseli(0, z1(k, t)*y1(k, t)/noise_level, 1))...
            ); 
        
        end
    end

end

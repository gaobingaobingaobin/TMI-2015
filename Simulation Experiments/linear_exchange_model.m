classdef linear_exchange_model
    %CLASS LINEAR_EXCHANGE_MODEL describes multi-site model of exchange, relaxation, perfusion, etc. 
    %   Allows model description to easily be passed to functions 
    %   For use with Flip Angle Design scripts
    %   
    %   Flip Angle Design Toolbox 
    %   John Maidens (maidens@eecs.berkeley.edu) 
    %   June 2014 
    
    properties
        A     % continuous-time state dynamics matrix 
        B     % continuous-time input matrix 
        C     % output matrix
        D     % feedthrough matrix 
        u     % input function 
        x0    % initial condition 
        TR    % repetition time 
        N     % number of acquisitions 
        n     % dimension of state
        ni    % dimension of input 
        no    % dimension of the output 
        
        parameters_of_interest                % vector of symbolic variables 
        parameters_of_interest_nominal_values % vector of corresponding nominal values
        nuisance_parameters                   % vector of symbolic variables
        nuisance_parameters_nominal_values    % vector of corresponding nominal values
        known_parameters                      % vector of symbolic variables 
        known_parameter_values                % vector of corresponding nominal values
        noise_type                            % string selection noise type
        noise_parameters                      % parameters defining noise shape 
        flip_angle_input_matrix               % matrix specifying linear transformation of flip angles 
        nq                                    % number of independent flip angles to be chosen at each time 
        
        discretized = false % flag to determine whether discretization has occured 
        Ad_sym         % discretized dynamics matrix (symbolic)
        Bd_sym         % discretized input matrix (symbolic)
        Ad_nom         % discretized dynamics matrix evaluated at nominal parameter values
        Bd_nom         % discretized input matrix evaluated at nominal parameter values
        C_nom          % output matrix evaluated at nominal parameter values 
        D_nom          % feedthrough matrix evaluated at nominal parameter values 
        u_fun          % matlab function handle for computing numeric values of input
        x0_nom         % numeric value of initial condition 
        
        sensitivities_computed = false % flag to determine whether sensitivities were computed
        sensitivity_Ad % sensitivity of Ad with respect to model parameters
        sensitivity_Bd % sensitivity of Ad with respect to model parameters
        sensitivity_C  % sensitivity of C with respect to model parameters
        sensitivity_D  % sensitivity of D with respect to model parameters 
        sensitivity_u  % sensitivity of u with respect to model parameters
        sensitivity_x0 % sensitivity of x0 with respect to model parameters
    end
    
    methods
        function model = set.A(model, A_val)  % set method for A       
            if size(A_val, 1) ~= size(A_val, 2)
                error('flipAngleDesignTools:systemSizeMismatch', 'A matrix must be square')
            else
                model.A = A_val;
                model.n = size(A_val, 1);
            end
        end 
        function model = set.B(model, B_val)  % set method for B       
            if model.n ~= size(B_val, 1)
                error('flipAngleDesignTools:systemSizeMismatch', 'B matrix must have number of rows equal to the size of A')
            else
                model.B = B_val;
                model.ni = size(B_val, 2); 
            end
        end
        function model = set.C(model, C_val)  % set method for C       
            if model.n ~= size(C_val, 2)
                error('flipAngleDesignTools:systemSizeMismatch', 'C matrix must have number of columns equal to the size of A')
            else
                model.C = C_val;
                model.no = size(C_val, 1); 
            end
        end
        function model = set.D(model, D_val)  % set method for D       
            if model.no ~= size(D_val, 1)
                error('flipAngleDesignTools:systemSizeMismatch', 'D matrix must have number of rows equal to the number of rows of C')
            elseif model.ni ~= size(D_val, 2)
                error('flipAngleDesignTools:systemSizeMismatch', 'D matrix must have number of columns equal to the number of columns of B')
            else
                model.D = D_val;
            end
        end
        function model = set.flip_angle_input_matrix(model, M) % set method for flip_angle_input_matrix
            if size(M, 1) ~= model.n 
                error('flipAngleDesignTools:systemSizeMismatch', 'linear_exchange_model.flip_angle_input_matrix must have a number of rows equal to the number of inputs plus the number of outputs') 
            else
                model.flip_angle_input_matrix = M; 
                model.nq = size(M, 2); 
            end
        end
        function model = set.noise_parameters(model, values)    % set method for noise parameter values
            if ~isa(values, 'numeric')
                error('flipAngleDesignTools:systemSizeMismatch', 'linear_exchange_model.noise_parameters must have numeric values (unknown symbolic values not supported)') 
            elseif size(values, 1) ~= 1
                error('flipAngleDesignTools:systemSizeMismatch', 'linear_exchange_model.noise_parameters must be a row vector') 
            elseif size(values, 2) ~= model.no
                error('flipAngleDesignTools:systemSizeMismatch', 'linear_exchange_model.noise_parameters must be a row vector with number of entries equal to number of rows of C') 
            else
                model.noise_parameters = values; 
            end 
        end
        function model = set.x0(model, x0_val)% set method for x0      
            if size(x0_val, 1) ~= model.n 
                error('flipAngleDesignTools:systemSizeMismatch', 'x0 must be a column vector with number of rows equal to size of A') 
            elseif size(x0_val, 2) ~= 1
                error('flipAngleDesignTools:systemSizeMismatch', 'x0 must be a column vector with number of rows equal to size of A') 
            else
                model.x0 = x0_val; 
            end
            
        end            
        function model = discretize(model)    % compute discretization 
            
            % compute symbolic system discretization  
            model.Ad_sym = expm(model.TR*model.A);
            model.Bd_sym = model.A\((model.Ad_sym - eye(size(model.A)))*model.B);

            % evaluate system matrices, input and initial state at nominal parameter values 
            model.Ad_nom = double(subs(model.Ad_sym, [model.parameters_of_interest, model.nuisance_parameters, model.known_parameters], ...
                [model.parameters_of_interest_nominal_values, model.nuisance_parameters_nominal_values, model.known_parameter_values]));
            model.Bd_nom = double(subs(model.Bd_sym, [model.parameters_of_interest, model.nuisance_parameters, model.known_parameters], ...
                [model.parameters_of_interest_nominal_values, model.nuisance_parameters_nominal_values, model.known_parameter_values]));
            model.C_nom = double(subs(model.C, [model.parameters_of_interest, model.nuisance_parameters, model.known_parameters], ...
                [model.parameters_of_interest_nominal_values, model.nuisance_parameters_nominal_values, model.known_parameter_values]));
            model.D_nom = double(subs(model.D, [model.parameters_of_interest, model.nuisance_parameters, model.known_parameters], ...
                [model.parameters_of_interest_nominal_values, model.nuisance_parameters_nominal_values, model.known_parameter_values]));
            % model.u_fun = matlabFunction(subs(model.u,[model.parameters_of_interest, model.nuisance_parameters, model.known_parameters], ...
            %     [model.parameters_of_interest_nominal_values, model.nuisance_parameters_nominal_values, model.known_parameter_values]));
            model.x0_nom = double(subs(model.x0, [model.parameters_of_interest, model.nuisance_parameters, model.known_parameters], ...
                [model.parameters_of_interest_nominal_values, model.nuisance_parameters_nominal_values, model.known_parameter_values]));
                        
            % set flag 
            model.discretized = true; 
            
        end
        function model = sensitivities(model) % compute sentitivites   
            
            % compute discretized model (if necessary) 
            if ~model.discretized 
                model = discretize(model); 
            end
            
            display('===== Computing model sensitivities =====')
            
            % define vectors of parameters 
            parameters = [model.parameters_of_interest, model.nuisance_parameters];
            parameters_all = [parameters, model.known_parameters]; 
            parameters_all_vals = [model.parameters_of_interest_nominal_values, model.nuisance_parameters_nominal_values, model.known_parameter_values]; 

            % differentiate Ad with respect to parameters 
            model.sensitivity_Ad = zeros(size(model.Ad_sym, 1), size(model.Ad_sym, 2), length(parameters)); 
            for i=1:length(parameters)
                model.sensitivity_Ad(:, :, i) = double(subs(diff(model.Ad_sym, parameters(i)), parameters_all, parameters_all_vals)); 
            end
            
            % differentiate Bd with respect to parameters 
            model.sensitivity_Bd = zeros(size(model.Bd_sym, 1), size(model.Bd_sym, 2), length(parameters)); 
            for i=1:length(parameters)
                model.sensitivity_Bd(:, :, i) = double(subs(diff(model.Bd_sym, parameters(i)), parameters_all, parameters_all_vals)); 
            end
            
            % differentiate C with respect to parameters 
            model.sensitivity_C = zeros(size(model.C, 1), size(model.C, 2), length(parameters)); 
            for i=1:length(parameters)
                model.sensitivity_C(:, :, i) = double(subs(diff(model.C, parameters(i)), parameters_all, parameters_all_vals)); 
            end
            
            % differentiate D with respect to parameters 
            model.sensitivity_D = zeros(size(model.D, 1), size(model.D, 2), length(parameters)); 
            for i=1:length(parameters)
                model.sensitivity_D(:, :, i) = double(subs(diff(model.D, parameters(i)), parameters_all, parameters_all_vals)); 
            end
            
            % differentiate u with respect to parameters 
            model.sensitivity_u = zeros(model.ni, model.N, length(parameters)); 
            % syms t
   
            for i = 1:length(parameters)
                model.sensitivity_u(:, :, i) = double(subs(diff(model.u, parameters(i)), parameters_all, parameters_all_vals)); 
            end

            
            % differentiate x0 with respect to parameters 
            model.sensitivity_x0 = zeros(model.n, length(parameters)); 
            for i=1:length(parameters)
                model.sensitivity_x0(:, i) = double(subs(diff(model.x0, parameters(i)), parameters_all, parameters_all_vals)); 
            end
            
            % set flag 
            model.sensitivities_computed = true; 
            
            display('done.')
            
        end
    end
    
end


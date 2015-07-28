function y_mean = trajectories_2param(parameters, known_parameters, Pi, alpha, u)

% number of discretization points within profile
M = length(Pi) - 1;

%% Specify system model

% parameters
kTRANS = parameters(1);
kPL    = parameters(2);
R1P    = known_parameters(1);
R1L    = known_parameters(2);
B1     = 1; 

% compute B1-modified flip angles 
alpha = B1*alpha; 


% define system matrices for differential eq.
%   dx/dt = A*x(t) + B*u(t)
%    y(t) = C*x(t) + D*u(t)

% two-site exchange model with input feedthrough
A = [ -kPL-R1P  0  ;
    kPL     -R1L];

B = [kTRANS; 0];

C = eye(2);

D = [0; 0];

% define repetition time
TR = 2;

% define number of acquisitions
N = 30;

% discretize system
sys = ss(A, B, C, D);
sysd = c2d(sys, TR);
Ad = sysd.A;
Bd = sysd.B;

%% Compute trajectories

x = zeros(2, N, M+1);
y = zeros(2, N-1, M+1);
for i=1:M+1
    for t=1:N
        x(:, t+1, i) = Ad * cos(alpha(:, t).*(Pi(i))) .* x(:, t, i) + Bd*u(t);
        y(:, t, i)   = sin(alpha(:, t).*(Pi(i))) .* x(:, t, i);
    end
end

pyr_output_surface = reshape(y(1, :, :), N, M+1);
la_output_surface  = reshape(y(2, :, :), N, M+1);

%% Compare z-averaged signal with "ideal" central signal

y_mean  = [sum(pyr_output_surface, 2) sum(la_output_surface, 2)];

end


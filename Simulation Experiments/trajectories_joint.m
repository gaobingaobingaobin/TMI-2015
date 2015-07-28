function y = trajectories_joint(parameters, known_parameters, alpha, u)


%% Specify system model

% parameters
kTRANS = known_parameters(1);
kPL    = parameters(1);
R1P    = parameters(2);
R1L    = parameters(3);

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

x = zeros(2, N);
y = zeros(2, N-1);

for t=1:N
    x(:, t+1) = Ad * diag(cos(alpha(:, t))) * x(:, t) + Bd*u(t);
    y(:, t)   = diag(sin(alpha(:, t))) * x(:, t);
end

end


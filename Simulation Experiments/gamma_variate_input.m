function u  = gamma_variate_input( params, thetas )

alpha_1_in = params(1);
beta_1_in = params(2); 
A0_in = params(3); 
t0_in = params(4);

t = 2*(1:25);
u_true = A0_in* 1000* (max( t - t0_in, zeros(size(t)))).^alpha_1_in .*exp(-(t - t0_in)/beta_1_in);
u = sin(thetas) .* u_true';

end


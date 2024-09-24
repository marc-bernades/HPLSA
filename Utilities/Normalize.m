function [u_star_0, T_star_0, rho_star_0, P_star_0, mu_star_0, kappa_star_0, c_p_star_0, c_v_star_0, E_star_0, e_star_0, lambda_star_0, y_star] = Normalize(norm, u_0, T_0, rho_0, P_0, mu_0, kappa_0, c_p_0, c_v_0, E_0, e_0, lambda_0, y)

u_star_0      = u_0/norm.u;
T_star_0      = T_0/norm.T;
rho_star_0    = rho_0/norm.rho;
P_star_0      = P_0/norm.P;
mu_star_0     = mu_0/norm.mu;
kappa_star_0  = kappa_0/norm.kappa;
c_p_star_0    = c_p_0/norm.c_p;
c_v_star_0    = c_v_0/norm.c_v;
E_star_0      = E_0/norm.E;
e_star_0      = e_0/norm.e;
lambda_star_0 = lambda_0/norm.lambda;
y_star        = y/norm.y;

end
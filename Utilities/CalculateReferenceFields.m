function [rho_b,u_b,T_b,mu_b,kappa_b,c_p_b,c_v_b,E_b,e_b,sos_b] = CalculateReferenceFields(y,rho,u,T,mu,kappa,c_p,c_v,E,e,sos)


% Reference index at velocity max
[~,idx] = max(abs(u));

% Reference quantities
rho_b   = rho(idx);
u_b     = u(idx);
T_b     = T(idx);
mu_b    = mu(idx);
kappa_b = kappa(idx);
c_p_b   = c_p(idx);
c_v_b   = c_v(idx);
E_b     = E(idx);
e_b     = e(idx);
sos_b   = sos(idx);


end


function [rho_b,u_b,T_b,mu_b,kappa_b,c_p_b,c_v_b,E_b,e_b,sos_b] = CalculateBulkFields(y,rho,u,T,mu,kappa,c_p,c_v,E,e,sos)

rho_b   = 0;
u_b     = 0;
T_b     = 0;
mu_b    = 0;
kappa_b = 0;
c_p_b   = 0;
c_v_b   = 0;
E_b     = 0;
e_b     = 0;
sos_b   = 0;

y_total = abs(y(end) - y(1));

for jj = 1:length(y)

    if jj == 1
        dy = y(jj+1) - y(jj);
    elseif jj == length(y)
        dy = y(jj) - y(jj-1);
    else
        dy = abs(0.5*(y(jj+1) - y(jj-1)));
    end
    rho_b   = rho_b   + rho(jj).*dy/y_total;
    u_b     = u_b     + u(jj).*dy/y_total;
    T_b     = T_b     + T(jj).*dy/y_total;
    mu_b    = mu_b    + mu(jj).*dy/y_total;
    kappa_b = kappa_b + kappa(jj).*dy/y_total;
    c_p_b   = c_p_b   + c_p(jj).*dy/y_total;
    c_v_b   = c_v_b   + c_v(jj).*dy/y_total;
    E_b     = E_b     + E(jj).*dy/y_total;
    e_b     = e_b     + e(jj).*dy/y_total;
    sos_b   = sos_b   + sos(jj).*dy/y_total;

end


% % Standard version
% y_total = abs(y(end-1) - y(2));
% 
% for jj = 2:length(y)-1
% 
%     dy = abs(0.5*(y(jj+1) - y(jj-1)));
%     rho_b   = rho_b   + rho(jj).*dy/y_total;
%     u_b     = u_b     + u(jj).*dy/y_total;
%     T_b     = T_b     + T(jj).*dy/y_total;
%     mu_b    = mu_b    + mu(jj).*dy/y_total;
%     kappa_b = kappa_b + kappa(jj).*dy/y_total;
%     c_p_b   = c_p_b   + c_p(jj).*dy/y_total;
%     c_v_b   = c_v_b   + c_v(jj).*dy/y_total;
%     E_b     = E_b     + E(jj).*dy/y_total;
%     e_b     = e_b     + e(jj).*dy/y_total;
%     sos_b   = sos_b   + sos(jj).*dy/y_total;
% 
% end



end
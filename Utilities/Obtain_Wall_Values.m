function [u_tau_bw, u_tau_tw, Re_tau_bw, Re_tau_tw, Re_tau_bw_max, Re_tau_tw_max, ...
    Cf_bw, Cf_tw, Nu_bw, Nu_tw, Pr_tau_bw, Pr_tau_tw, Ec_tau_bw, Ec_tau_tw, Ma_tau_bw, Ma_tau_tw, ...
    T_tau_bw, T_tau_tw] = Obtain_Wall_Values(x,y,z,avg_u,avg_rho,avg_T,avg_mu,avg_kappa,avg_c_p,avg_sos,u_b,rho_b,T_b,u_max, delta)

total_area_bw = 0;
mu_bw = 0;
rho_bw = 0;
T_bw = 0;
kappa_bw = 0;
c_p_bw = 0;
sos_bw = 0;

u_inner_bw = 0;
u_boundary_bw = 0;
T_inner_bw = 0;
T_boundary_bw = 0;

total_area_tw = 0;
mu_tw = 0;
rho_tw = 0;
T_tw = 0;
kappa_tw = 0;
c_p_tw = 0;
sos_tw = 0;
u_inner_tw = 0;
u_boundary_tw = 0;
T_inner_tw = 0;
T_boundary_tw = 0;


num_points_x = length(x(:,1,1));
num_points_y = length(y(1,:,1));
num_points_z = length(z(1,1,:));


for ii = 2:num_points_x-1
    for kk = 2:num_points_z-1

        % Bottom wall
        jj = 1;
        delta_x = 0.5*(x(ii+1,jj,kk) - x(ii-1,jj,kk));
        delta_z = 0.5*(z(ii,jj,kk+1) - z(ii,jj,kk-1));
        area_bw = delta_x*delta_z;
        total_area_bw = total_area_bw + area_bw;

        mu_bw    = mu_bw    + area_bw*0.5*(avg_mu(ii,jj,kk) + avg_mu(ii,jj+1,kk));
        rho_bw   = rho_bw   + area_bw*0.5*(avg_rho(ii,jj,kk) + avg_rho(ii,jj+1,kk));
        T_bw     = T_bw     + area_bw*0.5*(avg_T(ii,jj,kk) + avg_T(ii,jj+1,kk));
        kappa_bw = kappa_bw + area_bw*0.5*(avg_kappa(ii,jj,kk) + avg_kappa(ii,jj+1,kk));
        c_p_bw   = c_p_bw   + area_bw*0.5*(avg_c_p(ii,jj,kk) + avg_c_p(ii,jj+1,kk));
        sos_bw   = sos_bw   + area_bw*0.5*(avg_sos(ii,jj,kk) + avg_sos(ii,jj+1,kk));

        u_inner_bw    = u_inner_bw + area_bw*(avg_u(ii,jj+1,kk));
        u_boundary_bw = u_boundary_bw + area_bw*(avg_u(ii,jj,kk));
        T_inner_bw    = T_inner_bw + area_bw*(avg_T(ii,jj+1,kk));
        T_boundary_bw = T_boundary_bw + area_bw*(avg_T(ii,jj,kk));


        % Top wall
        jj = num_points_y;
        delta_x = 0.5*(x(ii+1,jj,kk) - x(ii-1,jj,kk));
        delta_z = 0.5*(z(ii,jj,kk+1) - z(ii,jj,kk-1));
        area_tw = delta_x*delta_z;
        total_area_tw = total_area_tw + area_tw;

        mu_tw    = mu_tw    + area_tw*0.5*(avg_mu(ii,jj,kk) + avg_mu(ii,jj-1,kk));
        rho_tw   = rho_tw   + area_tw*0.5*(avg_rho(ii,jj,kk) + avg_rho(ii,jj-1,kk));
        T_tw     = T_tw     + area_tw*0.5*(avg_T(ii,jj,kk) + avg_T(ii,jj-1,kk));
        kappa_tw = kappa_tw + area_tw*0.5*(avg_kappa(ii,jj,kk) + avg_kappa(ii,jj-1,kk));
        c_p_tw   = c_p_tw   + area_tw*0.5*(avg_c_p(ii,jj,kk) + avg_c_p(ii,jj-1,kk));
        sos_tw   = sos_tw   + area_bw*0.5*(avg_sos(ii,jj,kk) + avg_sos(ii,jj-1,kk));

        u_inner_tw    = u_inner_tw    + area_tw*(avg_u(ii,jj-1,kk));
        u_boundary_tw = u_boundary_tw + area_tw*(avg_u(ii,jj,kk));
        T_inner_tw    = T_inner_tw    + area_tw*(avg_T(ii,jj-1,kk));
        T_boundary_tw = T_boundary_tw + area_tw*(avg_T(ii,jj,kk));

    end


end


% Normalize weighting to total area
mu_bw        = mu_bw/total_area_bw;
rho_bw        = rho_bw/total_area_bw;
T_bw          = T_bw/total_area_bw;
kappa_bw      = kappa_bw/total_area_bw;
c_p_bw        = c_p_bw/total_area_bw;
u_inner_bw    = u_inner_bw/total_area_bw;
u_boundary_bw = u_boundary_bw/total_area_bw;
T_inner_bw    = T_inner_bw/total_area_bw;
T_boundary_bw = T_boundary_bw/total_area_bw;
sos_bw        = sos_bw/total_area_bw;

mu_tw        = mu_tw/total_area_tw;
rho_tw        = rho_tw/total_area_tw;
T_tw          = T_tw/total_area_tw;
kappa_tw      = kappa_tw/total_area_tw;
c_p_tw        = c_p_tw/total_area_tw;
u_inner_tw    = u_inner_tw/total_area_tw;
u_boundary_tw = u_boundary_tw/total_area_tw;
T_inner_tw    = T_inner_tw/total_area_tw;
T_boundary_tw = T_boundary_tw/total_area_tw;
sos_tw        = sos_tw/total_area_tw;


% delta y
delta_y_bw = 0.5*(y(1,2,1)   - y(1,1,1));
delta_y_tw = 0.5*(y(1,end,1) - y(1,end-1,1));

% Tau wall
tau_bw = mu_bw*(u_inner_bw - u_boundary_bw)/delta_y_bw;
tau_tw = mu_tw*(u_inner_tw - u_boundary_tw)/delta_y_tw;

% u tau
u_tau_bw = sqrt(tau_bw/rho_bw);
u_tau_tw = sqrt(tau_tw/rho_tw);

% Reynolds tau
Re_tau_bw = u_tau_bw*delta*rho_bw/mu_bw;
Re_tau_tw = u_tau_tw*delta*rho_tw/mu_tw;

% Reynolds tau umax
Re_tau_bw_max = u_max*delta*rho_bw/mu_bw;
Re_tau_tw_max = u_max*delta*rho_tw/mu_tw;

% Skin friction
Cf_bw = tau_bw/(0.5*rho_b*u_b^2);
Cf_tw = tau_tw/(0.5*rho_b*u_b^2);

% Nusselt number tau
Nu_bw = - delta*((T_inner_bw - T_boundary_bw)/delta_y_bw)/(T_bw - T_b);
Nu_tw =   delta*((T_boundary_tw - T_inner_tw)/delta_y_tw)/(T_tw - T_b);


% T tau wall
T_tau_bw =   kappa_bw*(T_inner_bw - T_boundary_bw)/delta_y_bw/(rho_bw*c_p_bw*u_tau_bw);
T_tau_tw = - kappa_tw*(T_inner_tw - T_boundary_tw)/delta_y_tw/(rho_tw*c_p_tw*u_tau_tw);


% Prandtl wall
Pr_tau_bw = mu_bw.*c_p_bw/kappa_bw;
Pr_tau_tw = mu_tw.*c_p_tw/kappa_tw;

% Eckert Wall at u_max
Ec_tau_bw = u_tau_bw^2./(c_p_bw*T_bw);
Ec_tau_tw = u_tau_tw^2./(c_p_tw*T_tw);

% Ma at umax
Ma_tau_bw = u_tau_bw/sos_bw;
Ma_tau_tw = u_tau_tw/sos_tw;



end





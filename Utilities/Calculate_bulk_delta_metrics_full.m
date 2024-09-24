function [rho_b, mu_b, P_b, T_b, u_b, c_p_b, kappa_b, sos_b, volume, delta] = Calculate_bulk_delta_metrics_full(x,y,z,avg_rho, avg_mu, avg_P, avg_T, avg_u, avg_c_p, avg_kappa, avg_sos, num_points_x, num_points_y, num_points_z, delta_h, dx, dy, dz)

delta = zeros(size(x)); 
volume_total = 0;
rho_b = 0;
mu_b  = 0;
P_b   = 0;
T_b   = 0;
u_b   = 0;
c_p_b     = 0;
kappa_b   = 0;
sos_b     = 0;

for ii = 2:num_points_x-1
    for jj = 2:num_points_y-1
        for kk = 2:num_points_z-1


            volume          = dx(ii,jj,kk)*dy(ii,jj,kk)*dz(ii,jj,kk);

            delta(ii,jj,kk) = volume.^(1/3);

            % Bulk values
            volume_total    = volume_total + volume;
            rho_b = rho_b + avg_rho(ii,jj,kk)*volume;
            mu_b  = mu_b  + avg_mu(ii,jj,kk)*volume;
            P_b   = P_b   + avg_P(ii,jj,kk)*volume;
            T_b   = T_b   + avg_T(ii,jj,kk)*volume;
            u_b   = u_b   + avg_u(ii,jj,kk)*volume;
            c_p_b     = c_p_b   + avg_c_p(ii,jj,kk)*volume;
            kappa_b   = kappa_b + avg_kappa(ii,jj,kk)*volume;
            sos_b     = sos_b   + avg_sos(ii,jj,kk)*volume;

        end
    end
end

% Normalize bulk values
rho_b = rho_b/volume_total;
mu_b  = mu_b/volume_total;
P_b   = P_b/volume_total;
T_b   = T_b/volume_total;
u_b   = u_b/volume_total;
c_p_b   = c_p_b/volume_total;
kappa_b = kappa_b/volume_total;
sos_b   = sos_b/volume_total;


end
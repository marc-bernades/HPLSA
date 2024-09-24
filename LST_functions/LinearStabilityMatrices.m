function [Lt,Lx,Ly,Lz,Lq,Vxx,Vxy,Vxz,Vyy,Vyz,Vzz] = LinearStabilityMatrices(idx,u_0,rho_0,T_0,P_0,e_0,mu_0,kappa_0,lambda_0,Re,Pr,Ec,J_t,J_t2,Dy,q_vars,F_hat_0)

%% Eigen-problem initialization
% Build Lagrangian Matrices
nx = q_vars; ny = q_vars;
Lt  = zeros(nx,ny);
Lx  = zeros(nx,ny);
Ly  = zeros(nx,ny);
Lz  = zeros(nx,ny);
Lq  = zeros(nx,ny);

Vxx = zeros(nx,ny);
Vxy = zeros(nx,ny);
Vxz = zeros(nx,ny);
Vyy = zeros(nx,ny);
Vyz = zeros(nx,ny);
Vzz = zeros(nx,ny);


%% Assess matrices > Input fields are base flow
% Ref Ren (Pecnick) JFM 2019
Lt(1,1) = 1;
Lt(2,1) = u_0;
Lt(2,2) = rho_0;
Lt(3,3) = Lt(2,2);
Lt(4,4) = Lt(2,2);
Lt(5,1) = e_0 + rho_0*J_t.de_drho(idx);
Lt(5,5) = rho_0*J_t.de_dT(idx);

Lx(1,1) = u_0;
Lx(1,2) = rho_0;
Lx(2,1) = u_0*u_0 + J_t.dP_drho(idx);
Lx(2,2) = 2*rho_0*u_0;
Lx(2,3) = -1/Re*Dy.dmu_dy(idx);
Lx(2,5) = J_t.dP_dT(idx);
Lx(3,1) = -1/Re*J_t2.dmu_drho(idx)*Dy.du_dy(idx);
Lx(3,2) = -1/Re*Dy.dlambda_dy(idx);
Lx(3,3) = rho_0*u_0;
Lx(4,4) = Lx(3,3);
Lx(3,5) = -1/Re*J_t2.dmu_dT(idx)*Dy.du_dy(idx);
Lx(5,1) = e_0*u_0 + rho_0*u_0*J_t.de_drho(idx);
Lx(5,2) = rho_0*e_0 + P_0;
Lx(5,3) = -2/Re*mu_0*Dy.du_dy(idx);
Lx(5,5) = rho_0*u_0*J_t.de_dT(idx);

Ly(1,3) = rho_0;
Ly(2,1) = -1/Re*J_t2.dmu_drho(idx)*Dy.du_dy(idx);
Ly(2,2) = -1/Re*Dy.dmu_dy(idx);
Ly(2,3) = rho_0*u_0;
Ly(2,5) = -1/Re*J_t2.dmu_dT(idx)*Dy.du_dy(idx);
Ly(3,1) = J_t.dP_drho(idx);
Ly(3,3) = -2/Re*Dy.dmu_dy(idx) - 1/Re*Dy.dlambda_dy(idx);
Ly(3,5) = J_t.dP_dT(idx);
Ly(4,4) = -1/Re*Dy.dmu_dy(idx);
Ly(5,1) = -1/(Re*Pr*Ec)*(J_t2.dkappa_drho(idx)*Dy.dT_dy(idx));
Ly(5,2) = -2/Re*mu_0*Dy.du_dy(idx);
Ly(5,3) = rho_0*e_0 + P_0;
Ly(5,5) = -1/(Re*Pr*Ec)*(Dy.dkappa_dy(idx) + J_t2.dkappa_dT(idx)*Dy.dT_dy(idx));

Lz(1,4) = rho_0;
Lz(2,4) = rho_0*u_0;
Lz(3,4) = -1/Re*Dy.dlambda_dy(idx);
Lz(4,1) = J_t.dP_drho(idx);
Lz(4,3) = -1/Re*Dy.dmu_dy(idx);
Lz(4,5) = J_t.dP_dT(idx);
Lz(5,4) = rho_0*e_0 + P_0;

Lq(1,3) = Dy.drho_dy(idx);
Lq(2,1) = -1/Re*J_t2.dmu_drho(idx)*Dy.du2_dy2(idx) - 1/Re*Dy.du_dy(idx)*(J_t2.d2mu_d2rho(idx)*Dy.drho_dy(idx) + J_t2.d2mu_drhodT(idx)*Dy.dT_dy(idx));
Lq(2,3) = rho_0*Dy.du_dy(idx) + u_0*Dy.drho_dy(idx);
Lq(2,5) = -1/Re*J_t2.dmu_dT(idx)*Dy.du2_dy2(idx)   - 1/Re*Dy.du_dy(idx)*(J_t2.d2mu_d2T(idx)*Dy.dT_dy(idx) + J_t2.d2mu_dTdrho(idx)*Dy.drho_dy(idx));
Lq(3,1) = J_t.d2P_d2rho(idx)*Dy.drho_dy(idx) + J_t.d2P_drhodT(idx)*Dy.dT_dy(idx);
Lq(3,5) = J_t.d2P_d2T(idx)*Dy.dT_dy(idx)     + J_t.d2P_drhodT(idx)*Dy.drho_dy(idx);
Lq(5,1) = -1/(Re*Pr*Ec)*(Dy.dT2_dy2(idx)*J_t2.dkappa_drho(idx) + (J_t2.d2kappa_d2rho(idx)*Dy.drho_dy(idx) + J_t2.d2kappa_drhodT(idx)*Dy.dT_dy(idx))*Dy.dT_dy(idx)) - 1/Re*J_t2.dmu_drho(idx)*(Dy.du_dy(idx)^2);
Lq(5,2) = 1.0*-F_hat_0/Re;
Lq(5,3) = e_0*Dy.drho_dy(idx) + rho_0*Dy.de_dy(idx);
Lq(5,5) = -1/(Re*Pr*Ec)*(Dy.dT2_dy2(idx)*J_t2.dkappa_dT(idx) + (J_t2.d2kappa_d2T(idx)*Dy.dT_dy(idx) + J_t2.d2kappa_drhodT(idx)*Dy.drho_dy(idx))*Dy.dT_dy(idx)) - 1/Re*J_t2.dmu_dT(idx)*(Dy.du_dy(idx)^2);

Vxx(2,2) = -(2*mu_0 + lambda_0)/Re;
Vyy(3,3) = Vxx(2,2);
Vzz(4,4) = Vxx(2,2);

Vxx(3,3) = -mu_0/Re;
Vxx(4,4) = Vxx(3,3);
Vyy(2,2) = Vxx(3,3);
Vyy(4,4) = Vxx(3,3);
Vzz(2,2) = Vxx(3,3);
Vzz(3,3) = Vxx(3,3);

Vxx(5,5) = - kappa_0/(Re*Pr*Ec);
Vyy(5,5) = Vxx(5,5);
Vzz(5,5) = Vxx(5,5);
Vxy(2,3) = - (mu_0 + lambda_0)/Re;
Vxy(3,2) = Vxy(2,3);
Vxz(2,4) = Vxy(2,3);
Vxz(4,2) = Vxy(2,3);
Vyz(3,4) = Vxy(2,3);
Vyz(4,3) = Vxy(2,3);
end
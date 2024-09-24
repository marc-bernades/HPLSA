function [A,B] = LST(delta,N,BF,Dy,omega,beta,alpha,bSolver,HP_model, Substance, Fluid, bScaling, varargin)

% Imaginary
img     = sqrt(-1);

% Unpack structure and load variables
struct_field_names = fieldnames(BF);
for k = 1:length(struct_field_names)
      eval([struct_field_names{k}, ' = BF.',struct_field_names{k}, ';']);
end

% Ovewrite Re for LST Re Sweep
if ~isempty(varargin)
    Re = varargin{1};
end

% Update normalized quantities
[~,D] = Chevysheb_collocation(N,delta/BF.norm.y);

% Perturbation vector
q_vars = 5;                 % Number state variables
N      = N + 1;             % Chevyshev collocation grid points

% Re-cast Chevysheb discretization matrix
D_tot = zeros(N*q_vars,N*q_vars);
for ii = 1:length(D(:,1))
    vec_tot = [];
    for jj = 1:length(D(1,:))
        % Build the D ii-line layout given current q_vars
        vec_temp = zeros(1,q_vars);
        vec_temp(1) = D(ii,jj);
        vec_tot = [vec_tot vec_temp]; 
    end
    % Allocate pattern for all q_vars for D ii-line
    for n_var = 1:q_vars
        D_tot(ii*q_vars - (q_vars - 1) + n_var - 1,:) = circshift(vec_tot,n_var-1);
    end
end

% Numeric Jacobian
J_t2 = Calculate_Numeric_Thermodynamic_Jacobian(bSolver,c_v,T,rho, Fluid, Substance, HP_model);

% Thermodynamic Jacobians
J_t = Jacobian_thermodynamics(bSolver,rho,T,P_0,J_t2,Substance, HP_model);

% Normalize reference scalings
norm = Reference_Scalings(u_b,T_b,rho_b,mu_b,kappa_b,c_p_b,c_v_b,E_b,e_b,P_b(1),u,T,rho,mu,kappa,c_p,c_v,lambda,Dy,delta,bScaling, BF.norm.u);
[u, T, rho, P_0, mu, kappa, c_p, c_v, E, e, lambda, y] = Normalize(norm, u, T, rho, P_0, mu, kappa, c_p, c_v, E, e, lambda, BF.y);

% Normalize Dy
Dy.du_dy      = Dy.du_dy*norm.y/norm.u;      
Dy.dT_dy      = Dy.dT_dy*norm.y/norm.T;     
Dy.dmu_dy     = Dy.dmu_dy*norm.y/norm.mu;    
Dy.dkappa_dy  = Dy.dkappa_dy*norm.y/norm.kappa;  
Dy.drho_dy    = Dy.drho_dy*norm.y/norm.rho;    
Dy.dlambda_dy = Dy.dlambda_dy*norm.y/norm.lambda;
Dy.de_dy      = Dy.de_dy*norm.y/norm.e;      
Dy.du2_dy2    = Dy.du2_dy2*norm.y*norm.y/norm.u;    
Dy.dT2_dy2    = Dy.dT2_dy2*norm.y*norm.y/norm.T;   

% Normalize Jacobians
J_t2.dcv_dT         = J_t2.dcv_dT/norm.c_v*norm.T;
J_t2.dmu_dT         = J_t2.dmu_dT/norm.mu*norm.T;
J_t2.dkappa_dT      = J_t2.dkappa_dT/norm.kappa*norm.T;
J_t2.d2mu_d2T       = J_t2.d2mu_d2T/norm.mu*norm.T^2;
J_t2.d2kappa_d2T    = J_t2.d2kappa_d2T/norm.kappa*norm.T^2;
J_t2.dmu_drho       = J_t2.dmu_drho/norm.mu*norm.rho;
J_t2.d2mu_d2rho     = J_t2.d2mu_d2rho/norm.mu*norm.rho^2;
J_t2.dkappa_drho    = J_t2.dkappa_drho/norm.kappa*norm.rho;
J_t2.d2kappa_d2rho  = J_t2.d2kappa_d2rho/norm.kappa*norm.rho^2;
J_t2.d2mu_drhodT    = J_t2.d2mu_drhodT/norm.mu*norm.rho*norm.T;
J_t2.d2kappa_drhodT = J_t2.d2kappa_drhodT/norm.kappa*norm.rho*norm.T;

J_t.dP_dT       = J_t.dP_dT/norm.P*norm.T;
J_t.dP_drho     = J_t.dP_drho/norm.P*norm.rho;
J_t.d2P_drhodT  = J_t.d2P_drhodT/norm.P*norm.rho*norm.T;
J_t.d2P_d2T     = J_t.d2P_d2T/norm.P*norm.T^2;
J_t.d2P_d2rho   = J_t.d2P_d2rho/norm.P*norm.rho^2;
J_t.de_dT       = J_t.de_dT/norm.e*norm.T;
J_t.de_drho     = J_t.de_drho/norm.e*norm.rho;

% Initialize cellarray
Lt_vec  = cell(1,N);
Lx_vec  = cell(1,N);
Ly_vec  = cell(1,N);
Lz_vec  = cell(1,N);
Lq_vec  = cell(1,N);
Vxx_vec = cell(1,N);
Vxy_vec = cell(1,N);
Vxz_vec = cell(1,N);
Vyy_vec = cell(1,N);
Vyz_vec = cell(1,N);
Vzz_vec = cell(1,N);


% Build Jacobian Matrices to n_var x N
for ii = 1:N

    % For each data point we build the matrices
    [Lt_i,Lx_i,Ly_i,Lz_i,Lq_i,Vxx_i,Vxy_i,Vxz_i,Vyy_i,Vyz_i,Vzz_i] = LinearStabilityMatrices(ii,u(ii),rho(ii),T(ii),P_0(ii),e(ii),mu(ii),kappa(ii),lambda(ii),Re,Pr,Ec,J_t,J_t2,Dy,q_vars,F_hat_0);

    % Build full matrices (qvar x N) x (qvar x N)
    Lt_vec{ii} = Lt_i;   Lx_vec{ii} = Lx_i;   Ly_vec{ii} = Ly_i;   Lz_vec{ii} = Lz_i; Lq_vec{ii} = Lq_i;
    Vxx_vec{ii} = Vxx_i; Vxy_vec{ii} = Vxy_i; Vxz_vec{ii} = Vxz_i;
    Vyy_vec{ii} = Vyy_i; Vyz_vec{ii} = Vyz_i;
    Vzz_vec{ii} = Vzz_i;

end


% Total matrix - Set diagonal each grid point L_i matrix to obtain qxN x
% qxN matrices
Lt  = blkdiag(Lt_vec{:});
Lx  = blkdiag(Lx_vec{:});
Ly  = blkdiag(Ly_vec{:});
Lz  = blkdiag(Lz_vec{:});
Lq  = blkdiag(Lq_vec{:});
Vxx = blkdiag(Vxx_vec{:});
Vxy = blkdiag(Vxy_vec{:});
Vxz = blkdiag(Vxz_vec{:});
Vyy = blkdiag(Vyy_vec{:});
Vyz = blkdiag(Vyz_vec{:});
Vzz = blkdiag(Vzz_vec{:});


%% Eigen value problem
% A = alpha*Lx - img*Ly*D_tot + beta*Lz - img*Lq + img*alpha.^2*Vxx + alpha*Vxy*D_tot + img*alpha*beta*Vxz - img*Vyy*D_tot^2 + beta*Vyz*D_tot + img*beta^2*Vzz;
A = alpha*Lx - img*Ly*D_tot + beta*Lz - img*Lq + img*alpha.^2*Vxx + alpha*Vxy*D_tot + img*alpha*beta*Vxz - img*D_tot^2*Vyy + beta*D_tot*Vyz + img*beta^2*Vzz;


B = Lt; % Temporal B
% A = -img*omega*Lt + Ly*D_tot + img*beta*Lz + Lq + D_tot.^2*Vyy +
% img*beta*D_tot*Vyz - beta^2*Vzz; % Spatial A

%% Boundary conditions
% Fluctuations u' = v' = w' = T' = 0 at walls (not rho') > Fluctuations to 0 on the respective lines q = rho, u, v, w, T
A(2:5,:) = 0; A(end-3:end,:) = 0;
% Diagonal of these lines should have large value on A
A(2,2) = 1E4; A(3,3) = 1E4; A(4,4) = 1E4; A(5,5) = 1E4;
A(end-3,end-3) = 1E4; A(end-2,end-2) = 1E4; A(end-1,end-1) = 1E4; A(end,end) = 1E4;
% Idem as A
B(2:5,:) = 0; B(end-3:end,:) = 0;
B(2,2) = 1E4; B(3,3) = 1E4; B(4,4) = 1E4; B(5,5) = 1E4;
B(end-3,end-3) = 1E4; B(end-2,end-2) = 1E4; B(end-1,end-1) = 1E4; B(end,end) = 1E4;


end
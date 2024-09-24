%% Load results and calculations
clear all; clc; close all;
img = sqrt(-1);

bSolver           = 'Real';
HP_model          = 'CoolProp';  % CoolProp, RefProp, Constant, Power, LowPressure, HighPressure   

% Select fluid and initial conditions
Fluid.Substance         = 'N2';
[~, T_c, P_c, ~, rho_c, ~,  ~,  ~,  ~,  ~,  ~, ~,  ~,  ~,  ~,  ~,  ~,  ~] = Substance_Library(Fluid.Substance);

% File name

% % Validation CO2
% name_file_load =  'CO2_Tbw_290_Ttw_290_Pb_8000000_PrEc_wall_isothermal_Ke_Validation';

% % Isothermal flow cases (I1-3)
name_file_load =  'N2_Tbw_94.6425_Ttw_94.6425_Pb_5093700_PrEc_bulk_isothermal_Ke';
% name_file_load =  'N2_Tbw_119.8805_Ttw_119.8805_Pb_5093700_PrEc_bulk_isothermal_Ke';
% name_file_load =  'N2_Tbw_119.8805_Ttw_119.8805_Pb_5093700_PrEc_wall_isothermal_Ke';
% name_file_load =  'N2_Tbw_189.285_Ttw_189.285_Pb_5093700_PrEc_bulk_isothermal_Ke';

% % Non-isothermal flow cases (NI-1-6)
% name_file_load =  'N2_Tbw_94.6425_Ttw_189.285_Pb_5093700_PrEc_bulk_Non_isothermal_ke';
% name_file_load =  'N2_Tbw_94.6425_Ttw_189.285_Pb_16979000_PrEc_bulk_Non_isothermal_Ke';
% name_file_load =  'N2_Tbw_113.571_Ttw_138.809_Pb_5093700_PrEc_bulk_Non_isothermal_Ke';
% name_file_load =  'N2_Tbw_94.6425_Ttw_189.285_Pb_101874_PrEc_bulk_Non_isothermal_Ke';
% name_file_load =  'N2_Tbw_94.6425_Ttw_189.285_Pb_6791600_PrEc_bulk_Non_isothermal_Ke';
% name_file_load =  'N2_Tbw_94.6425_Ttw_189.285_Pb_101874_PrEc_bulk_Non_isothermal_Ke';


% Unpack structure and load variables
Results = load(strcat('Results/EnergyBudget/', name_file_load, '.mat'));
% Results = load(strcat('Results/Modal/', name_file_load, '.mat'));
struct_field_names = fieldnames(Results.Data);
for k = 1:length(struct_field_names)
      eval([struct_field_names{k}, ' = Results.Data.',struct_field_names{k}, ';']);
end

% Check isothermal case
bIsothermal = ~contains(name_file_load,'non_isothermal'); 

% Perturbation variables
q_vars = 5;
delta  = 1;


%% Dataset selection
beta_target  = 0;
alpha_target = 1.0;
Re_target = 10000;
Br_target = 0.05;
[~,aa] = min(abs(alpha-alpha_target));
[~,bb] = min(abs(beta-beta_target));
[~,ii] = min(abs(N_target-Br_target));
[~,jj] = min(abs(n_LST_Sweep-Re_target));

% Unpack baseflow
struct_field_names = fieldnames(BF{ii});
for k = 1:length(struct_field_names)
      eval([struct_field_names{k}, ' = BF{ii}.',struct_field_names{k}, ';']);
end

% Numeric Jacobian
J_t2 = Calculate_Numeric_Thermodynamic_Jacobian(bSolver,c_v,T,rho, Fluid, Fluid.Substance, HP_model);

% Thermodynamic Jacobians
J_t = Jacobian_thermodynamics(bSolver,rho,T,P_0,J_t2,Fluid.Substance, HP_model);


% Normalize reference scalings
norm = Reference_Scalings(u_b,T_b,rho_b,mu_b,kappa_b,c_p_b,c_v_b,E_b,e_b,P_b(1),u,T,rho,mu,kappa,c_p,c_v,lambda,Dy,delta,bScaling, BF{ii}.norm.u);
[u, T, rho, P_0, mu, kappa, c_p, c_v, E, e, lambda, y] = Normalize(norm, u, T, rho, P_0, mu, kappa, c_p, c_v, E, e, lambda, BF{1}.y);

% Normalize Dy
Dy = Dy{ii};
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

% Order eigenvalues imaginary part decreasing
[Val_eigen_dec, i_sort] = sort(-imag(Val_eigen{aa,ii,jj} ));
Val_eigen_temp          = Val_eigen{aa,ii,jj};
Vec_eigen_temp          = Vec_eigen{aa,ii,jj};
Val_eigen{aa,ii,jj}     = Val_eigen_temp(i_sort);
Vec_eigen{aa,ii,jj}     = Vec_eigen_temp(:,i_sort);

% Modes of interest
Imag_max = 0.1;
Imag_min = -1;
Real_max = 0.8;
Real_min = 0.0;

Val_eigen_tg = zeros(length(Val_eigen{aa,ii,jj}),1);
Vec_eigen_tg = zeros(size(Vec_eigen{aa,ii,jj}));

num     = 0;
num_INF = 0;

for i = 1:length(Val_eigen{aa,ii,jj}) %N*q_vars
    if(isinf(real(Val_eigen{aa,ii,jj}(i)))||isinf(imag(Val_eigen{aa,ii,jj}(i))))
        num_INF=num_INF+1;
    else
        if(imag(Val_eigen{aa,ii,jj}(i))>Imag_min...
                && imag(Val_eigen{aa,ii,jj}(i))<Imag_max...
                && real(Val_eigen{aa,ii,jj}(i))<Real_max...
                && real(Val_eigen{aa,ii,jj}(i))>Real_min)
            num = num+1;
            Val_eigen_tg(num)   = Val_eigen{aa,ii,jj}(i);
            Vec_eigen_tg(:,num) = Vec_eigen{aa,ii,jj}(:,i);
        end
    end
end


% Val_eigen_tg(num+1:N*q_vars)   = [];
% Vec_eigen_tg(:,num+1:N*q_vars) = [];
Val_eigen_tg(num+1:length(Val_eigen{aa,ii,jj}))   = [];
Vec_eigen_tg(:,num+1:length(Val_eigen{aa,ii,jj})) = [];

% Obtain perturbation vectors
for n_var = 1:q_vars
    Val_eigen_var{n_var} = Val_eigen_tg(n_var:q_vars:end);
    Vec_eigen_var{n_var} = Vec_eigen_tg(n_var:q_vars:end,:); % Take all columns for each perturbation row
end

% kinetic energy budget most unstable mode
Re = n_LST_Sweep(jj);


idx = 1; % Most unstable mode first % for idx = 1:length(Val_eigen_tg)
imagc       = imag(Val_eigen_tg);

[idx_val,~] = find(imagc ~= 0);
if imagc(idx) == 0
    idx = idx_val(1);
end
disp("Most unstable mode is " + num2str(Val_eigen_tg(idx)))

Norm = 1; %max(abs(real(Vec_eigen_var{2}(:,idx))));

q_rho = Vec_eigen_var{1}(:,idx)/Norm;
q_u   = Vec_eigen_var{2}(:,idx)/Norm;
q_v   = Vec_eigen_var{3}(:,idx)/Norm;
q_w   = Vec_eigen_var{4}(:,idx)/Norm;
q_T   = Vec_eigen_var{5}(:,idx)/Norm;

%     d_q_rho_dy = CentralDerivative_d1_2ndOrder_1D(q_rho,y)'; %Dy.D*real((q_rho));
%     d_q_u_dy   = CentralDerivative_d1_2ndOrder_1D(q_u,y)'; %Dy.D*(real(q_u));
%     d_q_v_dy   = CentralDerivative_d1_2ndOrder_1D(q_v,y)'; %Dy.D*(real(q_v));
%     d_q_T_dy   = CentralDerivative_d1_2ndOrder_1D(q_T,y)'; %Dy.D*(real(q_w));
%     d2_q_u_d2y = CentralDerivative_d2_2ndOrder_1D(q_u,y)'; %Dy.D*Dy.D*(real(q_u));
%     d2_q_v_d2y = CentralDerivative_d2_2ndOrder_1D(q_v,y)'; %Dy.D*Dy.D*(real(q_v));

d_q_rho_dy = Dy.D*((q_rho));
d_q_u_dy   = Dy.D*((q_u));
d_q_v_dy   = Dy.D*((q_v));
d_q_T_dy   = Dy.D*((q_T));
d2_q_u_d2y = Dy.D*Dy.D*((q_u));
d2_q_v_d2y = Dy.D*Dy.D*((q_v));

Ke_r{aa,ii,jj,idx} = -1.0*real(img*Val_eigen_tg(idx).*trapz(BF{ii}.y,BF{ii}.rho./BF{ii}.norm.rho.*(q_u.*conj(q_u) + q_v.*conj(q_v) + 0*q_w.*conj(q_w))));
Fi_r{aa,ii,jj,idx} = -1.0*real(img*alpha(aa).*trapz(BF{ii}.y,BF{ii}.rho./BF{ii}.norm.rho.*BF{ii}.u./BF{ii}.norm.u.*(q_u.*conj(q_u) + q_v.*conj(q_v) + 0*q_w.*conj(q_w))));
P_r{aa,ii,jj,idx}  = -1.0*real(trapz(BF{ii}.y,BF{ii}.rho./BF{ii}.norm.rho.*Dy.du_dy.*(q_v.*conj(q_u))));
T_r{aa,ii,jj,idx}  = -1.0*real(trapz(BF{ii}.y,img*alpha(aa)*J_t.dP_drho.*q_rho.*conj(q_u) + img*alpha(aa)*J_t.dP_dT.*q_T.*conj(q_u) + J_t.dP_drho.*d_q_rho_dy.*conj(q_v) + J_t.dP_dT.*d_q_T_dy.*conj(q_v) + ...
    (J_t.d2P_d2rho.*Dy.drho_dy + J_t.d2P_drhodT.*Dy.dT_dy).*q_rho.*conj(q_v) + ...
    (J_t.d2P_d2T.*Dy.dT_dy + J_t.d2P_drhodT.*Dy.drho_dy).*q_T.*conj(q_v)));
V_r{aa,ii,jj,idx}  = 1/Re.*real(trapz(BF{ii}.y,-alpha(aa)^2.*(2*mu_star + lambda).*q_u.*conj(q_u) + mu_star.*d2_q_u_d2y.*conj(q_u) + img*alpha(aa)*(mu_star + lambda).*d_q_v_dy.*conj(q_u) + ...
    img*alpha(aa)*Dy.dmu_dy.*q_v.*conj(q_u) + J_t2.dmu_drho.*Dy.du_dy.*d_q_rho_dy.*conj(q_u) + Dy.dmu_dy.*d_q_u_dy.*conj(q_u) + J_t2.dmu_dT.*Dy.du_dy.*d_q_T_dy.*conj(q_u) + ...
    J_t2.dmu_drho.*Dy.du2_dy2.*q_rho.*conj(q_u) + Dy.du_dy.*(J_t2.d2mu_d2rho.*Dy.drho_dy + J_t2.d2mu_drhodT.*Dy.dT_dy).*q_rho.*conj(q_u) + ...
    J_t2.dmu_dT.*Dy.du2_dy2.*q_T.*conj(q_u) + Dy.du_dy.*(J_t2.d2mu_d2T.*Dy.dT_dy + J_t2.d2mu_drhodT.*Dy.drho_dy).*q_T.*conj(q_u) + ...
    - alpha(aa)^2.*mu_star.*q_v.*conj(q_v) + (2*mu_star + lambda).*d2_q_v_d2y.*conj(q_v) + img*alpha(aa)*(mu_star + lambda).*d_q_u_dy.*conj(q_v) + ...
    img*alpha(aa).*J_t2.dmu_drho.*Dy.du_dy.*q_rho.*conj(q_v) + img*alpha(aa)*Dy.dlambda_dy.*q_u.*conj(q_v) + img*alpha(aa).*J_t2.dmu_dT.*Dy.du_dy.*q_T.*conj(q_v) + (2*Dy.dmu_dy + Dy.dlambda_dy).*d_q_v_dy.*conj(q_v)));

% figure
% plot(y,rho.*Dy.du_dy); hold on
% plot(y,real(q_v.*conj(q_u))*20)
% figure
% plot(y,real(q_v.*conj(q_u))); hold on
% figure
% plot(y,real(q_v.*conj(q_u))*norm.u^2); hold on

disp("Ke = " + num2str(Ke_r{aa,ii,jj,idx}*1E4))
disp("P = " + num2str(P_r{aa,ii,jj,idx}*1E4))
disp("T = " + num2str(T_r{aa,ii,jj,idx}*1E4))
disp("V = " + num2str(V_r{aa,ii,jj,idx}*1E4))
disp("Sum = " + num2str((P_r{aa,ii,jj,idx} + T_r{aa,ii,jj,idx} + V_r{aa,ii,jj,idx})*1E4))
 

% end


% Plot spatial distribution
Ke_r_plot{aa,ii,jj,idx} = -1.0*real(img*Val_eigen_tg(idx).*BF{ii}.rho./BF{ii}.norm.rho.*(q_u.*conj(q_u) + q_v.*conj(q_v) + 0*q_w.*conj(q_w)));
Fi_r_plot{aa,ii,jj,idx} = -1.0*real(img*alpha(aa).*BF{ii}.rho./BF{ii}.norm.rho.*BF{ii}.u./BF{ii}.norm.u.*(q_u.*conj(q_u) + q_v.*conj(q_v) + 0*q_w.*conj(q_w)));
P_r_plot{aa,ii,jj,idx}  = -1.0*real(BF{ii}.rho./BF{ii}.norm.rho.*Dy.du_dy.*(q_v.*conj(q_u)));
T_r_plot{aa,ii,jj,idx}  = -1.0*real(img*alpha(aa)*J_t.dP_drho.*q_rho.*conj(q_u) + img*alpha(aa)*J_t.dP_dT.*q_T.*conj(q_u) + J_t.dP_drho.*d_q_rho_dy.*conj(q_v) + J_t.dP_dT.*d_q_T_dy.*conj(q_v) + ...
    (J_t.d2P_d2rho.*Dy.drho_dy + J_t.d2P_drhodT.*Dy.dT_dy).*q_rho.*conj(q_v) + ...
    (J_t.d2P_d2T.*Dy.dT_dy + J_t.d2P_drhodT.*Dy.drho_dy).*q_T.*conj(q_v));
V_r_plot{aa,ii,jj,idx}  = 1/Re.*real(-alpha(aa)^2.*(2*mu_star + lambda).*q_u.*conj(q_u) + mu_star.*d2_q_u_d2y.*conj(q_u) + img*alpha(aa)*(mu_star + lambda).*d_q_v_dy.*conj(q_u) + ...
    img*alpha(aa)*Dy.dmu_dy.*q_v.*conj(q_u) + J_t2.dmu_drho.*Dy.du_dy.*d_q_rho_dy.*conj(q_u) + Dy.dmu_dy.*d_q_u_dy.*conj(q_u) + J_t2.dmu_dT.*Dy.du_dy.*d_q_T_dy.*conj(q_u) + ...
    J_t2.dmu_drho.*Dy.du2_dy2.*q_rho.*conj(q_u) + Dy.du_dy.*(J_t2.d2mu_d2rho.*Dy.drho_dy + J_t2.d2mu_drhodT.*Dy.dT_dy).*q_rho.*conj(q_u) + ...
    J_t2.dmu_dT.*Dy.du2_dy2.*q_T.*conj(q_u) + Dy.du_dy.*(J_t2.d2mu_d2T.*Dy.dT_dy + J_t2.d2mu_drhodT.*Dy.drho_dy).*q_T.*conj(q_u) + ...
    - alpha(aa)^2.*mu_star.*q_v.*conj(q_v) + (2*mu_star + lambda).*d2_q_v_d2y.*conj(q_v) + img*alpha(aa)*(mu_star + lambda).*d_q_u_dy.*conj(q_v) + ...
    img*alpha(aa).*J_t2.dmu_drho.*Dy.du_dy.*q_rho.*conj(q_v) + img*alpha(aa)*Dy.dlambda_dy.*q_u.*conj(q_v) + img*alpha(aa).*J_t2.dmu_dT.*Dy.du_dy.*q_T.*conj(q_v) + (2*Dy.dmu_dy + Dy.dlambda_dy).*d_q_v_dy.*conj(q_v));



%% Plot multiple setups
% I-2 Br = 0.1 vs. 0.5
save_label = 'I_2_Br_0_1_0_5';
name_file_load_plot{1} = 'N2_Tbw_119.8805_Ttw_119.8805_Pb_5093700_PrEc_bulk_isothermal_Ke';
name_file_load_plot{2} = 'N2_Tbw_119.8805_Ttw_119.8805_Pb_5093700_PrEc_bulk_isothermal_Ke';
[K_plot{1},P_plot{1},T_plot{1},V_plot{1},~,~,~,~,~,~,~,~,~,~,~] = CalculateEnergybudget(name_file_load_plot{1},0, 1, 10000, 0.1, bSolver,HP_model,Fluid.Substance);
[K_plot{2},P_plot{2},T_plot{2},V_plot{2},~,~,~,~,~,~,~,~,~,~,~] = CalculateEnergybudget(name_file_load_plot{2},0, 1, 10000, 0.5, bSolver,HP_model,Fluid.Substance);
test_label{1} = '$I-2 \thinspace Br = 0.1$';
test_label{2} = '$I-2 \thinspace Br = 0.5$';
plot_EnergyBudgetSpatial(y,delta,K_plot,P_plot,T_plot,V_plot,test_label,save_label)

% I-1 vs I-3 at Br = 0.5
save_label = 'I_1_I_3_Br_0_5';
name_file_load_plot{1} = 'N2_Tbw_94.6425_Ttw_94.6425_Pb_5093700_PrEc_bulk_isothermal_Ke';
name_file_load_plot{2} = 'N2_Tbw_189.285_Ttw_189.285_Pb_5093700_PrEc_bulk_isothermal_Ke';
[K_plot{1},P_plot{1},T_plot{1},V_plot{1},~,~,~,~,~,~,~,~,~,~,~] = CalculateEnergybudget(name_file_load_plot{1},0, 1, 10000, 0.5, bSolver,HP_model,Fluid.Substance);
[K_plot{2},P_plot{2},T_plot{2},V_plot{2},~,~,~,~,~,~,~,~,~,~,~] = CalculateEnergybudget(name_file_load_plot{2},0, 1, 10000, 0.5, bSolver,HP_model,Fluid.Substance);
test_label{1} = '$I-1 \thinspace Br = 0.5$';
test_label{2} = '$I-3 \thinspace Br = 0.5$';
plot_EnergyBudgetSpatial(y,delta,K_plot,P_plot,T_plot,V_plot,test_label,save_label)


%% Plot single Non-isothermal
% NI-1 at Br = 0.1
save_label = 'NI_1_Br_0_1';
name_file_load_plot{1} = 'N2_Tbw_94.6425_Ttw_189.285_Pb_5093700_PrEc_bulk_Non_isothermal_Ke';
[K_plot{1},P_plot{1},T_plot{1},V_plot{1},~,~,~,~,~,~,~,~,~] = CalculateEnergybudget(name_file_load_plot{1},0, 0.6, 10000, 0.1, bSolver,HP_model,Fluid.Substance);
plot_EnergyBudgetSpatial_Single(y,delta,K_plot,P_plot,T_plot,V_plot,save_label)

% NI-3 at Br = 0.1
save_label = 'NI_3_Br_0_1';
name_file_load_plot{1} = 'N2_Tbw_113.571_Ttw_138.809_Pb_5093700_PrEc_bulk_Non_isothermal_Ke';
[K_plot{1},P_plot{1},T_plot{1},V_plot{1},~,~,~,~,~,~,~,~,~,~,~] = CalculateEnergybudget(name_file_load_plot{1},0, 0.8, 10000, 0.1, bSolver,HP_model,Fluid.Substance);
plot_EnergyBudgetSpatial_Single(y,delta,K_plot,P_plot,T_plot,V_plot,save_label)


%% Instability mechanism
% I-2 Br = 0.1 vs. 0.5
save_label = 'I_2_Br_0_1_0_5';
name_file_load_plot{1} = 'N2_Tbw_119.8805_Ttw_119.8805_Pb_5093700_PrEc_bulk_isothermal_Ke';
name_file_load_plot{2} = 'N2_Tbw_119.8805_Ttw_119.8805_Pb_5093700_PrEc_bulk_isothermal_Ke';
[~,~,~,~,W1{1},W2{1},W3{1},P_V{1},~,~,~,~,~,~,~] = CalculateEnergybudget(name_file_load_plot{1},0, 1, 10000, 0.1, bSolver,HP_model,Fluid.Substance);
[~,~,~,~,W1{2},W2{2},W3{2},P_V{2},~,~,~,~,~,~,~] = CalculateEnergybudget(name_file_load_plot{2},0, 1, 10000, 0.5, bSolver,HP_model,Fluid.Substance);
test_label{1} = '$Br = 0.1$';
test_label{2} = '$Br = 0.5$';
plot_EnergyBudgetInstability(y,delta,W1,W2,W3,P_V,test_label,save_label)

% I-2 Br = 0.1 vs. 0.5
save_label = 'NI_1_5_Br_0_1';
name_file_load_plot{1} = 'N2_Tbw_94.6425_Ttw_189.285_Pb_5093700_PrEc_bulk_isothermal_Ke';
% name_file_load_plot{2} = 'N2_Tbw_113.571_Ttw_138.809_Pb_5093700_PrEc_bulk_isothermal_Ke';
name_file_load_plot{2} = 'N2_Tbw_94.6425_Ttw_189.285_Pb_6791600_PrEc_bulk_isothermal_Ke';
[~,~,~,~,W1{1},W2{1},W3{1},P_V{1},~,~,~,~,~,~,~] = CalculateEnergybudget(name_file_load_plot{1},0, 0.6, 10000, 0.1, bSolver,HP_model,Fluid.Substance);
[~,~,~,~,W1{2},W2{2},W3{2},P_V{2},~,~,~,~,~,~,~] = CalculateEnergybudget(name_file_load_plot{2},0, 0.6, 10000, 0.1, bSolver,HP_model,Fluid.Substance);
test_label{1} = '$NI-1$';
test_label{2} = '$NI-5$';
plot_EnergyBudgetInstability(y,delta,W1,W2,W3,P_V,test_label,save_label)


%% Compressibility effects
% I-2 Br = 0.1 vs. 0.5
save_label = 'I_2_Br_0_01_0_5';
name_file_load_plot{1} = 'N2_Tbw_119.8805_Ttw_119.8805_Pb_5093700_PrEc_bulk_isothermal_Ke';
name_file_load_plot{2} = 'N2_Tbw_119.8805_Ttw_119.8805_Pb_5093700_PrEc_bulk_isothermal_Ke';
clear q_rho,clear q_u, clear q_v,clear q_w,clear q_T
[~,~,~,~,~,~,~,~,q_rho{1},q_u{1},q_v{1},q_w{1},q_T{1},q_P{1},Ma_g{1}] = CalculateEnergybudget(name_file_load_plot{1},0, 1, 10000, 0.01, bSolver,HP_model,Fluid.Substance);
[~,~,~,~,~,~,~,~,q_rho{2},q_u{2},q_v{2},q_w{2},q_T{2},q_P{2},Ma_g{2}] = CalculateEnergybudget(name_file_load_plot{1},0, 1, 10000, 0.5, bSolver,HP_model,Fluid.Substance);
test_label{1} = '$Br = 0.01$';
test_label{2} = '$Br = 0.5$';
plot_EnergyBudgetCompressibility(y,delta,q_u,q_v,q_P,Ma_g,test_label,save_label)



% I-1 Br = 0.1 vs. 0.5
save_label = 'I_1_Br_0_01_0_5';
name_file_load_plot{1} = 'N2_Tbw_94.6425_Ttw_94.6425_Pb_5093700_PrEc_bulk_isothermal_Ke';
name_file_load_plot{2} = 'N2_Tbw_94.6425_Ttw_94.6425_Pb_5093700_PrEc_bulk_isothermal_Ke';
clear q_rho,clear q_u, clear q_v,clear q_w,clear q_T
[~,~,~,~,~,~,~,~,q_rho{1},q_u{1},q_v{1},q_w{1},q_T{1},q_P{1},Ma_g{1}] = CalculateEnergybudget(name_file_load_plot{1},0, 1, 10000, 0.01, bSolver,HP_model,Fluid.Substance);
[~,~,~,~,~,~,~,~,q_rho{2},q_u{2},q_v{2},q_w{2},q_T{2},q_P{2},Ma_g{2}] = CalculateEnergybudget(name_file_load_plot{1},0, 1, 10000, 0.5, bSolver,HP_model,Fluid.Substance);
test_label{1} = '$Br = 0.01$';
test_label{2} = '$Br = 0.5$';
plot_EnergyBudgetCompressibility(y,delta,q_u,q_v,q_P,Ma_g,test_label,save_label)



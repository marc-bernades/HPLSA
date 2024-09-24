% %% Linear Stability Theory Solver %%
clear all; close all; clc
% format long g
% Solver and model
bSolver           = 'Real';
HP_model          = 'CoolProp';  % CoolProp, RefProp, Constant, Power, LowPressure, HighPressure   
% pyversion('/usr/bin/python3.8')
% [v,e] = pyversion; system([e,' -m pip install --user -U CoolProp'])
% system([e,' -m pip install --user -U CoolProp'])

% RefProp initialization
if strcmp(HP_model,"RefProp")
    py.CoolProp.CoolProp.set_config_string(py.CoolProp.CoolProp.ALTERNATIVE_REFPROP_PATH,'/home/marc/Documents/REFPROP-cmake/')
    py.CoolProp.CoolProp.PropsSI("D","T",295,"P",8E6,"REFPROP::CO2");
end

% Fluid properties (ideal-gas)
Fluid.R_specific  = 188.9;
Fluid.gamma       = 1.289;
Fluid.mu_0        = 1;
Fluid.kappa_0     = 1;

%% Initial conditions
% Setup initial conditions
% Fluid.Substance         = 'N2';
% [~, T_c, P_c, ~, rho_c, ~,  ~,  ~,  ~,  ~,  ~, ~,  ~,  ~,  ~,  ~,  ~,  ~] = Substance_Library(Fluid.Substance);
% T_bw  = 0.8*T_c; %0.8*T_c; %0.9425*T_c;
% T_tw  = 0.8*T_c; %0.8*T_c; %0.9425*T_c;
% P_b   = 2*P_c; %2*P_c; %1.0844*P_c;

Fluid.Substance         = 'N2';
[~, T_c, P_c, ~, rho_c, ~,  ~,  ~,  ~,  ~,  ~, ~,  ~,  ~,  ~,  ~,  ~,  ~] = Substance_Library(Fluid.Substance);
T_bw  = 0.75*T_c; 
T_tw  = 1.5*T_c; 
P_b   = 0.03*P_c;
% T_bw  = 310; 
% T_tw  = 310; 
% P_b   = 8E6;

% Grid points (y)
N = 200;

% Reference legnth scale
delta = 1; % Half channel-height
L_y   = 2*delta;
y_0   = 0; % Bottom y-position

%% Reference scaling
bScaling = 'bulk'; % wall, bulk

%% Grid
q_vars = 5;
% Chevysheb collocation Matrix
[y,D] = Chevysheb_collocation(N,delta); %dy = D*y = CentralDerivative_d1_2ndOrder_1D(y,y)

% Uniform with boundary cells
% y   = (y_0 + -0.5*L_y/N):L_y/N:(L_y + y_0 + 0.5*L_y/N); 

% Uniform with cells at boundary
% y   = y_0:L_y/(N-1):(L_y + y_0); 

%% Baseflow
bTarget       = 'PrEc';
% N_target      = 0.1; % [1E-2 5E-2 0.1 0.25 0.5]; % Isothermal Transcritical
% N_target      = [1E-2 5E-2 0.1]; % Non-Isothermal Transcritical
N_target      = 5.6*1E-6; % Non-isothermal U_r = 1m/s

% N_target      = [0.01 0.03 0.05 0.07]; %1E-5; % CO_2 isothermal limit
% N_target      =  [1E-3 1E-2 5E-2 0.1 0.2 0.5]; % CO2 bulk subcritical / transcritical 0.5 not necessary / supercritical 5E-2 Re > 10000

% n_LST_Sweep   = 1000; %1000:1000:10000; %200:200:8000;% 1000:1000:10000;
% n_LST_Sweep   = [100, 500, 1000, 2000:1000:10000]; % For TG
n_LST_Sweep   = 500:200:7000;% 1000:1000:10000; % For Critical Reynolds TG

BF        = cell(1,length(N_target));
Dy        = cell(1,length(N_target));

for ii = 1:length(N_target)
    
    disp("Computing Baseflow " + bTarget + " = "  + num2str(N_target(ii)) + ", Sweep " + num2str(ii) + " out of " + num2str(length(N_target)))
    
    [BF{ii}, Dy{ii}] = Baseflow_LST_adimensional(y,D,N,delta,L_y,T_bw,T_tw,P_b,bSolver,HP_model, Fluid.Substance, Fluid,bTarget,bScaling,N_target(ii));
    
    disp(" Resulted Reynolds " + bScaling + " = " + num2str(BF{ii}.Re) + newline) 
end

%% Linear stability theory (LST)

% Wavenumbers
omega = 0; % Frequency only for spatial LST
% alpha = 0.6:0.2:1.6; % Streamwise wavenumber subrcritical
% alpha = 0.6:0.2:1.6; % Streamwise wavenumber transcritical
% alpha = 0.7:0.2:1.3; % 0.7:0.1:1.2 Streamwise wavenumber supercritical
% alpha = 0.4:0.2:1.2; % Non-isothermal
% alpha = 0.0:0.5:2.0; % Streamwise wavenumber (3D perturbation)
% beta  = 0.0:1.0:5.0; % Spanwise wavenumber (3D perturbation)
% beta  = 0;
% Single alpha study
alpha = 0;
beta  = 2;

% Initialize
A            = cell(length(alpha),length(beta),length(N_target),length(n_LST_Sweep));
B            = cell(length(alpha),length(beta),length(N_target),length(n_LST_Sweep));
Vec_eigen    = cell(length(alpha),length(beta),length(N_target),length(n_LST_Sweep));
D_Val_eigen  = cell(length(alpha),length(beta),length(N_target),length(n_LST_Sweep));
Val_eigen    = cell(length(alpha),length(beta),length(N_target),length(n_LST_Sweep));
G_tot        = cell(length(alpha),length(beta),length(N_target),length(n_LST_Sweep));
G            = cell(length(alpha),length(beta),length(N_target),length(n_LST_Sweep));
G_max        = zeros(length(alpha),length(beta),length(N_target),length(n_LST_Sweep));
E_norm       = cell(length(alpha),length(beta),length(N_target),(N+1)*5);
GG           = cell(length(alpha),length(beta),length(N_target),length(n_LST_Sweep));
GG_max       = zeros(length(alpha),length(beta),length(N_target),length(n_LST_Sweep));
t_all        = cell(length(alpha),length(beta),length(N_target),length(n_LST_Sweep));
q_in         = cell(length(alpha),length(beta),length(N_target),length(n_LST_Sweep));
q_out        = cell(length(alpha),length(beta),length(N_target),length(n_LST_Sweep));
Val_eigen_tg = cell(length(alpha),length(beta),length(N_target),length(n_LST_Sweep));
Vec_eigen_tg = cell(length(alpha),length(beta),length(N_target),length(n_LST_Sweep));


for bb = 1:length(beta)
    for aa = 1:length(alpha)
        for ii = 1:length(N_target)
            for jj = 1:length(n_LST_Sweep)

                k = sqrt(alpha(aa).^2 + beta(bb).^2);

                disp("Computing LST " + "beta" + " = "  + num2str(beta(bb))  + " and " + "alpha" + " = "  + num2str(alpha(aa)) + " and " + bTarget + " = "  + num2str(N_target(ii)) + " Re = " + num2str(n_LST_Sweep(jj)) + ", Sweep " + num2str((bb-1)*length(alpha)*length(N_target)*length(n_LST_Sweep) + (aa-1)*length(N_target)*length(n_LST_Sweep) + (ii-1)*length(n_LST_Sweep) + jj) + " out of " + num2str(length(beta)*length(alpha)*length(N_target)*length(n_LST_Sweep)))

                %% Operator
                [A{aa,bb,ii,jj},B{aa,bb,ii,jj}] = LST(delta,N,BF{ii},Dy{ii},omega,beta(bb),alpha(aa),bSolver,HP_model, Fluid.Substance, Fluid, bScaling, n_LST_Sweep(jj));

                %% EigenProblem
                % Option A - Eig unstable for transient growth
                % [Vec_eigen{aa,bb,ii,jj}, D_Val_eigen{aa,bb,ii,jj}] = eig(A{aa,bb,ii,jj},B{aa,bb,ii,jj},'qz');
                % Val_eigen{aa,bb,ii,jj}                             = eig(A{aa,bb,ii,jj},B{aa,bb,ii,jj},'qz');

                % Option B - Eigs iterative method
                % Pre-conditioning
                % A{aa,bb,ii,jj} = sparse(A{aa,bb,ii,jj});
                % B{aa,bb,ii,jj} = sparse(B{aa,bb,ii,jj});
                N_modes    = 350; %length(y)*q_vars; % Modes to compute
                opts.v0    = rand(length(y)*q_vars,1);
                opts.tol   = 1e-15;
                opts.maxit = 1E10;
                % opts.p     = 2*N_modes+1;
                sigma      = 0; % Shift value (e.g., around which you want to find eigenvalues)
                % parpool; % Start parallel pool
                [Vec_eigen{aa,bb,ii,jj}, D_Val_eigen{aa,bb,ii,jj}, flag] = eigs(A{aa,bb,ii,jj},B{aa,bb,ii,jj},N_modes,sigma,opts);
                Val_eigen{aa,bb,ii,jj}                             = eigs(A{aa,bb,ii,jj},B{aa,bb,ii,jj},N_modes,sigma,opts);
                % delete(gcp); % Close parallel pool
                if flag == 0
                    disp("All eigenvalues converged...")
                else
                    disp("Some eigenvalues did NOT converged...")
                end


                % Order eigenvalues imaginary part decreasing
                [Val_eigen_dec, i_sort] = sort(-imag(Val_eigen{aa,bb,ii,jj} ));
                Val_eigen_temp          = Val_eigen{aa,bb,ii,jj};
                Vec_eigen_temp          = Vec_eigen{aa,bb,ii,jj};
                Val_eigen{aa,bb,ii,jj}  = Val_eigen_temp(i_sort);
                Vec_eigen{aa,bb,ii,jj}  = Vec_eigen_temp(:,i_sort);


                % Growth rate modal-analysis
                G_tot{aa,bb,ii,jj} = imag(k.*Vec_eigen{aa,bb,ii,jj});
                G{aa,bb,ii,jj}     = imag(k*Val_eigen{aa,bb,ii,jj});

                try
                    G_max(aa,bb,ii,jj) = max(G{aa,bb,ii,jj}(~G{aa,bb,ii,jj}==0));
                catch
                    disp("Modes from alpha = " + num2str(alpha(aa)) + " and " + bTarget + " = "  + num2str(N_target(ii)) + " and Re = " + num2str(n_LST_Sweep(jj)) + " are null...")
                    G_max(aa,bb,ii,jj) = max(G{aa,bb,ii,jj});
                end

                % Spectrum check
                for n_var = 1:q_vars
                    Val_eigen_var{n_var} = Val_eigen{aa,bb,ii,jj}(n_var:q_vars:end);
                    Vec_eigen_var{n_var} = Vec_eigen{aa,bb,ii,jj}(n_var:q_vars:end,:); % Take all columns for each perturbation row
                end
                % plot_Spectrum(y, delta, N, alpha, N_target,n_LST_Sweep, Val_eigen{aa,bb,ii,jj},  Vec_eigen_var, 'test')

                %% Transient growth
                q_vars = 5;
                t_max = 400;
                N_t   = 200;

                % Spectrum selected range
                Imag_max = 0.1;
                Imag_min = -2.0;
                Real_max = 3.0;
                Real_min = -1.0;
                % Imag_max = 10^6*Imag_max;
                % Imag_min = 10^6*Imag_min;
                % Real_max = 10^6*Real_max;
                % Real_min = 10^6*Real_min;

                % Weightings
                m_d = 1;
                m_T = 1;
                % m_d = BF{ii}.T./(BF{ii}.rho.*BF{ii}.gamma.*BF{ii}.Ma.^2);
                % m_T = 1./(BF{ii}.gamma.*(BF{ii}.gamma - 1).*BF{ii}.T.*BF{ii}.Ma.^2);
                [GG{aa,bb,ii,jj}, t_all{aa,bb,ii,jj}, GG_max(aa,bb,ii,jj), q_in{aa,bb,ii,jj}, q_out{aa,bb,ii,jj}, ...
                 Val_eigen_tg{aa,bb,ii,jj}, Vec_eigen_tg{aa,bb,ii,jj}] = LST_TransientGrowth( length(y),y, q_vars, Val_eigen{aa,bb,ii,jj}, Vec_eigen{aa,bb,ii,jj}, ...
                                                                                                  Imag_max, Imag_min, Real_max, Real_min, N_t, t_max, m_d, m_T, N_modes );

                % Check alpha and beta = 0
                if alpha(aa) == 0 && beta(bb) == 0
                    GG{aa,bb,ii,jj}     = 0*GG{aa,bb,ii,jj};
                    GG_max(aa,bb,ii,jj) = 0;
                end

                % Spectrum check
                q_vars = 5;
                for n_var = 1:q_vars
                    Val_eigen_var{n_var} = Val_eigen_tg{aa,bb,ii,jj}(n_var:q_vars:end);
                    Vec_eigen_var{n_var} = Vec_eigen_tg{aa,bb,ii,jj}(n_var:q_vars:end,:); % Take all columns for each perturbation row
               end
                % plot_Spectrum(y, delta, N, alpha, N_target,n_LST_Sweep, Val_eigen_tg{aa,bb,ii,jj},  Vec_eigen_var, 'test_reduced')

%                 figure
%                 drawnow
%                 plot( t_all{aa,bb,ii,jj}, GG{aa,bb,ii,jj}); hold on
%                 aaa = 1;


            end
        end
    end
end

%% Save results
name_label = 'Non_Isothermal_TG_Ur_1_alpha_0_beta_2';
% name_label = 'Non_Isothermal_TG_Ur_1_alpha_0_beta_2';
% name_label = 'Isothermal_TG_CoolProp_2D_Re_Sweep_Test';
% name_label = 'Isothermal_Limit_TG_CoolProp';


LST_output_TransientGrowth(Fluid, T_bw, T_tw, T_c, P_b, P_c, N, y, D, delta, bScaling, bTarget,N_target, n_LST_Sweep, alpha, beta, BF,Dy,A,B,Vec_eigen,Val_eigen, G, G_max,  GG, t_all, GG_max, q_in, q_out, name_label);
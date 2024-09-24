% %% Linear Stability Theory Solver %%
clear all; close all; clc

% Solver and model
bSolver           = 'Real';
HP_model          = 'CoolProp';  % CoolProp, Constant, Power, LowPressure, HighPressure   
% pyversion('/usr/bin/python3.8')
% [v,e] = pyversion; system([e,' -m pip install --user -U CoolProp'])
% system([e,' -m pip install --user -U CoolProp'])

% RefProp initialization
if strcmp(HP_model,"RefProp")
    py.CoolProp.CoolProp.set_config_string(py.CoolProp.CoolProp.ALTERNATIVE_REFPROP_PATH,'/home/marc/Documents/REFPROP-cmake/')
    py.CoolProp.CoolProp.PropsSI("D","T",295,"P",8E6,"REFPROP::CO2");
end

% Fluid properties
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
T_bw  = 1.5*T_c; 
T_tw  = 1.5*T_c; 
P_b   = 1.5*P_c;
% T_bw  = 290; 
% T_tw  = 290; 
% P_b   = 8E6;

% Grid points (y)
N = 200;

% Reference legnth scale
delta = 1; % Half channel-height
L_y   = 2*delta;
y_0   = 0; % Bottom y-position

%% Reference scaling
bScaling = 'bulk'; % wall, bulk

%% Algorithm choice
bMODAL = 1; % For Modal analysis, 0 for energy-based analysis

%% Grid
% Chevysheb collocation Matrix
[y,D] = Chevysheb_collocation(N,delta); %dy = D*y = CentralDerivative_d1_2ndOrder_1D(y,y)

% Uniform with boundary cells
% y   = (y_0 + -0.5*L_y/N):L_y/N:(L_y + y_0 + 0.5*L_y/N); 

% Uniform with cells at boundary
% y   = y_0:L_y/(N-1):(L_y + y_0); 

%% Baseflow
bTarget       = 'PrEc';
N_target      = 0.01; %[1E-2 5E-2 0.1 0.25 0.5]; % [1E-4 5E-2 0.1 0.25 0.5]; %[1E-4 1E-2 5E-2 0.1 0.25 0.5]; % Isothermal Transcritical
% N_target      = [1E-2 5E-2 0.1]; %[1E-4 1E-2 5E-2 0.1]; % Non-Isothermal Transcritical
% N_target      = 5.6E-6; % U_r = 1m/s NI-5


% N_target      = 1E-5; % CO_2 isothermal limit
% N_target      =  [1E-3 1E-2 5E-2 0.1 0.2 0.5]; % CO2 bulk subcritical / transcritical 0.5 not necessary / supercritical 5E-2 Re > 10000
% N_target = [0.01 0.03 0.05 0.07];

n_LST_Sweep   = [100, 500, 1000, 2000:1000:10000];


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
beta  = 0; % Spanwise wavenumber (2D perturbation)
alpha = 0.6:0.2:1.6; % Streamwise wavenumber subcritical
% alpha = 0.6:0.2:1.6; % Streamwise wavenumber transcritical
% alpha = 0.7:0.1:1.2; % Streamwise wavenumber supercritical
% alpha = 0.4:0.2:1.2; % Non-isothermal
% alpha = 1.0; % Ke

% alpha = 0.6:0.2:1.4; % CO2 bulk

% Initialize
A            = cell(length(alpha),length(N_target),length(n_LST_Sweep));
B            = cell(length(alpha),length(N_target),length(n_LST_Sweep));
Vec_eigen    = cell(length(alpha),length(N_target),length(n_LST_Sweep));
D_Val_eigen  = cell(length(alpha),length(N_target),length(n_LST_Sweep));
Val_eigen    = cell(length(alpha),length(N_target),length(n_LST_Sweep));
G_tot        = cell(length(alpha),length(N_target),length(n_LST_Sweep));
G            = cell(length(alpha),length(N_target),length(n_LST_Sweep));
G_max        = zeros(length(alpha),length(N_target),length(n_LST_Sweep));


for aa = 1:length(alpha)
    for ii = 1:length(N_target)
        for jj = 1:length(n_LST_Sweep)


            disp("Computing LST " + "alpha" + " = "  + num2str(alpha(aa)) + " and " + bTarget + " = "  + num2str(N_target(ii)) + " Re = " + num2str(n_LST_Sweep(jj)) + ", Sweep " + num2str((aa-1)*length(N_target)*length(n_LST_Sweep) + (ii-1)*length(n_LST_Sweep) + jj) + " out of " + num2str(length(alpha)*length(N_target)*length(n_LST_Sweep)))

            [A{aa,ii,jj},B{aa,ii,jj}] = LST(delta,N,BF{ii},Dy{ii},omega,beta,alpha(aa),bSolver,HP_model, Fluid.Substance, Fluid, bScaling, n_LST_Sweep(jj));

            if bMODAL == 1
                %% MODAL ANALYSIS - Option A - Eigen values
                [Vec_eigen{aa,ii,jj}, D_Val_eigen{aa,ii,jj}] = eig(A{aa,ii,jj},B{aa,ii,jj});
                Val_eigen{aa,ii,jj}                          = eig(A{aa,ii,jj},B{aa,ii,jj});

            else

                %% ENERGY BUDGET - Option B - Eigs iterative method
                % Pre-conditioning
                % A{aa,ii,jj} = sparse(A{aa,ii,jj});
                % B{aa,ii,jj} = sparse(B{aa,ii,jj});
                q_vars     = 5;
                N_modes    = length(y)*q_vars; % Modes to compute
                rng('default');
                opts.v0    = rand(length(y)*q_vars,1);
                opts.tol   = 1e-15;
                opts.maxit = 1E10;
                % opts.p     = 2*N_modes+1;
                sigma      = 0.2; % Shift value (e.g., around which you want to find eigenvalues)
                % parpool; % Start parallel pool
                [Vec_eigen{aa,ii,jj}, D_Val_eigen{aa,ii,jj}, flag] = eigs(A{aa,ii,jj},B{aa,ii,jj},N_modes,sigma,opts);
                Val_eigen{aa,ii,jj}                                = eigs(A{aa,ii,jj},B{aa,ii,jj},N_modes,sigma,opts);
                % delete(gcp); % Close parallel pool
                if flag == 0
                    disp("All eigenvalues converged...")
                else
                    disp("Some eigenvalues did NOT converged...")
                end
            end


            % Growth rate
            G_tot{aa,ii,jj} = imag(alpha(aa).*Vec_eigen{aa,ii,jj});
            G{aa,ii,jj}     = imag(alpha(aa).*Val_eigen{aa,ii,jj});

            try
                G_max(aa,ii,jj) = max(G{aa,ii,jj}(~G{aa,ii,jj}==0));
            catch
                disp("Modes from alpha = " + num2str(alpha(aa)) + " and " + bTarget + " = "  + num2str(N_target(ii)) + " and Re = " + num2str(n_LST_Sweep(jj)) + " are null...")
                G_max(aa,ii,jj) = max(G{aa,ii,jj});
            end


        end
    end
end

%% Save results
TestFolder = 'Modal'; % Modal, Non_Modal, EnergyBudget, ECCOMAS, Verification, Test
name_label = 'isothermal_Br_0_01';
% name_label = 'isothermal_Ke';
% name_label = 'non_isothermal_Ur_1';

LST_output(Fluid, T_bw, T_tw, T_c, P_b, P_c, N, y, D, delta, bScaling, bTarget,N_target, n_LST_Sweep, alpha, beta, BF,Dy,Vec_eigen,Val_eigen, G, G_max, TestFolder, name_label);
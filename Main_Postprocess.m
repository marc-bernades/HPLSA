%% Load results and calculations
clear all; clc; close all;
% Select fluid and initial conditions
Fluid.Substance         = 'N2';
[~, T_c, P_c, ~, rho_c, ~,  ~,  ~,  ~,  ~,  ~, ~,  ~,  ~,  ~,  ~,  ~,  ~] = Substance_Library(Fluid.Substance);

% File name
% name_file_load = 'Verification/CO2_Tbw_290_Ttw_290_Pb_8000000_PrEc_wall_1E-5_RG_CP_Re_10000';
% name_file_load = 'CO2_Tbw_290_Ttw_290_Pb_8000000_PrEc_wall_isothermal';
% name_file_load = 'Verification/CO2_Tbw_290_Ttw_290_Pb_8000000_PrEc_bulk_isothermal_limit';

% name_file_load =  'N2_Tbw_100.952_Ttw_100.952_Pb_6791600_PrEc_bulk_1E-5_RG_CP_Re_10000';
% name_file_load =  'N2_Tbw_94.6425_Ttw_94.6425_Pb_6791600_PrEc_bulk_isothermal';
% name_file_load =  'N2_Tbw_119.8805_Ttw_119.8805_Pb_6791600_PrEc_bulk_isothermal';

% % Isothermal flow cases (I1-3)
% name_file_load =  'N2_Tbw_94.6425_Ttw_94.6425_Pb_5093700_PrEc_bulk_isothermal';
% name_file_load =  'N2_Tbw_119.8805_Ttw_119.8805_Pb_5093700_PrEc_bulk_isothermal';
% name_file_load =  'N2_Tbw_189.285_Ttw_189.285_Pb_5093700_PrEc_bulk_isothermal';
% name_file_load =  'N2_Tbw_94.6425_Ttw_94.6425_Pb_5093700_PrEc_bulk_isothermal_Br_0_01';
% name_file_load =  'N2_Tbw_119.8805_Ttw_119.8805_Pb_5093700_PrEc_bulk_isothermal_Br_0_01';
name_file_load =  'N2_Tbw_189.285_Ttw_189.285_Pb_5093700_PrEc_bulk_isothermal_Br_0_01';


% % Non-Isothermal flow cases (NI-1-6)
% name_file_load =  'N2_Tbw_94.6425_Ttw_189.285_Pb_5093700_PrEc_bulk_non_isothermal';
% name_file_load =  'N2_Tbw_94.6425_Ttw_189.285_Pb_16979000_PrEc_bulk_non_isothermal';
% name_file_load =  'N2_Tbw_113.571_Ttw_138.809_Pb_5093700_PrEc_bulk_non_isothermal';
% name_file_load =  'N2_Tbw_94.6425_Ttw_189.285_Pb_101874_PrEc_bulk_non_isothermal';
% name_file_load =  'N2_Tbw_94.6425_Ttw_189.285_Pb_6791600_PrEc_bulk_non_isothermal_Ur_1';
% name_file_load =  'N2_Tbw_94.6425_Ttw_189.285_Pb_101874_PrEc_bulk_non_isothermal_Ur_1';


% Others
% name_file_load =  'N2_Tbw_94.6425_Ttw_189.285_Pb_33958_PrEc_bulk_non_isothermal';
% name_file_load =  'N2_Tbw_94.6425_Ttw_94.6425_Pb_101874_PrEc_bulk_non_isothermal';

% name_file_load =  'CO2_Tbw_290_Ttw_290_Pb_8000000_PrEc_wall_Non_Modal';

% name_file_load =  'EnergyBudget/N2_Tbw_119.8805_Ttw_119.8805_Pb_5093700_PrEc_wall_Non_isothermal_Ke';

% Check isothermal setup boolean
bIsothermal = ~contains(name_file_load,'non_isothermal'); 

% name_file_load = strcat(Fluid.Substance, '_Tbw_', num2str(T_bw), '_Ttw_', num2str(T_tw), '_Pb_', num2str(P_b), '_', bTarget, '_', bScaling, '_', name_label);
Results = load(strcat('Results/Modal/', name_file_load, '.mat'));

% Unpack structure and load variables
struct_field_names = fieldnames(Results.Data);
for k = 1:length(struct_field_names)
      eval([struct_field_names{k}, ' = Results.Data.',struct_field_names{k}, ';']);
end

%% Dimensionless numbers
% plot_DimensionlessNumbers(BF, N_target, name_file_load)

%% Plot Baseflow
for ii = 1:length(BF)
    disp("Base flow Br = " + num2str(BF{ii}.Pr*BF{ii}.Ec) + " and Ma = " + num2str(BF{ii}.Ma) + " u = " + num2str(BF{ii}.norm.u) + newline)
end
plot_Baseflow_LST(y,delta,BF,N_target,T_c, name_file_load);

% Obtion to plot only one Re
% PrEc = 0.1;
% [~,pos] = min(abs(N_target - PrEc));
% plot_Baseflow_LST(y,delta,BF,N_target,name_file_load, pos);

% % Compare to baseflows
% T_bw2  = 0.75*T_c;
% T_tw2  = 1.5*T_c;
% P_b2   = 2*P_c;
% delta2 = 1E-4;
% 
% % File name
% name_file_load_2 = strcat(Fluid.Substance, '_Tbw_', num2str(T_bw2/T_c), '_Ttw_', num2str(T_tw2/T_c), '_Pb_', num2str(P_b2/P_c), '_', bTarget, '_', bScaling);
% plot_Baseflow_LST_Comparison(y,delta,BF,N_target, PrEc, name_file_load, name_file_load_2);


%% Modes
% Select Alpha and Reynolds
for aa = 1:length(alpha)
    for ii = 1:length(N_target)
        for jj = 1:length(n_LST_Sweep)
            realc = real(Val_eigen{aa,ii,jj});
            imagc = imag(Val_eigen{aa,ii,jj});
            Cond = abs(realc) <= 10;
            disp(['Alpha = ', num2str(alpha(aa)), ', ', bTarget, ' = ', num2str(N_target(ii)), ' Re = ', num2str(n_LST_Sweep(jj)) , ': Unstable mode = ', num2str(alpha(aa)*max(imagc(Cond))), ' total of = ', num2str(sum(imagc(Cond)>0))])
        end
    end
end

%% Select alpha and N_plot for Spectrum and Perturbation
alpha_plot = 1.0;
PrEc_plot  = 0.01;
Re_plot    = 10000;

% Breakdown eigen vector for each variable
q_vars = 5;
Val_eigen_var = cell(1,q_vars);
Vec_eigen_var = cell(1,q_vars);

[~,aa] = min(abs(alpha-alpha_plot));
[~,ii] = min(abs(N_target-PrEc_plot));
[~,jj] = min(abs(n_LST_Sweep-Re_plot));

for n_var = 1:q_vars
    Val_eigen_var{n_var} = Val_eigen{aa,ii,jj}(n_var:q_vars:end);
    Vec_eigen_var{n_var} = Vec_eigen{aa,ii,jj}(n_var:q_vars:end,:); % Take all columns for each perturbation row
end

%% Spectrum at given alpha and Re
plot_Spectrum(y, delta, N, alpha(aa), N_target(ii), n_LST_Sweep(jj), Val_eigen{aa,ii,jj},  Vec_eigen_var, name_file_load)

%% Perturbation at given mode
realc = real(Val_eigen{aa,ii,jj});
imagc = imag(Val_eigen{aa,ii,jj});
Cond   = (abs(realc) <= 10) & (realc < 0.5);
Val_eigen_target   = max(imagc(Cond));
plot_Perturbation(y, delta, N, alpha(aa), N_target(ii),  n_LST_Sweep(jj), Val_eigen{aa,ii,jj}, Val_eigen_target, Vec_eigen_var, name_file_load)
plot_Perturbation_inset(y, delta, N, alpha(aa), N_target(ii),  n_LST_Sweep(jj), Val_eigen{aa,ii,jj}, Val_eigen_target, Vec_eigen_var, name_file_load)

if bIsothermal
    N_plot = N_target(N_target>1E-3 & N_target<=0.5); % Isothermal N_target<=0.5 sub, N_target<=0.5 trans, N_target<=0.05 super
else
    N_plot = N_target(N_target>=1E-2 & N_target<=0.1); % Non-isothermal
end
plot_Perturbation_Sweep(y, delta, N, alpha, N_plot, n_LST_Sweep, alpha_plot, Re_plot, Val_eigen, Vec_eigen, Val_eigen_var, Vec_eigen_var, name_file_load)

Val_eigen_target   = -0.138227; %max(imagc(Cond)); % Most unstable mode > Max %-0.138227;
plot_Perturbation(y, delta, N, alpha(aa), N_target(ii),  n_LST_Sweep(jj), Val_eigen{aa,ii,jj}, Val_eigen_target, Vec_eigen_var, name_file_load)

%% Stability map
if bIsothermal
    N_plot = N_target(N_target<=0.5); % Isothermal N_target<=0.5 sub, N_target<=0.5 trans, N_target<=0.05 super
else
    if length(N_target) > 1
        N_plot = N_target(N_target>=0.01 & N_target<=0.1); % Non-isothermal
    else
        N_plot = N_target;
    end
end
plot_StabilityDiagram(alpha, N_plot, N_target, n_LST_Sweep, G_max, name_file_load, bIsothermal)

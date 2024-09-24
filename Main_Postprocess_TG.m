%% Load results and calculations
clear all; clc; close all;
% Select fluid and initial conditions
Fluid.Substance         = 'N2';
[~, T_c, P_c, ~, rho_c, ~,  ~,  ~,  ~,  ~,  ~, ~,  ~,  ~,  ~,  ~,  ~,  ~] = Substance_Library(Fluid.Substance);

%% File name
% % Iosthermal flow cases validation CO2 wall-norm
% name_file_load =  'CO2_Tbw_290_Ttw_290_Pb_8000000_PrEc_wall_Isothermal_TG_CoolProp';
% name_file_load =  'CO2_Tbw_300_Ttw_300_Pb_8000000_PrEc_wall_Isothermal_TG_CoolProp';
% name_file_load =  'CO2_Tbw_310_Ttw_310_Pb_8000000_PrEc_wall_Isothermal_TG_CoolProp';
name_file_load =  'CO2_Tbw_290_Ttw_290_Pb_8000000_PrEc_wall_Isothermal_Limit_TG_CoolProp';

%% Scaling with lambda to mu
% % Isothermal flow cases (I1-3)
% name_file_load =  'N2_Tbw_94.6425_Ttw_94.6425_Pb_5093700_PrEc_bulk_Isothermal_TG';
% name_file_load =  'N2_Tbw_119.8805_Ttw_119.8805_Pb_5093700_PrEc_bulk_Isothermal_TG';
% name_file_load =  'N2_Tbw_189.285_Ttw_189.285_Pb_5093700_PrEc_bulk_Isothermal_TG';

% % Non-Isothermal flow cases (NI1-6)
% name_file_load =  'N2_Tbw_94.6425_Ttw_189.285_Pb_5093700_PrEc_bulk_Non_Isothermal_TG';
% name_file_load =  'N2_Tbw_94.6425_Ttw_189.285_Pb_16979000_PrEc_bulk_Non_Isothermal_TG';
% name_file_load =  'N2_Tbw_113.571_Ttw_138.809_Pb_5093700_PrEc_bulk_Non_Isothermal_TG';
% name_file_load =  'N2_Tbw_94.6425_Ttw_189.285_Pb_6791600_PrEc_bulk_Non_Isothermal_TG_Ur_1';
% name_file_load =  'N2_Tbw_94.6425_Ttw_189.285_Pb_101874_PrEc_bulk_Non_Isothermal_TG_Ur_1';

%% Scaling lamda* w/mu_b Re sweep
% % Isothermal flow cases (I1-3)
% name_file_load =  'N2_Tbw_94.6425_Ttw_94.6425_Pb_5093700_PrEc_bulk_Isothermal_TG_2D_Re_b';
% name_file_load =  'N2_Tbw_119.8805_Ttw_119.8805_Pb_5093700_PrEc_bulk_Isothermal_TG_2D_Re_b';
% name_file_load =  'N2_Tbw_189.285_Ttw_189.285_Pb_5093700_PrEc_bulk_Isothermal_TG_2D_Re_b';

% % Non-Isothermal flow cases (NI1-6)
% name_file_load =  'N2_Tbw_94.6425_Ttw_189.285_Pb_5093700_PrEc_bulk_Non_Isothermal_TG_2D_Re_b';
% name_file_load =  'N2_Tbw_94.6425_Ttw_189.285_Pb_16979000_PrEc_bulk_Non_Isothermal_TG_2D_Re_b';
% name_file_load =  'N2_Tbw_113.571_Ttw_138.809_Pb_5093700_PrEc_bulk_Non_Isothermal_TG_2D_Re_b';
% name_file_load =  'N2_Tbw_113.571_Ttw_138.809_Pb_101874_PrEc_bulk_Non_Isothermal_TG_2D_Re_b';
% name_file_load =  'N2_Tbw_94.6425_Ttw_189.285_Pb_6791600_PrEc_bulk_Non_Isothermal_TG_2D_Ur_1_Re_b';
% name_file_load =  'N2_Tbw_94.6425_Ttw_189.285_Pb_101874_PrEc_bulk_Non_Isothermal_TG_2D_Ur_1_Re_b';

%% Scaling lamda* Specific sweeps
% % Isothermal flow cases (I1-3)
% name_file_load =  'N2_Tbw_94.6425_Ttw_94.6425_Pb_5093700_PrEc_bulk_Isothermal_TG_alpha_0_beta_2';
% name_file_load =  'N2_Tbw_94.6425_Ttw_94.6425_Pb_5093700_PrEc_bulk_Isothermal_TG_alpha_1_beta_1';
% name_file_load =  'N2_Tbw_119.8805_Ttw_119.8805_Pb_5093700_PrEc_bulk_Isothermal_TG_alpha_0_beta_2';
% name_file_load =  'N2_Tbw_119.8805_Ttw_119.8805_Pb_5093700_PrEc_bulk_Isothermal_TG_alpha_1_beta_1';
% name_file_load =  'N2_Tbw_189.285_Ttw_189.285_Pb_5093700_PrEc_bulk_Isothermal_TG_alpha_0_beta_2';
% name_file_load =  'N2_Tbw_189.285_Ttw_189.285_Pb_5093700_PrEc_bulk_Isothermal_TG_alpha_1_beta_1';

% % Non-Isothermal flow cases (NI1-6)
% name_file_load =  'N2_Tbw_94.6425_Ttw_189.285_Pb_5093700_PrEc_bulk_Non_Isothermal_TG_alpha_0_beta_2';
% name_file_load =  'N2_Tbw_94.6425_Ttw_189.285_Pb_16979000_PrEc_bulk_Non_Isothermal_TG_alpha_0_beta_2';
% name_file_load =  'N2_Tbw_113.571_Ttw_138.809_Pb_5093700_PrEc_bulk_Non_Isothermal_TG_alpha_0_beta_2';
% name_file_load =  'N2_Tbw_94.6425_Ttw_189.285_Pb_101874_PrEc_bulk_Non_Isothermal_TG_alpha_0_beta_2';
% name_file_load =  'N2_Tbw_94.6425_Ttw_189.285_Pb_6791600_PrEc_bulk_Non_Isothermal_TG_Ur_1_alpha_0_beta_2';
% name_file_load =  'N2_Tbw_94.6425_Ttw_189.285_Pb_101874_PrEc_bulk_Non_Isothermal_TG_Ur_1_alpha_0_beta_2';


% Unpack structure and load variables
Results = load(strcat('Results/Non_Modal/', name_file_load, '.mat'));
struct_field_names = fieldnames(Results.Data);
for k = 1:length(struct_field_names)
      eval([struct_field_names{k}, ' = Results.Data.',struct_field_names{k}, ';']);
end

% Check isothermal case
bIsothermal = ~contains(name_file_load,'Non_Isothermal'); 

% Select Br to plot I- All, NI- only up to Br <= 0.1
if bIsothermal || contains(name_file_load,'Ur_1')
    N_plot = N_target(N_target<=0.5); % Isothermal N_target<=0.5 sub, N_target<=0.5 trans, N_target<=0.05 super
else
    N_plot = N_target(N_target>1E-5 & N_target<=0.1); % Non-isothermal
end

% Perturbation variables
q_vars = 5;
% Channel half height
delta  = 1;

%% Map Transient Unstable points
GG_unstable = zeros(size(G_max));
for ii = 1:length(N_target)
    for aa = 1:length(alpha)
        for jj = 1:length(n_LST_Sweep)
            for bb = 1:length(beta)
                GG_gradient = gradient(GG{aa,bb,ii,jj});
                GG_smooth   = smooth(GG{aa,bb,ii,jj},10);
                Cond_max    = islocalmax(GG_smooth);
                Cond_unstable = max(GG_smooth(Cond_max)) < GG_smooth(end);
                if isempty(Cond_unstable) %No maxmimum
                    Cond_unstable = true;
                end
                Cond_1      = (GG_smooth(end) > max(GG_smooth)); % I-3
                Cond_2      = (GG_smooth(end) > (GG_smooth(10))); % I-1-2
                Cond_3      = (GG_smooth(end) > (GG_smooth(25))); % NI-1
%                 if Cond_3 && (max(islocalmax(GG{aa,bb,ii,jj})) == 0 || (min(GG{aa,bb,ii,jj}(islocalmax(GG{aa,bb,ii,jj}))) < GG{aa,bb,ii,jj}(end)|| GG_gradient(end) > 10^-1))
                if (max(islocalmax(GG{aa,bb,ii,jj})) == 0) || Cond_unstable
                    GG_unstable(aa,bb,ii,jj) = 1;
                end
            end
        end
    end
end

%% Transient growth
beta_target  = 2.0;
alpha_target = 0.0;
Re_target = 4500;
Br_target = 0.1;
[~,aa] = min(abs(alpha-alpha_target));
[~,bb] = min(abs(beta-beta_target));
[~,ii] = min(abs(N_target-Br_target));
[~,jj] = min(abs(n_LST_Sweep-Re_target));

% Sweep beta
figure
title(strcat('$','Re = ', num2str(n_LST_Sweep(jj)), ', Br = ',  num2str(N_target(ii)), ', \alpha = ',  num2str(alpha(aa)),'$'),'interpreter','latex')
for bb = 1:length(beta)
    plot(t_all{aa,bb,ii,jj}, GG{aa,bb,ii,jj}, 'DisplayName', strcat('$','\beta = ', num2str(beta(bb)), '$')); hold on
end

xlabel('${t}$','interpreter','latex')
ylabel('${G}$','interpreter','latex')

legend('interpreter','latex','Location','northeast','box','off')

set(gca,'linewidth',1.5)
set(gca,'fontsize',14)

% Sweep Reynolds
figure
for jj = 1:length(n_LST_Sweep)
    plot(t_all{aa,bb,ii,jj}, GG{aa,bb,ii,jj}, 'DisplayName', strcat('$','Re = ', num2str(n_LST_Sweep(jj)), '$')); hold on
end

xlabel('${t}$','interpreter','latex')
ylabel('${G}$','interpreter','latex')

legend('interpreter','latex','Location','northeast','box','off')

set(gca,'linewidth',1.5)
set(gca,'fontsize',14)

return

%% Modal counter-part
for n_var = 1:q_vars
    Val_eigen_var{n_var} = Val_eigen{aa,bb,ii,jj}(n_var:q_vars:end);
    Vec_eigen_var{n_var} = Vec_eigen{aa,bb,ii,jj}(n_var:q_vars:end,:); % Take all columns for each perturbation row
end

% Spectrum at given alpha and Re
plot_Spectrum(y, delta, N, alpha(aa), N_target(ii), n_LST_Sweep(jj), Val_eigen{aa,bb,ii,jj},  Vec_eigen_var, name_file_load)

% Select Br and Beta
if length(n_LST_Sweep) > 1
    beta_target = 0;
    Br_target   = 0.00001;
    [~,ii] = min(abs(N_target-Br_target));
    [~,bb] = min(abs(beta-beta_target));
    plot_StabilityDiagram(alpha, N_target, n_LST_Sweep, squeeze(GG_max(:,bb,:,:)), name_file_load, bIsothermal)
end


%% Perturbation at given mode
realc = real(Val_eigen{aa,bb,ii,jj});
imagc = imag(Val_eigen{aa,bb,ii,jj});
Cond   = (abs(realc) <= 10) & (realc < 0.5);
Val_eigen_target   = max(imagc(Cond));
plot_Perturbation(y, delta, N, alpha(aa), N_target(ii),  n_LST_Sweep(jj), Val_eigen{aa,bb,ii,jj}, Val_eigen_target, Vec_eigen_var, name_file_load)
plot_Perturbation_inset(y, delta, N, alpha(aa), N_target(ii),  n_LST_Sweep(jj), Val_eigen{aa,ii,jj}, Val_eigen_target, Vec_eigen_var, name_file_load)


%% Growth rate diagram
% Select Br and Re
Re_target = 1000;
[~,jj] = min(abs(n_LST_Sweep-Re_target));
% [~,ii] = min(abs(N_target-Br_target));

for ii = 1:length(N_target)
    if ~isempty(find(N_target(ii) == N_plot))
        Br_target = N_target(ii);
        GR = squeeze(GG_max(:,:,ii,jj));
        GR_min = 0;
        GR_max = 200; % CO2 300-450-200 (sub-,trans-,super-)
        n_steps = 4; % CO2
%         GR_max = 200; % I-1-3
%         GR_max = 1000; % I-2
%         GR_max = 5000; % NI-1
%         GR_max = 750; % NI-2
%         GR_max = 900; % NI-3
%         GR_max = 3000; % NI-5
%         GR_max = 400; % NI-6
        Threshold = GR_max;
        GR(GR > Threshold) = Threshold;
        GR(GR < GR_min) = 0;
%         n_steps = 5; %Plot colorbar
        % Set GR = 0 for alpha = 0, beta = 0
        % if alpha(aa) == 0 && beta(bb) == 0
        %     GR(1,1) = 0;
        % end
        plot_GrowthRate(alpha, beta, Br_target, Re_target, GR, name_file_load, GR_min, GR_max,n_steps)
    end
end

%% Growth rate diagram at Re-alpha
% Select Br and Re
beta_target = 0;
[~,bb] = min(abs(beta-beta_target));
% [~,ii] = min(abs(N_target-Br_target));

for ii = 1:length(N_target)
    if ~isempty(find(N_target(ii) == N_plot))
        Br_target = N_target(ii);
        GR = squeeze(GG_max(:,bb,ii,:));
        
        % Choose max and min GR for plots
%         GR_max = 5000; GR_min = 200;  n_steps = 6; % I-1
%         GR_max = 5000; GR_min = 200; n_steps = 6; % I-2
%         GR_max = 80;  GR_min = 0;  n_steps = 5;% I-3
%         GR_max = 1000; GR_min = 200;  n_steps = 5;% NI-1
%           GR_max = 1000; GR_min = 200;  n_steps = 5;% NI-2
          GR_max = 1000; GR_min = 200;  n_steps = 5;% NI-3
%           GR_max = 80; GR_min = 0;  n_steps = 5; % NI-4
%           GR_max = 1000; GR_min = 200;  n_steps = 5; % NI-5
%           GR_max = 80; GR_min = 0;  n_steps = 5; % NI-6
        Threshold = GR_max;
        GR(GR > Threshold) = Threshold;
        GR(GR < GR_min) = 0;
%         n_steps = 7; %Plot colorbar
        % Set GR = 0 for alpha = 0, beta = 0
        % if alpha(aa) == 0 && beta(bb) == 0
        %     GR(1,1) = 0;
        % end
        plot_GrowthRate_Re_alpha(alpha, beta_target, Br_target, n_LST_Sweep, GR, squeeze(GG_unstable(:,bb,ii,:)), name_file_load, GR_min, GR_max,n_steps,bIsothermal)

        % Identify local max > modal stable or unstable
        disp(" Transient growth unstable or stable modes for Br = " + num2str(Br_target) + "...")
        for aa = 1:length(alpha)
            for jj = 1:length(n_LST_Sweep)
                if GG_unstable(aa,bb,ii,jj) == 1
                    disp(" Mode at Re = " + num2str(n_LST_Sweep(jj)) + " and alpha = " + alpha(aa) + " is UNSTABLE")
                    disp(" ")
                end
            end
        end
    end
end


%% Optimal perturbation
beta_target  = 2;
alpha_target = 0;
Re_target = 1000;
Br_target = 0.50;
[~,aa] = min(abs(alpha-alpha_target));
[~,bb] = min(abs(beta-beta_target));
[~,ii] = min(abs(N_target-Br_target));
[~,jj] = min(abs(n_LST_Sweep-Re_target));

pos_in  = 4; % Normalize w
pos_out = 1; % Normalize u post_out = 2;

% Optimal perturbation and response
plot_OptimalPerturbation(BF{ii}.y, delta, q_vars, q_in{aa,bb,ii,jj}, q_out{aa,bb,ii,jj}, pos_in, pos_out,name_file_load, N_target(ii), alpha(aa), beta(bb))

%% Optimal perturbation pattern
beta_target  = 2;
alpha_target = 0;
Re_target = 1000;
Br_target = 0.10;
[~,aa] = min(abs(alpha-alpha_target));
[~,bb] = min(abs(beta-beta_target));
[~,ii] = min(abs(N_target-Br_target));
[~,jj] = min(abs(n_LST_Sweep-Re_target));

pos_in  = 4; % Normalize w
pos_out = 2; % Normalize u

% Optimal perturbation pattern
for n_var = 1:q_vars
    Vec_eigen_var{n_var} = q_in{aa,bb,ii,jj}(n_var:q_vars:end); % Take all columns for each perturbation row
end

Norm = max(abs(Vec_eigen_var{pos_in}));
z      = linspace(0,1,200); % = Beta Z / (2pi)
% ph     = 0.20; % Isothermal
% ph     = 0.05; % Non-Isothermal NI-1
ph     = 0.45; % Non-Isothermal NI-5

v      = Vec_eigen_var{3}/Norm; 
w      = Vec_eigen_var{4}/Norm;

% Optimal response pattern
for n_var = 1:q_vars
    Vec_eigen_var{n_var} = q_out{aa,bb,ii,jj}(n_var:q_vars:end); % Take all columns for each perturbation row
end

Norm = max(abs(Vec_eigen_var{pos_out}));

rho = Vec_eigen_var{1}/Norm; 
u   = Vec_eigen_var{2}/Norm; 
T   = Vec_eigen_var{5}/Norm; 

% Plot Spanwise (alpa = 0)
plot_OptimalPerturbationPattern(z,y,v,w,rho,u,T,ph,name_file_load,N_target(ii),alpha(aa),beta(bb))



%% Optimal perturbation multiple sets
% Isothermal flow cases (I1-3)
% setup{1} =  'N2_Tbw_94.6425_Ttw_94.6425_Pb_5093700_PrEc_bulk_Isothermal_TG';
% setup{2} =  'N2_Tbw_119.8805_Ttw_119.8805_Pb_5093700_PrEc_bulk_Isothermal_TG';
% setup{3} =  'N2_Tbw_189.285_Ttw_189.285_Pb_5093700_PrEc_bulk_Isothermal_TG';
% setup_name = {'I-1','I-2','I-3'};
% name_save  = 'N2_Isothermal';

% % Non-Isothermal flow cases (NI1-3)
% setup{1} =  'N2_Tbw_94.6425_Ttw_189.285_Pb_5093700_PrEc_bulk_Non_Isothermal_TG';
% setup{2} =  'N2_Tbw_94.6425_Ttw_189.285_Pb_16979000_PrEc_bulk_Non_Isothermal_TG';
% setup{3} =  'N2_Tbw_113.571_Ttw_138.809_Pb_5093700_PrEc_bulk_Non_Isothermal_TG';
% setup_name = {'NI-1','NI-2','NI-3'};
% name_save  = 'N2_Non_Isothermal';

% % Non-Isothermal flow cases (N4-5)
setup{1} =  'N2_Tbw_94.6425_Ttw_189.285_Pb_6791600_PrEc_bulk_Non_Isothermal_TG_Ur_1';
setup{2} =  'N2_Tbw_94.6425_Ttw_189.285_Pb_101874_PrEc_bulk_Non_Isothermal_TG_Ur_1';
setup_name = {'NI-5','NI-6'};
name_save  = 'N2_Non_Isothermal_Ur_1';

% This should be based on current setup for input positions
beta_target  = 2;
alpha_target = 0;
Re_target = 1000;
% Br_target = 0.50; % Isothermal
Br_target = 0.10; % Non-isothermal

pos_in  = 4; % Normalize w
% pos_out = 2; % Normalize u isothermal
pos_out = 1; % Normalize rho Non-isothermal

plot_OptimalPerturbation_MultipleCases(setup, setup_name, delta, q_vars,pos_in, pos_out,alpha_target,beta_target,Re_target,Br_target,name_save);




%% Transient growth over time
beta_target  = 2.0;
alpha_target = 0.0;
Re_target    = 1000;

plot_TransientGrowth_Time(t_all, GG, alpha, alpha_target, beta, beta_target, n_LST_Sweep, Re_target, N_target, N_plot, name_file_load)

% Inset plot
beta_target2  = 1.0;
alpha_target2 = 1.0;
plot_TransientGrowth_Time(t_all, GG, alpha, alpha_target, beta, beta_target, n_LST_Sweep, Re_target, N_target, N_plot, name_file_load, alpha_target2, beta_target2)



function plot_Perturbation_inset(y, delta, N, alpha, N_target, n_LST_Sweep, Val_eigen, Val_eigen_target, Vec_eigen_var, name_file_load)

% Plot every i_plot
i_plot   = 2;
i_smooth = 3;

% Plot eigen vector with respect to y for each perturbation
[~,idx]            = min(abs(imag(Val_eigen) - Val_eigen_target));

% Normalize by u'
Norm = max(abs(Vec_eigen_var{2}(:,idx)));

f = figure;
hold on; grid on; box on;
plot(y(1:i_plot:end)/delta,smooth(abs(Vec_eigen_var{1}(1:i_plot:end,idx))/Norm,i_smooth),'LineWidth',2,'LineStyle',':','color',[0.4660 0.6740 0.1880]);
plot(y(1:i_plot:end)/delta,smooth(abs(Vec_eigen_var{2}(1:i_plot:end,idx))/Norm,i_smooth),'LineWidth',2,'LineStyle','-','color',[0 0.4470 0.7410])
plot(y(1:i_plot:end)/delta,smooth(abs(Vec_eigen_var{3}(1:i_plot:end,idx))/Norm,i_smooth),'LineWidth',2,'LineStyle','--','color',[0.8500 0.3250 0.0980])
plot(y(1:i_plot:end)/delta,smooth(abs(Vec_eigen_var{5}(1:i_plot:end,idx))/Norm,i_smooth),'LineWidth',2,'LineStyle','-.','color',[0.4940 0.1840 0.5560])
xlabel('${y/\delta}$','interpreter','latex')
ylabel('$|q^{\prime}|$','interpreter','latex')

legend('$|\rho^{\prime}|$','$|u^{\prime}|$','$|v^{\prime}|$', '$|T^{\prime}|$','interpreter','latex')
legend('Location','northeast','box','off','AutoUpdate','off')

set(gca,'linewidth',1.5)
set(gca,'fontsize',14)


% axes('position',[.2 .65 .4 .2]) % NI-5 Re = 10000
% axes('position',[.25 .4 .35 .3]) % NI-5 Re = 4000
axes('position',[.38 .56 .35 .3]) % NI-6 Re = 10000
% axes('position',[.4 .58 .35 .3]) % NI-6 Re = 4000

hold on; grid off; box on;
plot(y(1:i_plot:end)/delta,smooth(abs(Vec_eigen_var{1}(1:i_plot:end,idx))/Norm,i_smooth),'LineWidth',2,'LineStyle',':','color',[0.4660 0.6740 0.1880]);
plot(y(1:i_plot:end)/delta,smooth(abs(Vec_eigen_var{2}(1:i_plot:end,idx))/Norm,i_smooth),'LineWidth',2,'LineStyle','-','color',[0 0.4470 0.7410])
plot(y(1:i_plot:end)/delta,smooth(abs(Vec_eigen_var{3}(1:i_plot:end,idx))/Norm,i_smooth),'LineWidth',2,'LineStyle','--','color',[0.8500 0.3250 0.0980])
plot(y(1:i_plot:end)/delta,smooth(abs(Vec_eigen_var{5}(1:i_plot:end,idx))/Norm,i_smooth),'LineWidth',2,'LineStyle','-.','color',[0.4940 0.1840 0.5560])
% ylim([0 0.004]) % NI-5 Re = 10000
% ylim([0 0.02]) % NI-5 Re = 4000
ylim([0 0.007]) % NI-6 Re = 10000
% ylim([0 0.02]) % NI-6 Re = 4000

set(gca,'linewidth',1.5)
set(gca,'fontsize',10)


exportgraphics(f,strcat('Figures/',name_file_load,'Inset_alpha_',num2str(alpha), '_Re_',num2str(n_LST_Sweep),'.jpeg'),'Resolution',300)


% saveas(gca,strcat('Figures/', name_file_load,'_Spectrum_PrEc_',num2str(N_target),'_alpha_',num2str(alpha), '_Re_',num2str(n_LST_Sweep),'_N_', num2str(N)),'epsc')
% saveas(gca,strcat('Figures/', name_file_load,'_Spectrum_PrEc_',num2str(N_target),'_alpha_',num2str(alpha), '_Re_',num2str(n_LST_Sweep),'_N_', num2str(N)),'png')



end
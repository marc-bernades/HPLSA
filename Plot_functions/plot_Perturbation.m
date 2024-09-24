function plot_Perturbation(y, delta, N, alpha, N_target, n_LST_Sweep, Val_eigen, Val_eigen_target, Vec_eigen_var, name_file_load)



% Plot eigen vector with respect to y for each perturbation
[~,idx]            = min(abs(imag(Val_eigen) - Val_eigen_target));

% Normalize by u'
Norm = max(abs((Vec_eigen_var{2}(:,idx))));

f = figure;
hold on; grid on; box on;
plot(y/delta,abs(Vec_eigen_var{1}(:,idx))/Norm,'LineWidth',2,'LineStyle',':','color',[0.4660 0.6740 0.1880]);
plot(y/delta,abs(Vec_eigen_var{2}(:,idx))/Norm,'LineWidth',2,'LineStyle','-','color',[0 0.4470 0.7410])
plot(y/delta,abs(Vec_eigen_var{3}(:,idx))/Norm,'LineWidth',2,'LineStyle','--','color',[0.8500 0.3250 0.0980])
plot(y/delta,abs(Vec_eigen_var{5}(:,idx))/Norm,'LineWidth',2,'LineStyle','-.','color',[0.4940 0.1840 0.5560])
xlabel('${y/\delta}$','interpreter','latex')
ylabel('$|q^{\prime}|$','interpreter','latex')

legend('$|\rho^{\prime}|$','$|u^{\prime}|$','$|v^{\prime}|$', '$|T^{\prime}|$','interpreter','latex')
legend('Location','northeast','box','off')

set(gca,'linewidth',1.5)
set(gca,'fontsize',14)

exportgraphics(f,strcat('Figures/',name_file_load,'_alpha_',num2str(alpha), '_Re_',num2str(n_LST_Sweep),'.jpeg'),'Resolution',300)


% saveas(gca,strcat('Figures/', name_file_load,'_Spectrum_PrEc_',num2str(N_target),'_alpha_',num2str(alpha), '_Re_',num2str(n_LST_Sweep),'_N_', num2str(N)),'epsc')
% saveas(gca,strcat('Figures/', name_file_load,'_Spectrum_PrEc_',num2str(N_target),'_alpha_',num2str(alpha), '_Re_',num2str(n_LST_Sweep),'_N_', num2str(N)),'png')


% G_max_var     = cell(1,q_vars);
% for n_var = 1:q_vars
%     for aa = 1:length(alpha)
%         for ii = 1:length(N_target)
%             G_max_var{n_var}(aa,ii)       = max(G{aa,ii}(n_var:5:end));
%         end
%     end
% end

end
function plot_EnergyBudgetCompressibility(y,delta,q_u,q_v,q_P,Ma_g,test_label,save_label)


Color_map = {[0 0.4470 0.7410],[0.4660 0.6740 0.1880],[0.8500 0.3250 0.0980],...
    [0.4940 0.1840 0.5560], [0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840], ...
    [0.9290 0.6940 0.1250], [1, 0, 0], [0, 0, 1], [0, 1, 0], [1, 1, 0], [0, 0.5, 0.5]};


figure;
hold on; grid on; box on;

% Breakdown eigen vector for each variable
for ii = 1:length(Ma_g)

    plot(y/delta,real(Ma_g{ii}),'LineWidth',2,'LineStyle','-','color', Color_map{ii});
    plot(y/delta,real(q_u{ii}*5),'LineWidth',2,'LineStyle',':','color', Color_map{ii});
    plot(y/delta,real(q_v{ii}*5),'LineWidth',2,'LineStyle','-.','color',Color_map{ii});
    plot(y/delta,real(q_P{ii}*5),'LineWidth',2,'LineStyle','--','color',Color_map{ii});

end

xlabel('${y/\delta}$','interpreter','latex')
ylabel('${Ma}_g \thinspace , \thinspace \thinspace |q^{\prime}|$','interpreter','latex')
% ylim([-1 0.5]) % Isothermal
% ylim([-1.5 0.5]) % Non-Isothermal


set(gca,'linewidth',1.5)
set(gca,'fontsize',14)

% Legend perturbation vector and Br sweep
% legend('$|\rho^{\prime}|$','$|u^{\prime}|$','$|v^{\prime}|$', '$|T^{\prime}|$','interpreter','latex','fontsize',14)
% legend('Location','best','box','off')

h(1) = plot(NaN,'-','color','k','LineWidth',2);
h(2) = plot(NaN,':','color','k','LineWidth',2);
h(3) = plot(NaN,'-.','color','k','LineWidth',2);
h(4) = plot(NaN,'--','color','k','LineWidth',2);

for kk = 1:length(Ma_g)
    h(kk+4) = plot(NaN,'-','color',Color_map{kk},'LineWidth',2);
end
clegend = test_label;
clengend_final = {'${Ma}_g$','$u^{\prime} \times 5$','$v^{\prime} \times 5$', '$P^{\prime} \times 5$',clegend{:}};
legend(h,clengend_final,'interpreter','latex','fontsize',14,'Location','southwest',Box='off',NumColumns=3)


set(gca,'linewidth',1.5)
set(gca,'fontsize',14)

exportgraphics(gca,strcat('Figures/EnergyBudget/Budget_Compressibility_',save_label,'.jpeg'),'Resolution',300)


% G_max_var     = cell(1,q_vars);
% for n_var = 1:q_vars
%     for aa = 1:length(alpha)
%         for ii = 1:length(N_target)
%             G_max_var{n_var}(aa,ii)       = max(G{aa,ii}(n_var:5:end));
%         end
%     end
% end

end
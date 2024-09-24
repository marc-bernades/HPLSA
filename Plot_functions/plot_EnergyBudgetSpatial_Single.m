function plot_EnergyBudgetSpatial_Single(y,delta,K,P,T,V,save_label)

Color_map = {[0 0.4470 0.7410],[0.4660 0.6740 0.1880],[0.8500 0.3250 0.0980],...
    [0.4940 0.1840 0.5560], [0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840], ...
    [0.9290 0.6940 0.1250], [1, 0, 0], [0, 0, 1], [0, 1, 0], [1, 1, 0], [0, 0.5, 0.5]};



figure;
ax = axes;
hold on; grid on; box on;

% Plot LEFT
plot(y/delta,K{1},'LineWidth',2,'LineStyle','-','color', 'k');
plot(y/delta,P{1},'LineWidth',2,'LineStyle',':','color', Color_map{1});
plot(y/delta,T{1},'LineWidth',2,'LineStyle','-.','color',Color_map{2});
plot(y/delta,V{1},'LineWidth',2,'LineStyle','--','color',Color_map{3});
% ylim([-10 10]*10^-4)



clengend_final = {'$K$','$P$', '$T$','$V$'};
legend(clengend_final,'interpreter','latex','fontsize',14,'Location','north',Box='on')

xlabel('${y/\delta}$','interpreter','latex')
ylabel('$E$','interpreter','latex')
xlim([0 2])

set(gca,'linewidth',1.5)
set(gca,'fontsize',14)



exportgraphics(gca,strcat('Figures/EnergyBudget/Budget_',save_label,'.jpeg'),'Resolution',300)


end
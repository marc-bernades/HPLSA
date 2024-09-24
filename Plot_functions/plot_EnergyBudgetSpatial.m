function plot_EnergyBudgetSpatial(y,delta,K,P,T,V,test_label,save_label)

Color_map = {[0 0.4470 0.7410],[0.4660 0.6740 0.1880],[0.8500 0.3250 0.0980],...
    [0.4940 0.1840 0.5560], [0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840], ...
    [0.9290 0.6940 0.1250], [1, 0, 0], [0, 0, 1], [0, 1, 0], [1, 1, 0], [0, 0.5, 0.5]};



figure;
ax = axes;
hold on; grid on; box on;

% Plot LEFT
yyaxis left
idx_left = floor(length(K{1})/2);
plot(y(1:idx_left)/delta,K{1}(1:idx_left),'LineWidth',2,'LineStyle','-','color', 'k');
plot(y(1:idx_left)/delta,P{1}(1:idx_left),'LineWidth',2,'LineStyle',':','color', Color_map{1});
plot(y(1:idx_left)/delta,T{1}(1:idx_left),'LineWidth',2,'LineStyle','-.','color',Color_map{2});
plot(y(1:idx_left)/delta,V{1}(1:idx_left),'LineWidth',2,'LineStyle','--','color',Color_map{3});
% ylim([-10 10]*10^-4)
ylabel('$E$','interpreter','latex')


% Plot RIGHT
yyaxis right
idx_right = floor(length(K{2})/2) + 1;
plot(y(idx_right:end)/delta,K{2}(idx_right:end),'LineWidth',2,'LineStyle','-','color', 'k');
plot(y(idx_right:end)/delta,P{2}(idx_right:end),'LineWidth',2,'LineStyle',':','color', Color_map{1});
plot(y(idx_right:end)/delta,T{2}(idx_right:end),'LineWidth',2,'LineStyle','-.','color',Color_map{2});
plot(y(idx_right:end)/delta,V{2}(idx_right:end),'LineWidth',2,'LineStyle','--','color',Color_map{3});
ylim('auto')


% Middle line
plot([1 1],ylim,'--','color',[0.5 0.5 0.5],'LineWidth',2);


clengend_final = {'$K$','$P$', '$T$','$V$'};
legend(clengend_final,'interpreter','latex','fontsize',14,'Location','north',Box='on')

xlabel('${y/\delta}$','interpreter','latex')
xlim([0 2])

set(gca,'linewidth',1.5)
set(gca,'fontsize',14)

yyaxis right; ylimr = get(gca,'Ylim');ratio = ylimr(1)/ylimr(2);
yyaxis left; yliml = get(gca,'Ylim');
if yliml(2)*ratio<yliml(1)
    set(gca,'Ylim',[yliml(2)*ratio yliml(2)])
else
    set(gca,'Ylim',[yliml(1) yliml(1)/ratio])
end

% Exponent of left axis
ax.YAxis(1).Exponent = -2;
% Set left axis legend at top for I2 Br comparison
ax.YTick(end+1)      = yliml(1)/ratio; 
% Set left axis legend at top for I-1 vs. I-3
% Should be equivalent but they are not
ax.YTick(end+1)      = ax.YTick(end) + ax.YTick(end) - ax.YTick(end-1);

ax.YAxis(1).Color = [0 0 0];
ax.YAxis(2).Color = [0 0 0];

% Set labels each side
limits = ylim;
text(0.25,limits(2)*0.8,test_label{1},'interpreter','latex','FontSize',14)
text(1.25,limits(2)*0.8,test_label{2},'interpreter','latex','FontSize',14)


exportgraphics(gca,strcat('Figures/EnergyBudget/Budget_',save_label,'.jpeg'),'Resolution',300)


end

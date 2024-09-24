function plot_GrowthRate_Re_alpha(alpha, beta_target, N_target, n_LST_Sweep, G_max, GG_unstable, name_file_load, GR_min, GR_max, n_steps, bIsothermal, varargin)

% Contourplot interpolation
% f = figure; hold on; grid on; box on;

alpha_q          = linspace(min(alpha),max(alpha),500);
Re_q             = linspace(min(n_LST_Sweep),max(n_LST_Sweep),500);
[Re,alpha]       = meshgrid(n_LST_Sweep,alpha);
[Re_q,alpha_q]   = meshgrid(Re_q,alpha_q);
G_max_q          = griddata(Re,alpha,G_max,Re_q,alpha_q,'cubic'); %G_max;
G_max_q(G_max_q < 0) = 0;
GG_unstable_q    = griddata(Re,alpha,GG_unstable,Re_q,alpha_q,'linear'); %GG_unstable;

% Contour infinity energy growth
% [C,h1]=contourf(Re,alpha,GG_unstable,[1.0 0.0],'--','linewidth',3,'color','k'); % Plot max GMAX single level
close;
f = figure; hold on; grid on; box on;

% Contourplot
% [c,h]=contourf(alpha_q,beta_q,G_max_q,[linspace(min(min(G_max_q)),round(max(max(G_max_q)),1),500)],'HandleVisibility','off');
[c,h]=contour(Re,alpha,G_max,[linspace(GR_min,GR_max,10)],'HandleVisibility','off','LineWidth',2);
% set(h, 'edgecolor','none'); 
% set(gca,'XScale', 'log')
colorbar;
colormap('jet')
c = colorbar('Ticks',[linspace(GR_min,GR_max,n_steps)]);
c.Ruler.Exponent        = floor(log10(GR_max));                   % set the desired exponent
c.Ruler.TickLabelFormat = '%0.1f';       % fix up ugly default %g formatting
clim([GR_min GR_max])
cbh = findall(f, 'Type', 'ColorBar');
cTH = get(cbh,'Title');
set(cTH,'String',['$','G_{max}','$'],'Interpreter','latex','fontsize',14);
pbaspect([1 1.5 1])

% Plot unstable envelope
% x = C(1,2:end);
% y = C(2,2:end);
% idx_x = (x <=  Re(end,end)    & x>=  Re(1,1));
% idx_y = (y <=  alpha(end,end) & y>=  alpha(1,1));
% idx = (idx_x & idx_y);
% plot(x(idx),y(idx),'o','Color',[0.5 0.5 0.5])
% fill(x(idx),y(idx),[0.5 0.5 0.5])

% Alternative
Cond = GG_unstable == 1;
scatter(Re(Cond),alpha(Cond),50,'o','MarkerEdgeColor','k','MarkerFaceColor',[0.5 0.5 0.5],'LineWidth',1)

% data = C';
% n = size(data,1);
% data_interp = interp1(1:n,data,linspace(1,n,11*(n-1)+1),'makima');
% plot(data_interp(2:end,1),data_interp(2:end,2),'x','Color','y','MarkerSize',3)
% fill(data_interp(2:end,1),data_interp(2:end,2),[0.5 0.5 0.5])

set(gca,'fontsize',12)
set(gca,'linewidth',1.5)
xlabel('${Re}$','interpreter','latex')
ylabel('${\alpha}$','interpreter','latex')
xlim([min(n_LST_Sweep) max(n_LST_Sweep)])
if bIsothermal == 1
    ylim([0.6 1.6]); yticks([0.6 0.8 1.0 1.2 1.4 1.6]);
else
    ylim([0.4 1.2]); yticks([0.4 0.6 0.8 1.0 1.2]);
end
set(cTH,'Units','normalized','Position',[0.0 1.010 0])

xlim([min(n_LST_Sweep) max(n_LST_Sweep)]);
% xticklabels({'2','4','6','8','10'})
ax = gca; ax.XAxis.Exponent = 4;

set(gca,'linewidth',1.5)
set(gca,'fontsize',16)

exportgraphics(gca,strcat('Figures/',name_file_load,'_GrowthRate_Re_Alpha_','Br_',num2str(N_target),'_Beta_',num2str(beta_target),'.jpeg'),'Resolution',300)






end
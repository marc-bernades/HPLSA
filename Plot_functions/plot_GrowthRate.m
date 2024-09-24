function plot_GrowthRate(alpha, beta, N_target, n_LST_Sweep, G_max, name_file_load, GR_min, GR_max, n_steps, varargin)

% Contourplot interpolation
f = figure; hold on; grid on; box on;

alpha_q          = linspace(min(alpha),max(alpha),500);
beta_q           = linspace(min(beta),max(beta),500);
[beta,alpha]     = meshgrid(beta,alpha);
[beta_q,alpha_q] = meshgrid(beta_q,alpha_q);
G_max_q          = griddata(alpha,beta,G_max,alpha_q,beta_q,'cubic'); %G_max;
G_max_q(G_max_q<0) = 0;

% Contourplot
% [c,h]=contourf(alpha_q,beta_q,G_max_q,[linspace(min(min(G_max_q)),round(max(max(G_max_q)),1),500)],'HandleVisibility','off');
[c,h]=contourf(alpha_q,beta_q,G_max_q,[linspace(GR_min,GR_max,25)],'HandleVisibility','off');
set(h, 'edgecolor','none'); 
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

set(gca,'fontsize',12)
set(gca,'linewidth',1.5)
xlabel('${\alpha}$','interpreter','latex')
ylabel('${\beta}$','interpreter','latex')
ylim([0 5]); yticks([0.0 1.0 2.0 3.0 4.0 5.0]);
set(cTH,'Units','normalized','Position',[0.0 1.010 0])

if size(alpha,1)>1
    xlim([min(alpha(:,1)) max(alpha(:,1))]); %xticklabels({'2','4','6','8','10'})
else
    xlim([min(alpha) max(alpha)]); %xticklabels({'2','4','6','8','10'})
end

set(gca,'linewidth',1.5)
set(gca,'fontsize',16)

exportgraphics(gca,strcat('Figures/',name_file_load,'_GrowthRate_','Br_',num2str(N_target),'_Re_',num2str(n_LST_Sweep),'.jpeg'),'Resolution',300)






end
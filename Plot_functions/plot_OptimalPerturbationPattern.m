function plot_OptimalPerturbationPattern(z,y,v,w,rho,u,T,ph,name_file_load,Br,alpha,beta)

% [Z,Y] = meshgrid(z,y);
% figure
% contourf(Z,Y,real(v_convert)); colorbar
% figure
% contourf(Z,Y,real(w_convert)); colorbar

%% INPUT
% Convert to x,y,z,t space
z_convert  = z*2*pi/beta;   
v_convert  = (v.*exp(sqrt(-1).*beta.*(z_convert + ph)));
w_convert  = (w.*exp(sqrt(-1).*beta.*(z_convert + ph)));
vw_convert = sqrt((v_convert).^2 + real(w_convert).^2);
zz_plot    = z;
i_plot     = 10;

% quiver(zz_plot(1:i_plot:end),y(1:i_plot:end),real(w_convert(1:i_plot:end,1:i_plot:end)),real(v_convert(1:i_plot:end,1:i_plot:end)),'Color','k')
% Color by its intensity
figure
ax = axes('Color', [1 1 1]); 
hold(ax,'on')
box(ax,'on')
% Define the colormap for the quiver arrows.
% cmap can have any number of rows.
cmap        = jet(255); 
ax.Colormap = cmap; 
% Assign colors based on magnitude of vectors
vectorMagnitude = abs(vw_convert); 
% Scale magnitudes to rows of colormap
vecMagNorm = (vectorMagnitude-min(vectorMagnitude))./range(vectorMagnitude);
vecColorIdx = round(vecMagNorm * (size(cmap,1)-1)) + 1; 
% Plot the quiver data
for idx_z = 1:i_plot:numel(zz_plot)
    for idx_y = 1:i_plot:numel(y)
    quiver(ax, zz_plot(idx_z),y(idx_y),real(w_convert(idx_y,idx_z)),real(v_convert(idx_y,idx_z)), ...
        .2,'Color', cmap(vecColorIdx(idx_y,idx_z),:), 'LineWidth',1.0)
    end
end

% Set properties for the main axes
% axis equal
xlim(ax, [0 1])
ylim(ax, [0 2])
xlabel('${\beta z \thinspace / \thinspace (\delta \thinspace 2\pi)}$','interpreter','latex','fontsize',14)
ylabel('${y/\delta}$','interpreter','latex','fontsize',14)

% Add colorbar
cb = colorbar(ax); 
% Set colorbar range
caxis(ax, [floor(vectorMagnitude(1)), ceil(vectorMagnitude(2))])
cb.Ticks = linspace(floor(vectorMagnitude(1)),ceil(vectorMagnitude(2)),5);
% Label the colorbars
% ylabel(cb,'Vector magnitude')

set(gca,'linewidth',1.5)
set(gca,'fontsize',16)


cTH = get(cb,'Title');
set(cTH,'String',['$','\vert \textbf{u}^\prime \vert','$'],'Interpreter','latex','fontsize',14);
pbaspect([1 1.0 1])

exportgraphics(gca,strcat('Figures/',name_file_load,'_OptimalPerturbationPattern_IN_Br_',num2str(Br),'_Alpha_',num2str(alpha),'_Beta_',num2str(beta),'.jpeg'),'Resolution',300)


%% OUTPUT
% Convert to x,y,z,t space
rho_convert = real(rho.*exp(sqrt(-1).*beta.*(z_convert + ph)));
u_convert   = real(u.*exp(sqrt(-1).*beta.*(z_convert + ph)));
T_convert   = real(T.*exp(sqrt(-1).*beta.*(z_convert + ph)));

[Z,Y] = meshgrid(zz_plot,y);

% DENSITY
f = figure
[c,h]=contourf(Z,Y,rho_convert,[linspace(min(min(rho_convert)),max(max(rho_convert)),250)],'HandleVisibility','off');
set(h, 'edgecolor','none'); 
% set(gca,'XScale', 'log')
colorbar;
colormap('jet')
% c = colorbar('Ticks',[linspace(min(min(rho_convert)),max(max(rho_convert)),5)]);
% c.Ruler.Exponent        = floor(log10(GR_max));                   % set the desired exponent
% c.Ruler.TickLabelFormat = '%0.1f';       % fix up ugly default %g formatting

% clim([min(min(rho_convert)) max(max(rho_convert))])
cbh = findall(f, 'Type', 'ColorBar');
cTH = get(cbh,'Title');
set(cTH,'String',['$','\rho^\prime','$'],'Interpreter','latex','fontsize',14);
pbaspect([1 1.0 1])

set(gca,'fontsize',12)
set(gca,'linewidth',1.5)
xlabel('${\beta z \thinspace / \thinspace (\delta \thinspace 2\pi)}$','interpreter','latex','fontsize',14)
xlim([0 1]); xticks([0.0 0.25 0.5 0.75 1.0]);
ylabel('${y/\delta}$','interpreter','latex','fontsize',14)
ylim([0 2]); yticks([0.0 0.5 1.0 1.5 2.0]);
% set(cTH,'Units','normalized','Position',[0.0 1.010 0])

set(gca,'linewidth',1.5)
set(gca,'fontsize',16)

exportgraphics(gca,strcat('Figures/',name_file_load,'_OptimalPerturbationPattern_OUT_Density_Br_',num2str(Br),'_Alpha_',num2str(alpha),'_Beta_',num2str(beta),'.jpeg'),'Resolution',300)


% VELOCITY
f = figure
[c,h]=contourf(Z,Y,u_convert,[linspace(min(min(u_convert)),max(max(u_convert)),250)],'HandleVisibility','off');
set(h, 'edgecolor','none'); 
colorbar;
colormap('jet')

cbh = findall(f, 'Type', 'ColorBar');
cTH = get(cbh,'Title');
set(cTH,'String',['$','u^\prime','$'],'Interpreter','latex','fontsize',14);
pbaspect([1 1.0 1])
clim([-1 1]);
cbh.Ticks = [-1:0.25:1]

set(gca,'fontsize',12)
set(gca,'linewidth',1.5)
xlabel('${\beta z \thinspace / \thinspace (\delta \thinspace 2\pi)}$','interpreter','latex','fontsize',14)
xlim([0 1]); xticks([0.0 0.25 0.5 0.75 1.0]);
ylabel('${y^\star}$','interpreter','latex','fontsize',14)
ylim([0 2]); yticks([0.0 0.5 1.0 1.5 2.0]);
% set(cTH,'Units','normalized','Position',[0.0 1.010 0])

exportgraphics(gca,strcat('Figures/',name_file_load,'_OptimalPerturbationPattern_OUT_Velocity_X_Br_',num2str(Br),'_Alpha_',num2str(alpha),'_Beta_',num2str(beta),'.jpeg'),'Resolution',300)


% TEMPERATURE
f = figure
[c,h]=contourf(Z,Y,T_convert,[linspace(min(min(T_convert)),max(max(T_convert)),250)],'HandleVisibility','off');
set(h, 'edgecolor','none'); 
colorbar;
colormap('jet')

cbh = findall(f, 'Type', 'ColorBar');
cTH = get(cbh,'Title');
set(cTH,'String',['$','T^\prime','$'],'Interpreter','latex','fontsize',14);
pbaspect([1 1.0 1])

set(gca,'fontsize',12)
set(gca,'linewidth',1.5)
xlabel('${\beta z \thinspace / \thinspace (\delta \thinspace 2\pi)}$','interpreter','latex','fontsize',14)
xlim([0 1]); xticks([0.0 0.25 0.5 0.75 1.0]);
ylabel('${y^\star}$','interpreter','latex','fontsize',14)
ylim([0 2]); yticks([0.0 0.5 1.0 1.5 2.0]);
% set(cTH,'Units','normalized','Position',[0.0 1.010 0])

exportgraphics(gca,strcat('Figures/',name_file_load,'_OptimalPerturbationPattern_OUT_Temperature_Br_',num2str(Br),'_Alpha_',num2str(alpha),'_Beta_',num2str(beta),'.jpeg'),'Resolution',300)

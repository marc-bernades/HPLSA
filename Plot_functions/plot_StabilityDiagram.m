function plot_StabilityDiagram(alpha, N_target, N_full, n_LST_Sweep, G_max, name_file_load, bIsothermal, varargin)


Color_map = {[0 0.4470 0.7410],[0.4660 0.6740 0.1880],[0.8500 0.3250 0.0980],...
    [0.4940 0.1840 0.5560], [0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840], ...
    [0.9290 0.6940 0.1250], [1, 0, 0], [0, 0, 1], [0, 1, 0], [1, 1, 0], [0, 0.5, 0.5]};

% Incompressible critical Re
Re_c = 5850;

% Contourplot interpolation
figure; hold on; grid on; box on;
for i_N_target = 1:length(N_target)

    [~,ii] = min(abs(N_target(i_N_target)-N_full));

    %% For Growth-rate non-isothermal contour
    % Make interpolation
%     n_LST_Sweep_q    = linspace(min(n_LST_Sweep),max(n_LST_Sweep),500);
%     alpha_q          = linspace(min(alpha),max(alpha),500);
%     [n_LST_Sweep_q,alpha_q] = meshgrid(n_LST_Sweep_q,alpha_q);

%     G_max_q = griddata(n_LST_Sweep,alpha,squeeze(G_max(:,ii,:)),n_LST_Sweep_q,alpha_q,'cubic');
%     G_max_q = interp2(n_LST_Sweep,alpha,squeeze(G_max(:,ii,:)),n_LST_Sweep_q,alpha_q,'cubic');

    n_LST_Sweep_q    = n_LST_Sweep; %linspace(min(n_LST_Sweep),max(n_LST_Sweep),1000);
    alpha_q          = alpha; %linspace(min(alpha),max(alpha),1000);
    G_max_q          = squeeze(G_max(:,ii,:));

    % Isothermal case plot iosthermal limit separeatly
    if bIsothermal
        if ii == 1 && length(N_target) > 1
            [C,h]=contour(n_LST_Sweep_q,alpha_q,G_max_q,[0,0],'--','linewidth',3,'color','k'); % Plot max GMAX single level
            close;
            figure; hold on; grid on; box on;
            %clabel(C,h)
        else
            [c,h1]=contour(n_LST_Sweep_q,alpha_q,G_max_q,[0,0],'linewidth',3,'color',Color_map{ii-1});
            Re_disp = c(1,2:end);
            disp("Re_c = " + num2str(min(Re_disp(Re_disp~= 0))) + " for Br = " + num2str(N_target(ii)))
        end
    else
        % Non-isothermal
        [c1,h1]=contour(n_LST_Sweep_q,alpha_q,G_max_q,[0,0],'linewidth',3,'color',Color_map{i_N_target});
        Re_disp = c1(1,2:end);
        disp("Re_c = " + num2str(min(Re_disp(Re_disp~= 0))) + " for Br = " + num2str(N_target(i_N_target)))

        %% For Growth-rate non-isothermal contour
%         close;
%         f = figure; hold on; grid on; box on;
%         [c,h1]=contourf(n_LST_Sweep_q,alpha_q,G_max_q,[linspace(min(min(G_max_q)),round(max(max(G_max_q)),1),500)],'HandleVisibility','off');
%         set(h1, 'edgecolor','none');
%         % set(gca,'XScale', 'log')
%         colorbar;
%         colormap('jet')
%         % colorbar('Ticks',[0.75 1 1.25 1.5 1.75])
%         % clim([0 8])  NI-5
% %         clim([-0.1 0]); % NI-6
%         cbh = findall(f, 'Type', 'ColorBar');
%         cTH = get(cbh,'Title');
%         set(cTH,'String',['$','G_{max}','$'],'Interpreter','latex','fontsize',14);

        
    end


    clegend{ii} = strcat('$', 'Br', ' = ',num2str(N_target(i_N_target)), '$');
end

plot([Re_c Re_c],[min(alpha)-0.2,max(alpha)+0.2],'k--','linewidth',2)
% Isothermal plot isothermal-limit
if bIsothermal && length(N_target) > 1
    plot(C(1,2:end),C(2,2:end),'k-.','linewidth',3)
end

%% For Growth-rate non-isothermal contour
% Plot neutral curve in red
% plot(c1(1,2:end),c1(2,2:end),'r-.','linewidth',3)


set(gca,'fontsize',16)
set(gca,'linewidth',1.5)
xlabel('$Re \thinspace (\times 10^3)$','interpreter','latex')
ylabel('${\alpha}$','interpreter','latex')
if bIsothermal

    if length(N_target) > 1
        % Don't put legend Br = 0
        legend(clegend{2:end},'interpreter','latex', 'location','best','Box','on','NumColumns',2,'fontsize',14)
    else
        legend(clegend{1:end},'interpreter','latex', 'location','best','Box','on','NumColumns',2,'fontsize',14)
    end
    ylim([0.6 1.6]); yticks([0.6 0.8 1.0 1.2 1.4 1.6]);

else
    % Non-isothermal
%     legend(clegend{1:end},'interpreter','latex', 'location','best','Box','on','NumColumns',2,'fontsize',16)
    ylim([0.4 1.2]); yticks([0.4 0.6 0.8 1.0 1.2]);
end

xticklabels({'2','4','6','8','10'})
% saveas(gca,strcat('Figures/',name_file_load,'_StabilityDiagram'),'epsc')
% saveas(gca,strcat('Figures/',name_file_load,'_StabilityDiagram'),'jpeg')
exportgraphics(gca,strcat('Figures/',name_file_load,'_StabilityDiagram','.jpeg'),'Resolution',300)




end
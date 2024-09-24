function plot_Baseflow_LST(y,delta,BF,N_target,T_c,name_file_load, varargin)

Color_map = {[0 0.4470 0.7410],[0.4660 0.6740 0.1880],[0.8500 0.3250 0.0980],...
    [0.4940 0.1840 0.5560], [0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840], ...
    [0.9290 0.6940 0.1250], [1, 0, 0], [0, 0, 1], [0, 1, 0], [1, 1, 0], [0, 0.5, 0.5]};

if isempty(varargin)
    % Velocity
    figure; hold on; box on; grid on;
    for n_BF = 1:length(BF)

%         plot(y/delta,BF{n_BF}.u/BF{n_BF}.norm.u,'linewidth',2, 'color',Color_map{n_BF});
        plot(y/delta,BF{n_BF}.u,'linewidth',2, 'color',Color_map{n_BF});

        clabel{n_BF} = strcat('$', ' PrEc = ',num2str(N_target(n_BF)), '$');

    end


    xlabel('${{y}/{\delta}}$','interpreter','latex')
    ylabel('${{u}/{u_c}}$','interpreter','latex')
    % pbaspect([1.8 1 1])
    legend(clabel,'interpreter','latex', 'location','south')

    set(gca,'linewidth',1.5)
    set(gca,'fontsize',14)

    % set(gca, 'OuterPosition', [0,0,1,1])


%     saveas(gca,strcat('Figures/', name_file_load,'_Baseflow_u'),'epsc')
%     saveas(gca,strcat('Figures/', name_file_load,'_Baseflow_u'),'png')


    % Temperature
    figure; hold on; box on; grid on;
    for n_BF = 1:length(BF)

%         plot(y/delta,BF{n_BF}.T/BF{n_BF}.norm.T, 'linewidth',2,'color',Color_map{n_BF});
        plot(y/delta,BF{n_BF}.T, 'linewidth',2,'color',Color_map{n_BF});

    end


    xlabel('${{y}/{\delta}}$','interpreter','latex')
    ylabel('${{T}/{T_w}}$','interpreter','latex')
    legend(clabel,'interpreter','latex', 'location','best')

    set(gca,'linewidth',1.5)
    set(gca,'fontsize',14)


%     saveas(gca,strcat('Figures/', name_file_load,'_Baseflow_T'),'epsc')
%     saveas(gca,strcat('Figures/', name_file_load,'_Baseflow_T'),'png')


else
    % Only plot the selected PrEc
    ii = varargin{1};

    % Velocity
    figure; hold on; box on; grid on;
    plot(y/delta,BF{ii}.u/BF{ii}.norm.u,'linewidth',2, 'color',Color_map{1});
    clabel = strcat('$', 'PrEc = ',num2str(N_target(ii)), '$');

    xlabel('${{y}/{\delta}}$','interpreter','latex')
    ylabel('${{u}/{u_c}}$','interpreter','latex')
    legend(clabel,'interpreter','latex', 'location','south')
    set(gca,'linewidth',1.5)
    set(gca,'fontsize',14)


    saveas(gca,strcat('Figures/', name_file_load,'_Baseflow_u_'),'epsc')
    saveas(gca,strcat('Figures/', name_file_load,'_Baseflow_u_'),'png')


    % Temperature
    figure; hold on; box on; grid on;

    plot(y/delta,BF{ii}.T/BF{n_BF}.norm.T, 'linewidth',2,'color',Color_map{1});

    xlabel('${{y}/{\delta}}$','interpreter','latex')
    ylabel('${{T}/{T_w}}$','interpreter','latex')
    legend(clabel,'interpreter','latex', 'location','best')

    set(gca,'linewidth',1.5)
    set(gca,'fontsize',14)



    saveas(gca,strcat('Figures/', name_file_load,'_Baseflow_T_'),'epsc')
    saveas(gca,strcat('Figures/', name_file_load,'_Baseflow_T_'),'png')


end


end
function plot_OptimalPerturbation_MultipleCases(setup, setup_name, delta, q_vars,pos_in, pos_out,alpha_target,beta_target,Re_target,Br_target, name_save)

Color_map = {[0 0.4470 0.7410],[0.4660 0.6740 0.1880],[0.8500 0.3250 0.0980],...
    [0.4940 0.1840 0.5560], [0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840], ...
    [0.9290 0.6940 0.1250], [1, 0, 0], [0, 0, 1], [0, 1, 0], [1, 1, 0], [0, 0.5, 0.5]};

% Smooth and figure sample show
i_smooth = 3;
i_sample = 2;

%% Unpack structure and load variables
for idx = 1:length(setup)
    Results{idx} = load(strcat('Results/Non_Modal/', setup{idx}, '.mat'));
end

%% Q_IN
f = figure
hold on; grid on; box on;
for idx = 1:length(setup)
    [~,aa] = min(abs(Results{idx}.Data.alpha-alpha_target));
    [~,bb] = min(abs(Results{idx}.Data.beta-beta_target));
    [~,ii] = min(abs(Results{idx}.Data.N_target-Br_target));
    [~,jj] = min(abs(Results{idx}.Data.n_LST_Sweep-Re_target));

    % Plot q_in
    y         = Results{idx}.Data.y;
    q_in_plot = Results{idx}.Data.q_in{aa,bb,ii,jj};
    for n_var = 1:q_vars
        Vec_eigen_var{n_var} = q_in_plot(n_var:q_vars:end); % Take all columns for each perturbation row
    end

    Norm = max(abs(Vec_eigen_var{pos_in}));

%     plot(y/delta,abs(Vec_eigen_var{1})/Norm,'LineWidth',2,'LineStyle','-','color',[0 0.4470 0.7410]);
%     plot(y/delta,abs(Vec_eigen_var{2})/Norm,'LineWidth',2,'LineStyle',':','color',[0.4660 0.6740 0.1880])
    plot( y(1:i_sample:end)/delta,smooth(abs(Vec_eigen_var{3}(1:i_sample:end))/Norm,i_smooth),'LineWidth',2,'LineStyle','-','color',Color_map{idx});
    plot( y(1:i_sample:end)/delta,smooth(abs(Vec_eigen_var{4}(1:i_sample:end))/Norm,i_smooth),'LineWidth',2,'LineStyle',':','color',Color_map{idx});
%     plot(y/delta,abs(Vec_eigen_var{5})/Norm,'LineWidth',2,'LineStyle','--','color',[0.8500 0.3250 0.0980])

end

xlabel('${y/\delta}$','interpreter','latex')
ylabel('$|q^{\prime}|$','interpreter','latex')
pbaspect([1 1.2 1])


%     legend('$|\rho^{\prime}|$','$|u^{\prime}|$','$|v^{\prime}|$', '$|w^{\prime}|$','$|T^{\prime}|$','interpreter','latex')
% legend('Location','north','box','off')
%     ylim([0 1])

set(gca,'linewidth',1.5)
set(gca,'fontsize',16)
% f.Children.YTick = [f.Children.YTick 1.4]; %Update label top ylim

% Legend perturbation vector and Br sweep
% legend('$|\rho^{\prime}|$','$|u^{\prime}|$','$|v^{\prime}|$', '$|T^{\prime}|$','interpreter','latex','fontsize',16)
% legend('Location','best','box','off')

h(1) = plot( NaN,'-','color','k','LineWidth',2);
h(2) = plot( NaN,':','color','k','LineWidth',2);
% h(3) = plot(NaN,'-.','color','k','LineWidth',2);
% h(4) = plot(NaN,'--','color','k','LineWidth',2);
idx_add = length(h);

for kk = 1:length(setup)
    h(kk+idx_add) = plot( NaN,'-','color',Color_map{kk},'LineWidth',2);
    clegend{kk}   = strcat('$', setup_name{kk}, '$');
end
% clengend_final = {'$|\rho^{\prime}|$','$|u^{\prime}|$','$|v^{\prime}|$', '$|T^{\prime}|$',clegend{:}};
clengend_final = {'$|v^{\prime}|$', '$|w^{\prime}|$',clegend{:}};

legend(h,clengend_final,'interpreter','latex','fontsize',16,'Location','north',Box='off')


exportgraphics(gca,strcat('Figures/','Optimal_Input_',name_save,'_alpha_',num2str(alpha_target),'_beta_',num2str(beta_target),'_Br_',num2str(Br_target),'_Re_',num2str(Re_target),'.jpeg'),'Resolution',300)


clear h
clear clegend
clear clengend_final

%% Q_OUT
f = figure
hold on; grid on; box on;
for idx = 1:length(setup)
    [~,aa] = min(abs(Results{idx}.Data.alpha-alpha_target));
    [~,bb] = min(abs(Results{idx}.Data.beta-beta_target));
    [~,ii] = min(abs(Results{idx}.Data.N_target-Br_target));
    [~,jj] = min(abs(Results{idx}.Data.n_LST_Sweep-Re_target));

    % Plot q_out
    y          = Results{idx}.Data.y;
    q_out_plot = Results{idx}.Data.q_out{aa,bb,ii,jj};
    for n_var = 1:q_vars
        Vec_eigen_var{n_var} = q_out_plot(n_var:q_vars:end); % Take all columns for each perturbation row
    end

    Norm = max(abs(Vec_eigen_var{pos_out}));

    plot(y/delta,abs(Vec_eigen_var{1})/Norm,'LineWidth',2,'LineStyle','-','color',Color_map{idx});
    plot(y/delta,abs(Vec_eigen_var{2})/Norm,'LineWidth',2,'LineStyle',':','color',Color_map{idx});
%     plot( y/delta,abs(Vec_eigen_var{3})/Norm,'LineWidth',2,'LineStyle','-','color',Color_map{idx});
%     plot( y/delta,abs(Vec_eigen_var{4})/Norm,'LineWidth',2,'LineStyle',':','color',Color_map{idx});
    plot(y/delta,abs(Vec_eigen_var{5})/Norm,'LineWidth',2,'LineStyle','--','color',Color_map{idx})

end

xlabel('${y/\delta}$','interpreter','latex')
ylabel('$|q^{\prime}|$','interpreter','latex')
pbaspect([1 1.2 1])


%     legend('$|\rho^{\prime}|$','$|u^{\prime}|$','$|v^{\prime}|$', '$|w^{\prime}|$','$|T^{\prime}|$','interpreter','latex')
% legend('Location','north','box','off')


set(gca,'linewidth',1.5)
set(gca,'fontsize',16)

% ylim([0 1.6])
% f.Children.YTick = [0 0.4 0.8 1.2 1.6]; %Update label top ylim


% Legend perturbation vector and Br sweep
% legend('$|\rho^{\prime}|$','$|u^{\prime}|$','$|v^{\prime}|$', '$|T^{\prime}|$','interpreter','latex','fontsize',16)
% legend('Location','best','box','off')

h(1) = plot( NaN,'-','color','k','LineWidth',2);
h(2) = plot( NaN,':','color','k','LineWidth',2);
h(3) = plot( NaN,'--','color','k','LineWidth',2);
% h(4) = plot(NaN,'--','color','k','LineWidth',2);
idx_add = length(h);

for kk = 1:length(setup)
    h(kk+idx_add) = plot( NaN,'-','color',Color_map{kk},'LineWidth',2);
    clegend{kk}   = strcat('$', setup_name{kk}, '$');
end

clengend_final = {'$|\rho^{\prime}|$', '$|u^{\prime}|$','$|T^{\prime}|$',clegend{:}};
legend(h,clengend_final,'interpreter','latex','fontsize',16,'Location','north',Box='off')


exportgraphics(gca,strcat('Figures/','Optimal_Output_',name_save,'_alpha_',num2str(alpha_target),'_beta_',num2str(beta_target), ...
    '_Br_',num2str(Br_target),'_Re_',num2str(Re_target),'.jpeg'),'Resolution',300)

end
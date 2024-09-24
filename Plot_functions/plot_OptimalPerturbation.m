function plot_OptimalPerturbation(y, delta, q_vars, q_in, q_out, pos_in, pos_out, name_file_load, Br, alpha, beta)

% Smooth and figure sample show
i_smooth = 3;
i_sample = 2;

% Plot q_in
for n_var = 1:q_vars
    Vec_eigen_var{n_var} = q_in(n_var:q_vars:end); % Take all columns for each perturbation row
end

Norm = max(abs(Vec_eigen_var{pos_in}));

figure;
hold on; grid on; box on;
plot(y(1:i_sample:end)/delta,smooth(abs(Vec_eigen_var{1}(1:i_sample:end))/Norm,i_smooth),'LineWidth',2,'LineStyle',':','color',[0.4660 0.6740 0.1880]);
plot(y(1:i_sample:end)/delta,smooth(abs(Vec_eigen_var{2}(1:i_sample:end))/Norm,i_smooth),'LineWidth',2,'LineStyle','-','color',[0 0.4470 0.7410])
plot(y(1:i_sample:end)/delta,smooth(abs(Vec_eigen_var{3}(1:i_sample:end))/Norm,i_smooth),'LineWidth',2,'LineStyle','--','color',[0.8500 0.3250 0.0980])
plot(y(1:i_sample:end)/delta,smooth(abs(Vec_eigen_var{4}(1:i_sample:end))/Norm,i_smooth),'LineWidth',2,'LineStyle','-.','color','black')
plot(y(1:i_sample:end)/delta,smooth(abs(Vec_eigen_var{5}(1:i_sample:end))/Norm,i_smooth),'LineWidth',2,'LineStyle','-.','color',[0.4940 0.1840 0.5560])
xlabel('${y/\delta}$','interpreter','latex')
ylabel('$|q^{\prime}|$','interpreter','latex')
pbaspect([1 1.2 1])


legend('$|\rho^{\prime}|$','$|u^{\prime}|$','$|v^{\prime}|$', '$|w^{\prime}|$','$|T^{\prime}|$','interpreter','latex')
legend('Location','north','box','off')
ylim([0 1])

set(gca,'linewidth',1.5)
set(gca,'fontsize',14)



exportgraphics(gca,strcat('Figures/',name_file_load,'_OptimalPerturbation_IN_Br_',num2str(Br),'_Alpha_',num2str(alpha),'_Beta_',num2str(beta),'.jpeg'),'Resolution',300)


% Plot q_out
for n_var = 1:q_vars
    Vec_eigen_var{n_var} = q_out(n_var:q_vars:end); % Take all columns for each perturbation row
end

% Smooth and figure sample show
i_smooth = 1;
i_sample = 1;

Norm = max(abs(Vec_eigen_var{pos_out}));
delta = 1;
figure;
hold on; grid on; box on;
plot(y(1:i_sample:end)/delta,smooth(abs(Vec_eigen_var{1}(1:i_sample:end))/Norm,i_smooth),'LineWidth',2,'LineStyle',':','color',[0.4660 0.6740 0.1880]);
plot(y(1:i_sample:end)/delta,smooth(abs(Vec_eigen_var{2}(1:i_sample:end))/Norm,i_smooth),'LineWidth',2,'LineStyle','-','color',[0 0.4470 0.7410])
plot(y(1:i_sample:end)/delta,smooth(abs(Vec_eigen_var{3}(1:i_sample:end))/Norm,i_smooth),'LineWidth',2,'LineStyle','--','color',[0.8500 0.3250 0.0980])
plot(y(1:i_sample:end)/delta,smooth(abs(Vec_eigen_var{4}(1:i_sample:end))/Norm,i_smooth),'LineWidth',2,'LineStyle','-.','color','black')
plot(y(1:i_sample:end)/delta,smooth(abs(Vec_eigen_var{5}(1:i_sample:end))/Norm,i_smooth),'LineWidth',2,'LineStyle','-.','color',[0.4940 0.1840 0.5560])
xlabel('${y/\delta}$','interpreter','latex')
ylabel('$|q^{\prime}|$','interpreter','latex')
pbaspect([1 1.2 1])

legend('$|\rho^{\prime}|$','$|u^{\prime}|$','$|v^{\prime}|$', '$|w^{\prime}|$','$|T^{\prime}|$','interpreter','latex')
legend('Location','north','box','off')
ylim([0 1])

set(gca,'linewidth',1.5)
set(gca,'fontsize',14)

exportgraphics(gca,strcat('Figures/',name_file_load,'_OptimalPerturbation_OUT_Br_',num2str(Br),'_Alpha_',num2str(alpha),'_Beta_',num2str(beta),'.jpeg'),'Resolution',300)


end
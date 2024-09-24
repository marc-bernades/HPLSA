function plot_Spectrum(y, delta, N, alpha, N_target,n_LST_Sweep, Val_eigen,  Vec_eigen_var, name_file_load)

realc = real(Val_eigen);
imagc = imag(Val_eigen);
disp(['Most unstable mode = ', num2str(max(imagc)), ' total of = ', num2str(sum(imagc>0))])

% Split thermodynamic and dynamic modes
q_vars = 5;
Vec_eigen_var_int = zeros(1,q_vars);
bmode             = true(length(Vec_eigen_var{1}(1,:)),1);

% Sweep all eigen values
for n_mode = 1:length(Vec_eigen_var{1}(1,:))
    % Integration of the perturbation across y for each
    for n_var = 1:q_vars
        Vec_eigen_var_int(n_var) = abs(trapz(y/delta,abs(Vec_eigen_var{n_var}(:,n_mode))));
    end

    if Vec_eigen_var_int(1) < 1E-3 && Vec_eigen_var_int(5) < 1E-3
        bmode(n_mode) = false; % Dynamic (rho' = T' = 0)
    else
        bmode(n_mode) = true; % Thermo
    end
end

figure
hold on; grid on; box on
plot(realc(bmode),imagc(bmode),'o','color',[0 0.4470 0.7410],'MarkerSize',5,'MarkerFaceColor',[0 0.4470 0.7410])
plot(realc(~bmode),imagc(~bmode),'s','color','red','MarkerSize',5,'MarkerFaceColor','red')
legend('Thermo','Dynamic','Location','best')
try
xlim([ceil(min(realc)) ceil(max(realc(realc~=1E4)))])
catch
xlim([(min(realc))-0.1 (max(realc(realc~=1E4)))+0.1])
end
%     axis([0 1 -1 0])
xlabel('${\omega}_r$','interpreter','latex')
ylabel('${\omega}_i$','interpreter','latex')
set(gca,'fontsize',12)
set(gca,'linewidth',1.5)

axes('position',[.50 .2 .35 .3])

plot(realc(bmode),imagc(bmode),'o','color',[0 0.4470 0.7410],'MarkerSize',5,'MarkerFaceColor',[0 0.4470 0.7410])
hold on; grid on; box on;
plot(realc(~bmode),imagc(~bmode),'s','color','red','MarkerSize',5,'MarkerFaceColor','red')
xlim([0 1.2])
ylim([-1 ,0.2])
set(gca,'linewidth',1)
set(gca,'fontsize',10)


% saveas(gca,strcat('Figures/', name_file_load,'_Spectrum_PrEc_',num2str(N_target),'_alpha_',num2str(alpha), '_Re_',num2str(n_LST_Sweep),'_N_', num2str(N)),'epsc')
% saveas(gca,strcat('Figures/', name_file_load,'_Spectrum_PrEc_',num2str(N_target),'_alpha_',num2str(alpha), '_Re_',num2str(n_LST_Sweep),'_N_', num2str(N)),'png')

end
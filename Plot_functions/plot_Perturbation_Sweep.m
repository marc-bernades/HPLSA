function plot_Perturbation_Sweep(y, delta, N, alpha, N_target, n_LST_Sweep, alpha_plot, Re_plot, Val_eigen, Vec_eigen, Val_eigen_var, Vec_eigen_var, name_file_load)

Color_map = {[0 0.4470 0.7410],[0.4660 0.6740 0.1880],[0.8500 0.3250 0.0980],...
    [0.4940 0.1840 0.5560], [0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840], ...
    [0.9290 0.6940 0.1250], [1, 0, 0], [0, 0, 1], [0, 1, 0], [1, 1, 0], [0, 0.5, 0.5]};


figure;
hold on; grid on; box on;

% Breakdown eigen vector for each variable
q_vars = 5;
Val_eigen_var = cell(1,q_vars);
Vec_eigen_var = cell(1,q_vars);

for ii = 1:length(N_target)

[~,aa] = min(abs(alpha-alpha_plot));
[~,jj] = min(abs(n_LST_Sweep-Re_plot));

for n_var = 1:q_vars
    Val_eigen_var{n_var} = Val_eigen{aa,ii,jj}(n_var:q_vars:end);
    Vec_eigen_var{n_var} = Vec_eigen{aa,ii,jj}(n_var:q_vars:end,:); % Take all columns for each perturbation row
end

% Perturbation target
realc = real(Val_eigen{aa,ii,jj});
imagc = imag(Val_eigen{aa,ii,jj});
Cond   = abs(realc) <= 10 & realc < 0.5;
Val_eigen_target   = max(imagc(Cond));

% Plot eigen vector with respect to y for each perturbation
[~,idx]            = min(abs(imag(Val_eigen{aa,ii,jj}) - Val_eigen_target));

% Normalize by u'
Norm = max(abs(Vec_eigen_var{2}(:,idx)));

plot(y/delta,abs(Vec_eigen_var{1}(:,idx))/Norm,'LineWidth',2,'LineStyle','-','color', Color_map{ii});
plot(y/delta,abs(Vec_eigen_var{2}(:,idx))/Norm,'LineWidth',2,'LineStyle',':','color', Color_map{ii});
plot(y/delta,abs(Vec_eigen_var{3}(:,idx))/Norm,'LineWidth',2,'LineStyle','-.','color',Color_map{ii});
plot(y/delta,abs(Vec_eigen_var{5}(:,idx))/Norm,'LineWidth',2,'LineStyle','--','color',Color_map{ii});

end

xlabel('${y/\delta}$','interpreter','latex')
ylabel('$|q^{\prime}|$','interpreter','latex')
% ylim([0 10])

set(gca,'linewidth',1.5)
set(gca,'fontsize',14)

% Legend perturbation vector and Br sweep
% legend('$|\rho^{\prime}|$','$|u^{\prime}|$','$|v^{\prime}|$', '$|T^{\prime}|$','interpreter','latex','fontsize',14)
% legend('Location','best','box','off')

h(1) = plot(NaN,'-','color','k','LineWidth',2);
h(2) = plot(NaN,':','color','k','LineWidth',2);
h(3) = plot(NaN,'-.','color','k','LineWidth',2);
h(4) = plot(NaN,'--','color','k','LineWidth',2);

for kk = 1:length(N_target)
    h(kk+4) = plot(NaN,'-','color',Color_map{kk},'LineWidth',2);
    clegend{kk} = strcat('$', 'Br', ' = ',num2str(N_target(kk)), '$');
end
clengend_final = {'$|\rho^{\prime}|$','$|u^{\prime}|$','$|v^{\prime}|$', '$|T^{\prime}|$',clegend{:}};
legend(h,clengend_final,'interpreter','latex','fontsize',16,'Location','north',Box='off')

% ylim([0 2])

set(gca,'linewidth',1.5)
set(gca,'fontsize',16)

exportgraphics(gca,strcat('Figures/',name_file_load,'_alpha_',num2str(alpha(aa)), '_Re_',num2str(n_LST_Sweep(jj)),'.jpeg'),'Resolution',300)


% G_max_var     = cell(1,q_vars);
% for n_var = 1:q_vars
%     for aa = 1:length(alpha)
%         for ii = 1:length(N_target)
%             G_max_var{n_var}(aa,ii)       = max(G{aa,ii}(n_var:5:end));
%         end
%     end
% end

end
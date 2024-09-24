
function plot_Baseflow_LST_Comparison(y_base,delta_base,BF_base,N_target_base, PrEc, name_file_load, name_file_load_2)

Color_map = {[0 0.4470 0.7410],[0.4660 0.6740 0.1880],[0.4940 0.1840 0.5560],...
    [0.8500 0.3250 0.0980], [0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840], ...
    [0.9290 0.6940 0.1250], [1, 0, 0], [0, 0, 1], [0, 1, 0]};

% Unpack second setup
Results = load(strcat(name_file_load_2, '.mat'));

% Unpack structure and load variables
struct_field_names = fieldnames(Results.Data);
for k = 1:length(struct_field_names)
      eval([struct_field_names{k}, ' = Results.Data.',struct_field_names{k}, ';']);
end

% Position
[~,ii] = min(abs(N_target_base-PrEc));
[~,jj] = min(abs(N_target-PrEc));

% Velocity
figure; hold on; box on; grid on;

plot(y_base/delta_base,BF_base{ii}.u/BF_base{ii}.norm.u,'linewidth',2, 'color',Color_map{1});
clabel{1} = strcat('$', 'PrEc = ',num2str(N_target_base(ii)), '$');
plot(y/delta,BF{jj}.u/BF{jj}.norm.u,'linewidth',2, 'color',Color_map{2});
clabel{2} = strcat('$', 'PrEc = ',num2str(N_target(jj)), '$');
xlabel('${{y}/{\delta}}$','interpreter','latex')
ylabel('${{u}/{u_b}}$','interpreter','latex')
legend(clabel,'interpreter','latex', 'location','south')
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)
pbaspect([1.8 1 1])

saveas(gca,strcat('Figures/', name_file_load,'_Baseflow_u_PrEc_',num2str(N_target(ii)),'_Comparison'),'epsc')
saveas(gca,strcat('Figures/', name_file_load,'_Baseflow_u_PrEc_',num2str(N_target(ii)),'_Comparison'),'png')


% Temperature
figure; hold on; box on; grid on;
plot(y_base/delta_base,BF_base{ii}.T/BF_base{ii}.T_c, 'linewidth',2,'color',Color_map{1});
plot(y/delta,BF{jj}.T/T_c, 'linewidth',2,'color',Color_map{2});
xlabel('${{y}/{\delta}}$','interpreter','latex')
ylabel('${{T}/{T_c}}$','interpreter','latex')
legend(clabel,'interpreter','latex', 'location','south')
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)
pbaspect([1.8 1 1])

saveas(gca,strcat('Figures/',name_file_load,'_Baseflow_T_PrEc_',num2str(N_target(ii)),'_Comparison'),'epsc')
saveas(gca,strcat('Figures/',name_file_load,'_Baseflow_T_PrEc_',num2str(N_target(ii)),'_Comparison'),'png')


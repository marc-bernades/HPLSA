function plot_TransientGrowth_Time(t_all, GG, alpha, alpha_target, beta, beta_target, n_LST_Sweep, Re_target,N_target, N_plot, name_file_load, varargin)

Color_map = {[0 0.4470 0.7410],[0.4660 0.6740 0.1880],[0.8500 0.3250 0.0980],...
    [0.4940 0.1840 0.5560], [0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840], ...
    [0.9290 0.6940 0.1250], [1, 0, 0], [0, 0, 1], [0, 1, 0], [1, 1, 0], [0, 0.5, 0.5]};


[~,aa] = min(abs(alpha-alpha_target));
[~,bb] = min(abs(beta-beta_target));
[~,jj] = min(abs(n_LST_Sweep-Re_target));

f = figure;
hold on; grid on; box on;
idx = 0;

for ii = 1:length(N_target)
    if ~isempty(find(N_target(ii) == N_plot))
        idx = idx + 1;
        plot(smooth(t_all{aa,bb,ii,jj}(1:2:end),3), smooth(GG{aa,bb,ii,jj}(1:2:end),3), 'LineWidth',2,'LineStyle','-','color', Color_map{idx});
        clabel{idx} = strcat('$', 'Br = ',num2str(N_target(ii)), '$');
    end
end

xlabel('${t}$','interpreter','latex')
ylabel('$G(t)$','interpreter','latex')
% ylim([0 10])
legend(clabel,'interpreter','latex', 'location','south','box','off')
pbaspect([1 1.5 1])

set(gca,'linewidth',1.5)
set(gca,'fontsize',14)

% GG_max = round(max(GG{aa,bb,ii,jj}),-2);
% yticks([linspace(0,GG_max,5)]);


% Inset?
if ~isempty(varargin)
    [~,aa] = min(abs(alpha-varargin{1}));
    [~,bb] = min(abs(beta-varargin{2}));
    axes('position',[.45 0.35 .25 .25]) % NI-6 Re = 10000
    % axes('position',[.4 .58 .35 .3]) % NI-6 Re = 4000
    hold on; grid off; box on;

    idx = 0;
    for ii = 1:length(N_target)
        if ~isempty(find(N_target(ii) == N_plot))
            idx = idx + 1;
            plot(smooth(t_all{aa,bb,ii,jj}(1:2:end),3), smooth(GG{aa,bb,ii,jj}(1:2:end),3), 'LineWidth',2,'LineStyle','-','color', Color_map{idx});
            clabel{idx} = strcat('$', 'Br = ',num2str(N_target(ii)), '$');
        end
    end
    xticklabels({});

set(gca,'linewidth',1.5)
set(gca,'fontsize',10)
end



exportgraphics(f,strcat('Figures/',name_file_load,'_alpha_',num2str(alpha(aa)),'_beta_',num2str(beta(bb)), '_Re_',num2str(n_LST_Sweep(jj)),'.jpeg'),'Resolution',300)



end
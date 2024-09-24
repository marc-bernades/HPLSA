function plot_DimensionlessNumbers(BF, N_target, name_file_load)

Re = []; PrEc = []; Pr = []; Ec = []; Ma = [];
for jj = 1:length(N_target)
    Re(jj) = BF{jj}.Re;
    Pr(jj) = BF{jj}.Pr;
    Ec(jj) = BF{jj}.Ec;
    Ma(jj) = BF{jj}.Ma;
    PrEc(1,jj) = BF{jj}.Pr*BF{jj}.Ec;
end

figure
loglog(Re,PrEc,'o-')
xlabel('${Re}$','interpreter','latex')
ylabel('${PrEc}$','interpreter','latex')
set(gca,'fontsize',12)
set(gca,'linewidth',1.5)
grid on
saveas(gca,strcat('Figures/',name_file_load,'_Re_vs_PrEc'),'epsc')
saveas(gca,strcat('Figures/',name_file_load,'_Re_vs_PrEc'),'png')


disp("Re = "   + num2str(Re))
disp("Pr = "   + num2str(Pr))
disp("Ec = "   + num2str(Ec))
disp("PrEc = " + num2str(PrEc))
disp("Ma = "   + num2str(Ma))



end
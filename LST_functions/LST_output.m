function LST_output(Fluid, T_bw, T_tw, T_c, P_b, P_c, N, y, D, delta, bScaling, bTarget, N_target, n_LST_Sweep, alpha, beta,BF,Dy,Vec_eigen,Val_eigen, G, G_max, TestFolder , varargin)

Data.Fluid = Fluid;
Data.T_bw = T_bw;
Data.T_tw = T_tw;
Data.P_b = P_b;
Data.N = N;
Data.y = y;
Data.D = D;
Data.delta = delta;
Data.bScaling = bScaling;
Data.bTarget = bTarget;
Data.N_target = N_target;
Data.n_LST_Sweep = n_LST_Sweep;
Data.alpha = alpha;
Data.beta = beta;
Data.BF = BF;
Data.Dy = Dy;
Data.Val_eigen = Val_eigen;
Data.Vec_eigen = Vec_eigen;
Data.G = G;
Data.G_max = G_max;


% Name
if ~isempty(varargin)
    name_file_out = strcat(Fluid.Substance, '_Tbw_', num2str(T_bw), '_Ttw_', num2str(T_tw), '_Pb_', num2str(P_b), '_', bTarget, '_', bScaling, '_' , varargin{1});
else
    name_file_out = strcat(Fluid.Substance, '_Tbw_', num2str(T_bw), '_Ttw_', num2str(T_tw), '_Pb_', num2str(P_b), '_', bTarget, '_', bScaling);
end
% Write the table to a MAT file
save("Results/" + TestFolder + "/"+ name_file_out + ".mat",'Data', '-v7.3')
disp("File " + name_file_out + " saved succcessfully...")

end
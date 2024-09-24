function [BF, Dy]  = Baseflow_LST_adimensional(y,D,N,delta,L_y,T_bw,T_tw,P_b,bSolver,HP_model, Substance,Fluid,bTarget,bScaling,varargin)

%% Check Refprop
if strcmp(HP_model,'RefProp')
    Substance = strcat('REFPROP::',Substance);
end

%% Initialization
% Base flow
u_0     = y.*(L_y - y); % Normalized Laminar parabolic profile assumption
T_0     = T_bw + y./L_y*(T_tw - T_bw); %T_bw + T_bw.*y.*(1 - y); % 1*(T_bw + (T_tw - T_bw)*tanh(4*y)); % T_bw + y./L_y*(T_tw - T_bw); % Linear initialization
P_b     = zeros(size(y)) + P_b;
P_0     = P_b;

% Boundary conditions Dirichlet
u_0(1)   = 0;
u_0(end) = 0;
T_0(1)   = T_bw;
T_0(end) = T_tw;

% Calculate thermodynamic fields at initial state
if ~strcmp(HP_model,'CoolProp') && ~strcmp(HP_model,'RefProp')
    rho_0   = Calculate_Rho_from_TandP(    bSolver, T_0,P_0, Fluid, Substance );
else
    for nCP = 1:length(u_0)
        rho_0(nCP,1)    = py.CoolProp.CoolProp.PropsSI("D","T",T_0(nCP),"P",P_0(nCP),Substance);
    end
end

if ~strcmp(HP_model,'CoolProp') && ~strcmp(HP_model,'RefProp')
    e_0     = Calculate_e_from_TandPandRho(bSolver, T_0,rho_0,P_0,Fluid, Substance );
else
    for nCP = 1:length(u_0)
        e_0(nCP,1)    = py.CoolProp.CoolProp.PropsSI("U","T",T_0(nCP),"P",P_0(nCP),Substance);
    end
end

% mu and k from high-pressure model
if ~strcmp(HP_model,'CoolProp') && ~strcmp(HP_model,'RefProp')
    [mu_0,kappa_0] = Calculate_HighPressure_Transport_Coeff(bSolver, zeros(size(u_0)), zeros(size(u_0)), Fluid, Substance, T_0, rho_0, HP_model);
else
    for nCP = 1:length(T_0)
        kappa_0(nCP,1)     = py.CoolProp.CoolProp.PropsSI("L","T",T_0(nCP),"P",P_0(nCP),Substance);  	% Thermal conductivity at film temperature [W/mK]
        mu_0(nCP,1)        = py.CoolProp.CoolProp.PropsSI("V","T",T_0(nCP),"P",P_0(nCP),Substance);	 	% Dynamic viscosity at film temperature [kg/ms]
    end
end

% Heat capacities and sos
if ~strcmp(HP_model,'CoolProp') && ~strcmp(HP_model,'RefProp')
    [c_p_0,c_v_0,~]      = Calculate_SpecificHeatCapacities(bSolver, P_0,rho_0,T_0, Fluid,Substance);
    [sos_0, ~, ~, ~, ~]  = Calculate_sos(bSolver,rho_0,T_0,c_p_0,P_0,0,Fluid,Substance);
else
    for nCP = 1:length(u_0)
        c_p_0(nCP,1)    = py.CoolProp.CoolProp.PropsSI("C","T",T_0(nCP),"P",P_0(nCP),Substance);  	% Specific heat at constant pressure [J/kg/K]
        c_v_0(nCP,1)    = py.CoolProp.CoolProp.PropsSI("CVMASS","T",T_0(nCP),"P",P_0(nCP),Substance);  	% Specific heat at constant pressure [J/kg/K]
        sos_0(nCP,1)    = py.CoolProp.CoolProp.PropsSI("A","T",T_0(nCP),"P",P_0(nCP),Substance);  	% Speed of sound [m/s]
    end
end
%% Define targets
PrEc_0  = varargin{1}; % Target
if strcmp(bScaling,'wall')
    u_c     = sqrt(PrEc_0*kappa_0(1)*T_0(1)/mu_0(1)); % Normalization wall
else
    [~,~,T_b,mu_b,kappa_b,~,~,~,~,~] = CalculateBulkFields(y,rho_0,u_0,T_0,mu_0,kappa_0,c_p_0,c_v_0,e_0,e_0,sos_0);
    u_c     = sqrt(PrEc_0*kappa_b*T_b/mu_b); % Bulk
end

F_hat_0 = 2; %1/Re_0*mean(-CentralDerivative_d1_2ndOrder_1D(mu_0.*CentralDerivative_d1_2ndOrder_1D(u_0,y)',y)');
% F_hat_0 for baseflow laminar validation
% p       = [ -0.2558    1.4032   -3.0526    3.6797    1.2474];
% F_c     = 0.715; % Factor to adjust bulk velocity profile as per DNS
% F_hat_0 = F_c*polyval(p,y/delta);

u_0     = u_0*u_c;
E_0     = e_0 + 1/2*(u_0.^2);


%% Dimensionless parameters at initial state
% Bulk fields
[rho_b,u_b,T_b,mu_b,kappa_b,c_p_b,c_v_b,E_b,e_b,~] = CalculateBulkFields(y,rho_0,u_0,T_0,mu_0,kappa_0,c_p_0,c_v_0,E_0,e_0,sos_0);

% Second viscosity
lambda_0 = 0*mu_b - 2/3*mu_0;


%% Normalize parameters
% Normalize reference scalings
Dy.du_dy = D*u_0;
norm     = Reference_Scalings(u_b,T_b,rho_b,mu_b,kappa_b,c_p_b,c_v_b,E_b,e_b,P_b(1),u_0,T_0,rho_0,mu_0,kappa_0,c_p_0,c_v_0,lambda_0,Dy,delta,bScaling, u_c);

% Normalize variables
[u_star_0, T_star_0, ~, ~, mu_star_0, kappa_star_0, ~, ~, ~, ~] = Normalize(norm, u_0, T_0, rho_0, P_0, mu_0, kappa_0, c_p_0, c_v_0, E_0, e_0, lambda_0, y);


[y_star,D_star] = Chevysheb_collocation(N,delta/norm.y);

%% Base flow solver iterative integration method
%Initialize time
time  = 0;
n_max = 1E2;
dt    = 0.05;

% Norm error checks
u_L2_norm = zeros(1,n_max);
T_L2_norm = zeros(1,n_max);


% Iteration
for t = 1:n_max

    try

        % Calculate time step
%         disp("Computing iteration " + num2str(t) + " at t = " + num2str(time) + "s" + " dt = " + num2str(dt))

        % Calculate integrated quantities
        u_eq{1} = -F_hat_0.*cumtrapz(y_star,y_star./mu_star_0);
        u_eq{2} = cumtrapz(y_star,1./mu_star_0);
        C_u{1}  = -(u_eq{1}(1) - u_eq{1}(end))./(u_eq{2}(1) - u_eq{2}(end));
        C_u{2}  = -u_eq{1}(1) - C_u{1}.*u_eq{2}(1);

        u_star       = u_eq{1} + u_eq{2}.*C_u{1} + C_u{2};
        du_dy_star   = D_star*u_star; %CentralDerivative_d1_2ndOrder_1D(u,y)'; % This is equivalent to D*u

        T_eq{1} = -F_hat_0.*PrEc_0.*cumtrapz(y_star,cumtrapz(y_star,u_star)./kappa_star_0);
        T_eq{2} =  PrEc_0*cumtrapz(y_star,1./kappa_star_0);
        T_eq{3} = -PrEc_0*cumtrapz(y_star,u_star.*mu_star_0./kappa_star_0.*du_dy_star);
        C_T{1}  = (T_bw/norm.T - T_tw/norm.T - (T_eq{1}(1) - T_eq{1}(end)) - (T_eq{3}(1) - T_eq{3}(end)))./(T_eq{2}(1) - T_eq{2}(end));
        C_T{2}  = T_bw/norm.T - T_eq{1}(1) - C_T{1}*T_eq{2}(1) - T_eq{3}(1);

        T_star  = T_eq{1} + C_T{1}.*T_eq{2} + T_eq{3} + C_T{2};

        % Boundary conditions
        u_star(1) = 0; u_star(end) = 0;
        T_star(1) = T_bw/norm.T; T_star(end) = T_tw/norm.T;

        % Dimensionalized variables
        T        = T_star*norm.T;       

        % Update density with correct boundary condition fields
        if ~strcmp(HP_model,'CoolProp') && ~strcmp(HP_model,'RefProp')
            rho   = Calculate_Rho_from_TandP( bSolver,T,P_0, Fluid, Substance );
        else
            for nCP = 1:length(u_star)
                rho(nCP,1)    = py.CoolProp.CoolProp.PropsSI("D","T",T(nCP),"P",P_0(nCP),Substance);
            end
        end
        

        % High-pressure coefficients
        if ~strcmp(HP_model,'CoolProp') && ~strcmp(HP_model,'RefProp')
            [mu,kappa]        = Calculate_HighPressure_Transport_Coeff(bSolver, mu_0, kappa_0, Fluid, Substance, T, rho, HP_model);
            [c_p,c_v,gamma]   = Calculate_SpecificHeatCapacities(bSolver,P_0,rho,T, Fluid,Substance);
            [sos, ~, ~, ~, ~] = Calculate_sos(bSolver,rho,T,c_p,P_0,0,Fluid,Substance);
        else
            for nCP = 1:length(c_p_0)
                kappa(nCP,1)  = py.CoolProp.CoolProp.PropsSI("L","T",T(nCP),"P",P_0(nCP),Substance);  	% Thermal conductivity at film temperature [W/mK]
                mu(nCP,1)     = py.CoolProp.CoolProp.PropsSI("V","T",T(nCP),"P",P_0(nCP),Substance);	 	% Dynamic viscosity at film temperature [kg/ms]
                c_p(nCP,1)    = py.CoolProp.CoolProp.PropsSI("C","T",T(nCP),"P",P_0(nCP),Substance);  	% Specific heat at constant pressure [J/kg/K]
                c_v(nCP,1)    = py.CoolProp.CoolProp.PropsSI("CVMASS","T",T(nCP),"P",P_0(nCP),Substance); % Specific heat at constant pressure [J/kg/K]
                sos(nCP,1)    = py.CoolProp.CoolProp.PropsSI("A","T",T(nCP),"P",P_0(nCP),Substance);  	% Speed of sound [m/s]
            end
            gamma = c_p./c_v;
        end

        % Calculate reference velocity
        if strcmp(bScaling,'wall')
            u_c     = sqrt(PrEc_0*kappa(1)*T(1)/mu(1)); % Normalization wall
        else
            [~,u_b,T_b,mu_b,kappa_b,~,~,~,~,~] = CalculateBulkFields(y,rho,u_0,T,mu,kappa,c_p,c_v,sos,sos,sos); % Bulk (impose E = e = sos for utility)
            u_c     = sqrt(PrEc_0*kappa_b*T_b/mu_b); % Normalization bulk
        end


        % Update velocity dimensional
        u        = u_star.*u_c;
        Dy.du_dy = D*u;

        % Update energies
        ke   = 1/2*(u.^2);
        if ~strcmp(HP_model,'CoolProp') && ~strcmp(HP_model,'RefProp')
            e         = Calculate_e_from_TandPandRho(bSolver, T,rho,P_0,Fluid, Substance );
        else
            for nCP = 1:length(u)
                e(nCP,1)  = py.CoolProp.CoolProp.PropsSI("U","T",T(nCP),"P",P_0(nCP),Substance);
            end
        end
        E    = e + ke;
        
        % Update bulk values and dimensionless numbers
        [rho_b,u_b,T_b,mu_b,kappa_b,c_p_b,c_v_b,E_b,e_b,sos_b] = CalculateBulkFields(y,rho,u,T,mu,kappa,c_p,c_v,E,e,sos);
        [Re,Pr,Ec,Ma] = DimensionlessNumbers(rho_b,u_b,T_b,mu_b,kappa_b,c_p_b,sos_b,rho,u,T,mu,kappa,c_p,sos,delta,Dy.du_dy,bScaling, u_c);   
        
        % Second viscosity
        lambda = 0*mu_b - 2/3*mu;

        % Normalize
        norm     = Reference_Scalings(u_b,T_b,rho_b,mu_b,kappa_b,c_p_b,c_v_b,E_b,e_b,P_b(1),u,T,rho,mu,kappa,c_p,c_v,lambda,Dy,delta,bScaling, u_c);
        [u_star, T_star, ~, ~, mu_star, kappa_star, ~, ~, ~, ~, ~, ~] = Normalize(norm, u, T, rho, P_b, mu, kappa, c_p, c_v, E, e, lambda, y);

        % Error tolerance
        u_L2_norm(t) = sqrt(sum((u_star_0 - u_star).^2)./sum(u_star_0.^2));
        T_L2_norm(t) = sqrt(sum((T_star_0 - T_star).^2)./sum(T_star_0.^2));

        if u_L2_norm(t) < 1E-15 && T_L2_norm(t) < 1E-15
            break
        end


%         disp(['rhou norm = ', num2str(u_L2_norm(t)), ' and T norm = ', num2str(T_L2_norm(t))])
        

        % Update current state
        u_star_0     = u_star;
        T_star_0     = T_star;
        mu_star_0    = mu_star;
        kappa_star_0 = kappa_star;
     
        % Update time
        time = time + dt;

%         drawnow
%         plot(y,T); hold on

    catch CA
        CA.message
        break
    end

end

disp("Final iteration " + num2str(t) + " at t = " + num2str(time) + "s" + " dt = " + num2str(dt))
disp(['rhou norm = ', num2str(u_L2_norm(t)), ' and T norm = ', num2str(T_L2_norm(t))])


% Calculate 1st- and 2nd- central derivatves O(2)
Dy.D          = D;
Dy.du_dy      = D*u;      %CentralDerivative_d1_2ndOrder_1D(u,y);
Dy.dT_dy      = D*T;      %CentralDerivative_d1_2ndOrder_1D(T,y);
Dy.dmu_dy     = D*mu;     %CentralDerivative_d1_2ndOrder_1D(mu,y);
Dy.dkappa_dy  = D*kappa;  %CentralDerivative_d1_2ndOrder_1D(kappa,y);
Dy.drho_dy    = D*rho;    %CentralDerivative_d1_2ndOrder_1D(rho,y);
Dy.dlambda_dy = D*lambda; %CentralDerivative_d1_2ndOrder_1D(lambda,y);
Dy.de_dy      = D*e;      %CentralDerivative_d1_2ndOrder_1D(e,y);
Dy.du2_dy2    = D*D*u;    %CentralDerivative_d2_2ndOrder_1D(u,y);
Dy.dT2_dy2    = D*D*T;    %CentralDerivative_d2_2ndOrder_1D(T,y);

% Load params into structure
BF.u_star    = u_star;
BF.T_star    = T_star;
BF.mu_star   = mu_star;
BF.kappa_star= kappa_star;
BF.u         = u;
BF.T         = T;
BF.rho       = rho;
BF.P_0       = P_b;
BF.mu        = mu;
BF.kappa     = kappa;
BF.c_p       = c_p;
BF.c_v       = c_v;
BF.e         = e;
BF.E         = E;
BF.sos       = sos;
BF.ke        = ke;
BF.Re        = Re;
BF.Pr        = Pr;
BF.Ec        = Ec;
BF.Ma        = Ma;
BF.gamma     = gamma;
BF.lambda    = lambda;
BF.u_b       = u_b;
BF.T_b       = T_b;
BF.rho_b     = rho_b;
BF.mu_b      = mu_b;
BF.kappa_b   = kappa_b;
BF.c_p_b     = c_p_b;
BF.c_v_b     = c_v_b;
BF.E_b       = E_b;
BF.e_b       = e_b;
BF.sos_b     = sos_b;
BF.P_b       = P_b(1);
BF.norm      = norm;
BF.y         = y;
F_hat_0 = 2;
BF.F_hat_0   = F_hat_0;


end
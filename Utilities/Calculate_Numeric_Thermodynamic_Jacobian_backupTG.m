function J_t2 = Calculate_Numeric_Thermodynamic_Jacobian(bSolver,c_v,T,rho,Fluid, Substance, HP_model)

dcv_dT         = zeros(size(T));
dmu_dT         = zeros(size(T));
d2mu_d2T       = zeros(size(T));
dkappa_dT      = zeros(size(T));
d2kappa_d2T    = zeros(size(T));
dmu_drho       = zeros(size(T));
d2mu_d2rho     = zeros(size(T));
dkappa_drho    = zeros(size(T));
d2kappa_d2rho  = zeros(size(T));
dmu_drho_1     = zeros(size(T));
dmu_drho_2     = zeros(size(T));
d2mu_drhodT    = zeros(size(T));
dkappa_drho_1  = zeros(size(T));
dkappa_drho_2  = zeros(size(T));
d2kappa_drhodT = zeros(size(T));
dmu_dT_1       = zeros(size(T));
dmu_dT_2       = zeros(size(T));
d2mu_dTdrho    = zeros(size(T));
dkappa_dT_1    = zeros(size(T));
dkappa_dT_2    = zeros(size(T));
d2kappa_dTdrho = zeros(size(T));

stencil_weight = 10^-3;
delta_T        = 0.01; %0.01; % Nitrogen at bulk this needs to be 0.05 to obtain stable mode and 0.1 unstable mode
delta_rho      = 0.01;


for ii = 1:length(c_v)

    if ~strcmp(HP_model, 'CoolProp')

        %% Partial temperature derivatives
        T_prev = T(ii) - T(ii)*stencil_weight;
        T_next = T(ii) + T(ii)*stencil_weight;

        if strcmp(bSolver,'Real')
            P_prev   = Calculate_P_PengRobinson(rho(ii),T_prev,Substance);
            P_next   = Calculate_P_PengRobinson(rho(ii),T_next,Substance);
        else
            P_prev   = rho(ii)*T_prev*Fluid.R_specific;
            P_next   = rho(ii)*T_next*Fluid.R_specific;
        end

        [~,c_v_prev,~] = Calculate_SpecificHeatCapacities(bSolver, P_prev,rho(ii),T_prev, Fluid,Substance);
        [~,c_v_next,~] = Calculate_SpecificHeatCapacities(bSolver, P_next,rho(ii),T_next, Fluid,Substance);

        dcv_dT(ii) = (c_v_next - c_v_prev)./(T_next - T_prev);


        [mu_prev,kappa_prev]       = Calculate_HighPressure_Transport_Coeff(bSolver, 0, 0, Fluid, Substance, T_prev, rho(ii), HP_model);
        [mu_next,kappa_next]       = Calculate_HighPressure_Transport_Coeff(bSolver, 0, 0, Fluid, Substance, T_next, rho(ii), HP_model);
        [mu_current,kappa_current] = Calculate_HighPressure_Transport_Coeff(bSolver, 0, 0, Fluid, Substance, T(ii), rho(ii), HP_model);


        dmu_dT(ii)        = (mu_next - mu_prev)./(T_next - T_prev);
        dkappa_dT(ii)     = (kappa_next - kappa_prev)./(T_next - T_prev);
        d2mu_d2T(ii)      = (mu_next -2*mu_current + mu_prev)./(0.5*(T_next - T_prev)).^2;
        d2kappa_d2T(ii)   = (kappa_next -2*kappa_current + kappa_prev)./(0.5*(T_next - T_prev)).^2;

        %% Partial density derivatives
        rho_prev = rho(ii) - rho(ii)*stencil_weight;
        rho_next = rho(ii) + rho(ii)*stencil_weight;

        [mu_prev,kappa_prev]       = Calculate_HighPressure_Transport_Coeff(bSolver, 0, 0, Fluid, Substance, T(ii), rho_prev, HP_model);
        [mu_next,kappa_next]       = Calculate_HighPressure_Transport_Coeff(bSolver, 0, 0, Fluid, Substance, T(ii), rho_next, HP_model);
        [mu_current,kappa_current] = Calculate_HighPressure_Transport_Coeff(bSolver, 0, 0, Fluid, Substance, T(ii), rho(ii), HP_model);

        dmu_drho(ii)      = (mu_next - mu_prev)./(rho_next - rho_prev);
        d2mu_d2rho(ii)    = (mu_next - 2*mu_current + mu_prev)./(0.5*(rho_next - rho_prev))^2;

        dkappa_drho(ii)   = (kappa_next - kappa_prev)./(rho_next - rho_prev);
        d2kappa_d2rho(ii) = (kappa_next - 2*kappa_current + kappa_prev)./(0.5*(rho_next - rho_prev))^2;


        %% d2mu_drho at 2 different temperatures
        T_1 = T(ii) - T(ii)*stencil_weight;
        T_2 = T(ii) + T(ii)*stencil_weight;

        if strcmp(bSolver,'Real')
            P_1   = Calculate_P_PengRobinson(rho(ii),T_1,Substance);
            P_2   = Calculate_P_PengRobinson(rho(ii),T_2,Substance);
        else
            P_1   = rho(ii)*T_1*Fluid.R_specific;
            P_2   = rho(ii)*T_2*Fluid.R_specific;
        end

        % dmu / drho at T1
        rho_1 = Calculate_Rho_from_TandP( bSolver, T_1,P_1, Fluid, Substance );
        rho_prev = rho_1 - rho_1*stencil_weight;
        rho_next = rho_1 + rho_1*stencil_weight;

        [mu_prev,kappa_prev] = Calculate_HighPressure_Transport_Coeff(bSolver, 0, 0, Fluid, Substance, T_1, rho_prev, HP_model);
        [mu_next,kappa_next] = Calculate_HighPressure_Transport_Coeff(bSolver, 0, 0, Fluid, Substance, T_1, rho_next, HP_model);

        dmu_drho_1(ii)    = (mu_next - mu_prev)./(rho_next - rho_prev);
        dkappa_drho_1(ii) = (kappa_next - kappa_prev)./(rho_next - rho_prev);

        % dmu / drho at T2
        rho_2 = Calculate_Rho_from_TandP( bSolver, T_2,P_2, Fluid, Substance );
        rho_prev = rho_2 - rho_2*stencil_weight;
        rho_next = rho_2 + rho_2*stencil_weight;

        [mu_prev,kappa_prev] = Calculate_HighPressure_Transport_Coeff(bSolver, 0, 0, Fluid, Substance, T_2, rho_prev, HP_model);
        [mu_next,kappa_next] = Calculate_HighPressure_Transport_Coeff(bSolver, 0, 0, Fluid, Substance, T_2, rho_next, HP_model);

        dmu_drho_2(ii)    = (mu_next - mu_prev)./(rho_next - rho_prev);
        dkappa_drho_2(ii) = (kappa_next - kappa_prev)./(rho_next - rho_prev);


        d2mu_drhodT(ii)    = (dmu_drho_2(ii) - dmu_drho_1(ii))./(T_2 - T_1);
        d2kappa_drhodT(ii) = (dkappa_drho_2(ii) - dkappa_drho_1(ii))./(T_2 - T_1);

        %% d2mu_dT at 2 different densities
        rho_1 = rho(ii) - rho(ii)*stencil_weight;
        rho_2 = rho(ii) + rho(ii)*stencil_weight;

        if strcmp(bSolver,'Real')
            P_1   = Calculate_P_PengRobinson(rho_1,T(ii),Substance);
            P_2   = Calculate_P_PengRobinson(rho_2,T(ii),Substance);
        else
            P_1   = rho_1*T(ii)*Fluid.R_specific;
            P_2   = rho_2*T(ii)*Fluid.R_specific;
        end

        % dmu / dT at rho1
        T_1 = Calculate_T_fromPandRho( bSolver, P_1,rho_1,Fluid,Substance);
        T_prev = T_1 - T_1*stencil_weight;
        T_next = T_1 + T_1*stencil_weight;

        [mu_prev,kappa_prev] = Calculate_HighPressure_Transport_Coeff(bSolver, 0, 0, Fluid, Substance, T_prev, rho_1, HP_model);
        [mu_next,kappa_next] = Calculate_HighPressure_Transport_Coeff(bSolver, 0, 0, Fluid, Substance, T_next, rho_1, HP_model);

        dmu_dT_1(ii)      = (mu_next - mu_prev)./(T_next - T_prev);
        dkappa_dT_1(ii)   = (kappa_next - kappa_prev)./(T_next - T_prev);

        % dmu / dT at rho2s
        T_2 = Calculate_T_fromPandRho( bSolver, P_2,rho_2,Fluid,Substance);
        T_prev = T_2 - T_2*stencil_weight;
        T_next = T_2 + T_2*stencil_weight;

        [mu_prev,kappa_prev] = Calculate_HighPressure_Transport_Coeff(bSolver, 0, 0, Fluid, Substance, T_prev, rho_2, HP_model);
        [mu_next,kappa_next] = Calculate_HighPressure_Transport_Coeff(bSolver, 0, 0, Fluid, Substance, T_next, rho_2, HP_model);

        dmu_dT_2(ii)      = (mu_next - mu_prev)./(T_next - T_prev);
        dkappa_dT_2(ii)   = (kappa_next - kappa_prev)./(T_next - T_prev);


        d2mu_dTdrho(ii)    = (dmu_dT_2(ii) - dmu_dT_1(ii))./(rho_2 - rho_1);
        d2kappa_dTdrho(ii) = (dkappa_dT_2(ii) - dkappa_dT_1(ii))./(rho_2 - rho_1);

    else % CoolProp

        %% Partial temperature derivatives
        T_prev = T(ii) - delta_T; %T(ii)*stencil_weight;
        T_next = T(ii) + delta_T; %T(ii)*stencil_weight;

        c_v_prev    = py.CoolProp.CoolProp.PropsSI("CVMASS","T",T_prev,"D",rho(ii),Substance);  	% Specific heat at constant pressure [J/kg/K]
        c_v_next    = py.CoolProp.CoolProp.PropsSI("CVMASS","T",T_next,"D",rho(ii),Substance);  	% Specific heat at constant pressure [J/kg/K]

        dcv_dT(ii) = (c_v_next - c_v_prev)./(T_next - T_prev);

        kappa_prev     = py.CoolProp.CoolProp.PropsSI("L","T",T_prev,"D",rho(ii),Substance);  	% Thermal conductivity at film temperature [W/mK]
        mu_prev        = py.CoolProp.CoolProp.PropsSI("V","T",T_prev,"D",rho(ii),Substance);	 	% Dynamic viscosity at film temperature [kg/ms]
        kappa_next     = py.CoolProp.CoolProp.PropsSI("L","T",T_next,"D",rho(ii),Substance);  	% Thermal conductivity at film temperature [W/mK]
        mu_next        = py.CoolProp.CoolProp.PropsSI("V","T",T_next,"D",rho(ii),Substance);	 	% Dynamic viscosity at film temperature [kg/ms]
        kappa_current  = py.CoolProp.CoolProp.PropsSI("L","T",T(ii),"D",rho(ii),Substance);  	% Thermal conductivity at film temperature [W/mK]
        mu_current     = py.CoolProp.CoolProp.PropsSI("V","T",T(ii),"D",rho(ii),Substance);	 	% Dynamic viscosity at film temperature [kg/ms]
    
        dmu_dT(ii)        = (mu_next - mu_prev)./(T_next - T_prev);
        dkappa_dT(ii)     = (kappa_next - kappa_prev)./(T_next - T_prev);
        d2mu_d2T(ii)      = (mu_next -2*mu_current + mu_prev)./(0.5*(T_next - T_prev)).^2;
        d2kappa_d2T(ii)   = (kappa_next -2*kappa_current + kappa_prev)./(0.5*(T_next - T_prev)).^2;

        %% Partial density derivatives
        rho_prev = rho(ii) - delta_rho; %rho(ii)*stencil_weight;
        rho_next = rho(ii) + delta_rho; %rho(ii)*stencil_weight;

        kappa_prev     = py.CoolProp.CoolProp.PropsSI("L","T",T(ii),"D",rho_prev,Substance);  	% Thermal conductivity at film temperature [W/mK]
        mu_prev        = py.CoolProp.CoolProp.PropsSI("V","T",T(ii),"D",rho_prev,Substance);	 	% Dynamic viscosity at film temperature [kg/ms]
        kappa_next     = py.CoolProp.CoolProp.PropsSI("L","T",T(ii),"D",rho_next,Substance);  	% Thermal conductivity at film temperature [W/mK]
        mu_next        = py.CoolProp.CoolProp.PropsSI("V","T",T(ii),"D",rho_next,Substance);	 	% Dynamic viscosity at film temperature [kg/ms]
        kappa_current  = py.CoolProp.CoolProp.PropsSI("L","T",T(ii),"D",rho(ii),Substance);  	% Thermal conductivity at film temperature [W/mK]
        mu_current     = py.CoolProp.CoolProp.PropsSI("V","T",T(ii),"D",rho(ii),Substance);	 	% Dynamic viscosity at film temperature [kg/ms]
    
        dmu_drho(ii)      = (mu_next - mu_prev)./(rho_next - rho_prev);
        d2mu_d2rho(ii)    = (mu_next - 2*mu_current + mu_prev)./(0.5*(rho_next - rho_prev))^2;

        dkappa_drho(ii)   = (kappa_next - kappa_prev)./(rho_next - rho_prev);
        d2kappa_d2rho(ii) = (kappa_next - 2*kappa_current + kappa_prev)./(0.5*(rho_next - rho_prev))^2;


        %% d2mu_drho at 2 different temperatures
        T_1 = T(ii) - delta_T; %T(ii)*stencil_weight;
        T_2 = T(ii) + delta_T; %T(ii)*stencil_weight;

        % dmu / drho at T1
        P_1   = py.CoolProp.CoolProp.PropsSI("P","T",T(ii),"D",rho(ii),Substance);
        rho_1 = py.CoolProp.CoolProp.PropsSI("D","T",T_1,"P",P_1,Substance);

        rho_prev = rho_1 - delta_rho; %rho_1*stencil_weight;
        rho_next = rho_1 + delta_rho; %rho_1*stencil_weight;

        kappa_prev     = py.CoolProp.CoolProp.PropsSI("L","T",T_1,"D",rho_prev,Substance);  	% Thermal conductivity at film temperature [W/mK]
        mu_prev        = py.CoolProp.CoolProp.PropsSI("V","T",T_1,"D",rho_prev,Substance);	 	% Dynamic viscosity at film temperature [kg/ms]
        kappa_next     = py.CoolProp.CoolProp.PropsSI("L","T",T_1,"D",rho_next,Substance);  	% Thermal conductivity at film temperature [W/mK]
        mu_next        = py.CoolProp.CoolProp.PropsSI("V","T",T_1,"D",rho_next,Substance);	 	% Dynamic viscosity at film temperature [kg/ms]

        dmu_drho_1(ii)    = (mu_next - mu_prev)./(rho_next - rho_prev);
        dkappa_drho_1(ii) = (kappa_next - kappa_prev)./(rho_next - rho_prev);

        % dmu / drho at T2
        P_2   = py.CoolProp.CoolProp.PropsSI("P","T",T(ii),"D",rho(ii),Substance);
        rho_2 = py.CoolProp.CoolProp.PropsSI("D","T",T_2,"P",P_2,Substance);
        rho_prev = rho_2 - delta_rho; %rho_2*stencil_weight;
        rho_next = rho_2 + delta_rho; %rho_2*stencil_weight;

        kappa_prev     = py.CoolProp.CoolProp.PropsSI("L","T",T_2,"D",rho_prev,Substance);  	% Thermal conductivity at film temperature [W/mK]
        mu_prev        = py.CoolProp.CoolProp.PropsSI("V","T",T_2,"D",rho_prev,Substance);	 	% Dynamic viscosity at film temperature [kg/ms]
        kappa_next     = py.CoolProp.CoolProp.PropsSI("L","T",T_2,"D",rho_next,Substance);  	% Thermal conductivity at film temperature [W/mK]
        mu_next        = py.CoolProp.CoolProp.PropsSI("V","T",T_2,"D",rho_next,Substance);	 	% Dynamic viscosity at film temperature [kg/ms]

        dmu_drho_2(ii)    = (mu_next - mu_prev)./(rho_next - rho_prev);
        dkappa_drho_2(ii) = (kappa_next - kappa_prev)./(rho_next - rho_prev);

        d2mu_drhodT(ii)    = (dmu_drho_2(ii) - dmu_drho_1(ii))./(T_2 - T_1);
        d2kappa_drhodT(ii) = (dkappa_drho_2(ii) - dkappa_drho_1(ii))./(T_2 - T_1);

        %% d2mu_dT at 2 different densities
        rho_1 = rho(ii) - delta_rho; %rho(ii)*stencil_weight;
        rho_2 = rho(ii) + delta_rho; %rho(ii)*stencil_weight;

        % dmu / dT at rho1
        P_1 = py.CoolProp.CoolProp.PropsSI("P","T",T(ii),"D",rho(ii),Substance);
        T_1 = py.CoolProp.CoolProp.PropsSI("T","D",rho_1,"P",P_1,Substance);

        T_prev = T_1 - delta_T; %T_1*stencil_weight;
        T_next = T_1 + delta_T; %T_1*stencil_weight;

        kappa_prev     = py.CoolProp.CoolProp.PropsSI("L","T",T_prev,"D",rho_1,Substance);  	% Thermal conductivity at film temperature [W/mK]
        mu_prev        = py.CoolProp.CoolProp.PropsSI("V","T",T_prev,"D",rho_1,Substance);	 	% Dynamic viscosity at film temperature [kg/ms]
        kappa_next     = py.CoolProp.CoolProp.PropsSI("L","T",T_next,"D",rho_1,Substance);  	% Thermal conductivity at film temperature [W/mK]
        mu_next        = py.CoolProp.CoolProp.PropsSI("V","T",T_next,"D",rho_1,Substance);	 	% Dynamic viscosity at film temperature [kg/ms]

        dmu_dT_1(ii)      = (mu_next - mu_prev)./(T_next - T_prev);
        dkappa_dT_1(ii)   = (kappa_next - kappa_prev)./(T_next - T_prev);

        % dmu / dT at rho2s
        P_2 = py.CoolProp.CoolProp.PropsSI("P","T",T(ii),"D",rho(ii),Substance);
        T_2 = py.CoolProp.CoolProp.PropsSI("T","D",rho_2,"P",P_2,Substance);
        
        T_prev = T_2 - delta_T; %T_2*stencil_weight;
        T_next = T_2 + delta_T; %T_2*stencil_weight;

        kappa_prev     = py.CoolProp.CoolProp.PropsSI("L","T",T_prev,"D",rho_2,Substance);  	% Thermal conductivity at film temperature [W/mK]
        mu_prev        = py.CoolProp.CoolProp.PropsSI("V","T",T_prev,"D",rho_2,Substance);	 	% Dynamic viscosity at film temperature [kg/ms]
        kappa_next     = py.CoolProp.CoolProp.PropsSI("L","T",T_next,"D",rho_2,Substance);  	% Thermal conductivity at film temperature [W/mK]
        mu_next        = py.CoolProp.CoolProp.PropsSI("V","T",T_next,"D",rho_2,Substance);	 	% Dynamic viscosity at film temperature [kg/ms]

        dmu_dT_2(ii)      = (mu_next - mu_prev)./(T_next - T_prev);
        dkappa_dT_2(ii)   = (kappa_next - kappa_prev)./(T_next - T_prev);


        d2mu_dTdrho(ii)    = (dmu_dT_2(ii) - dmu_dT_1(ii))./(rho_2 - rho_1);
        d2kappa_dTdrho(ii) = (dkappa_dT_2(ii) - dkappa_dT_1(ii))./(rho_2 - rho_1);





    end
end

% Save in structure
J_t2.dcv_dT         = dcv_dT;
J_t2.dmu_dT         = dmu_dT;
J_t2.dkappa_dT      = dkappa_dT;
J_t2.d2mu_d2T       = d2mu_d2T;
J_t2.d2kappa_d2T    = d2kappa_d2T;
J_t2.dmu_drho       = dmu_drho;
J_t2.d2mu_d2rho     = d2mu_d2rho;
J_t2.dkappa_drho    = dkappa_drho;
J_t2.d2kappa_d2rho  = d2kappa_d2rho;
J_t2.d2mu_drhodT    = d2mu_drhodT;
J_t2.d2kappa_drhodT = d2kappa_drhodT;
J_t2.d2mu_dTdrho    = d2mu_dTdrho;
J_t2.d2kappa_dTdrho = d2kappa_dTdrho;


end
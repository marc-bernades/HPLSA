function J_t = Jacobian_thermodynamics(bSolver, rho,T,P, J_t2,Substance, HP_model)

% Peng-Robinson
[ a,b,R,dadT,d2adT2,NASA_coefficients ] = PengRobinson( T, Substance );

% N2 substance library
[MW, Tc, pc, p_inf, rhoc, vc_bar, omega, gamma, e_0, c_v, NASA_coefficients, ...
    mu_0, kappa_0, T_0, S_mu, S_kappa, dipole_moment, association_factor] = Substance_Library(Substance);

% Alpha calculation (part of a coefficient)
alpha   = a./(0.457236*(R*Tc)^2/pc);
a       = a./alpha;

% Peng Robsinson coefficients
if omega > 0.49 % Accentric factor
    c      = 0.379642 + 1.48503*omega - 0.164423*omega^2 + 0.016666*omega^3;
else
    c      = 0.37464 + 1.54226*omega - 0.26992*omega^2;
end

if strcmp(bSolver,'Real')

    if ~strcmp(HP_model,'CoolProp') && ~strcmp(HP_model,'RefProp')

        % Jacobians
        J_t.dP_dT       = rho*R./(1 - rho.*b) + c*sqrt(alpha./(T*Tc)).*(a.*rho.^2)./(1 + 2*b*rho - b.^2.*rho.^2);

        J_t.dP_drho     = R*T./(1 - rho*b).^2 - alpha.*(2*a.*rho + 2*a.*b.*rho.^2)./(1 + 2*b.*rho - b.^2.*rho.^2).^2;

        J_t.d2P_drhodT  = R./(1 - rho.*b).^2 + c*sqrt(alpha./(T*Tc)).*(2*a.*rho + 2*a.*b.*rho.^2)./(1 + 2*b.*rho - b.^2.*rho.^2).^2;

        J_t.d2P_d2T     = - c.*(1+c)./(2*sqrt(T.^3*Tc)).*a.*rho.^2./(1 + 2*b.*rho - b.^2.*rho.^2);

        J_t.d2P_d2rho   = 2*R*b.*T./(1-rho.*b).^3 - alpha.*2.*a.*(2*b.^3*rho.^3 + 3*b.^2*rho.^2 + 1)./(1 + 2*b.*rho - b.^2*rho.^2).^3;

        J_t.de_dT       = c_v + J_t2.dcv_dT.*T + a./(4.*sqrt(2).*b).*(-c.*(1+c).*sqrt(1./(T.*Tc)).*log((1+b.*(1 - sqrt(2)).*rho)./(1+b.*(1 + sqrt(2)).*rho)));

        J_t.de_drho     = - a./(1 + 2*b.*rho - b^2*rho.^2).*((1 + c)^2 - c*(1+c).*sqrt(T./Tc));


    else
        if strcmp(HP_model,'RefProp')
            Substance = strcat('REFPROP::',Substance);
        end

        for nCP = 1:length(rho)

            % Jacobians
            J_t.dP_dT(nCP,1)       = py.CoolProp.CoolProp.PropsSI("d(P)/d(T)|D","T",T(nCP),"P",P(nCP),Substance);

            J_t.dP_drho(nCP,1)     = py.CoolProp.CoolProp.PropsSI("d(P)/d(D)|T","T",T(nCP),"P",P(nCP),Substance);

            J_t.d2P_drhodT(nCP,1)  = py.CoolProp.CoolProp.PropsSI("d(d(P)/d(T)|D)/d(D)|T","T",T(nCP),"P",P(nCP),Substance);

%             J_t.d2P_d2T(nCP,1)     = py.CoolProp.CoolProp.PropsSI("d(d(P)/d(T)|D)/d(T)|D","T",T(nCP),"P",P(nCP),'CO2');
            J_t.d2P_d2T(nCP,1)     = py.CoolProp.CoolProp.PropsSI("d(d(P)/d(T)|D)/d(T)|D","T",T(nCP),"P",P(nCP),Substance);

%             J_t.d2P_d2rho(nCP,1)   = py.CoolProp.CoolProp.PropsSI("d(d(P)/d(D)|T)/d(D)|T","T",T(nCP),"P",P(nCP),'CO2');
            J_t.d2P_d2rho(nCP,1)   = py.CoolProp.CoolProp.PropsSI("d(d(P)/d(D)|T)/d(D)|T","T",T(nCP),"P",P(nCP),Substance);


            J_t.de_dT(nCP,1)       = py.CoolProp.CoolProp.PropsSI("d(U)/d(T)|D","T",T(nCP),"P",P(nCP),Substance);

            J_t.de_drho(nCP,1)     = py.CoolProp.CoolProp.PropsSI("d(U)/d(D)|T","T",T(nCP),"P",P(nCP),Substance);

        end

%         % Plot check
%         Name_var = fieldnames(J_t);
%         figure
%         for ii = 1:length(Name_var)
%             subplot(length(Name_var),1,ii)
%             plot(T,J_t.(Name_var{ii}));
%         end

    end
else
    % Jacobians Ideal gas
    J_t.dP_dT       = zeros(size(rho)) + rho*R;

    J_t.dP_drho     = zeros(size(rho)) + R*T;

    J_t.d2P_drhodT  = zeros(size(rho)) + R;

    J_t.d2P_d2T     = zeros(size(rho)) + 0;

    J_t.d2P_d2rho   = zeros(size(rho)) + 0;

    J_t.de_dT       = zeros(size(rho)) + c_v;

    J_t.de_drho     = zeros(size(rho));
end

end
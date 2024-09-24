function [Re,Pr,Ec,Ma] = DimensionlessNumbers(rho_b,u_b,T_b,mu_b,kappa_b,c_p_b,sos_b,rho,u,T,mu,kappa,c_p,sos,delta,du_dy,bScaling, varargin)

if strcmp(bScaling,'bulk')

    % Normalization metrics wall/center

    if isempty(varargin)
        [~,idx] = min(abs(du_dy-0)); % Only velocity at center
        u_c = u(idx);
    else
        u_c = varargin{1};
    end
    

    Re = rho_b.*u_c*(delta)./mu_b;
    Pr = mu_b.*c_p_b/kappa_b;
    Ec = u_c^2./(c_p_b*T_b);
    Ma = u_c/sos_b;

elseif strcmp(bScaling,'wall')

    % Normalization metrics wall/center

    if isempty(varargin)
        [~,idx] = min(abs(du_dy-0)); % Only velocity at center
        u_c = u(idx);
    else
        u_c = varargin{1};
    end
    Re = rho(1).*u_c*delta./mu(1);
    Pr = mu(1).*c_p(1)/kappa(1);
    Ec = u_c^2./(c_p(1)*T(1));
    Ma = u_c/sos(1);

else

    Re = rho_b.*u_b*(2*delta)./mu_b;
    Pr = mu_b.*c_p_b/kappa_b;
    Ec = u_b^2./(c_p_b*T_b);
    Ma = u_b/sos_b;
    
end


end
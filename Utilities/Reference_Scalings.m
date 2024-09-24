function norm = Reference_Scalings(u_b,T_b,rho_b,mu_b,kappa_b,c_p_b,c_v_b,E_b,e_b,P_b,u,T,rho,mu,kappa,c_p,c_v,lambda,Dy,delta,bScaling, varargin)

if strcmp(bScaling,'bulk')
    % Normalization metrics bulk
    % Normalization metrsics wall/center(max)
    if isempty(varargin)
        [~,idx] = min(abs(Dy.du_dy-0));
        u_c = u(idx);
    else
        u_c = varargin{1};
    end
    
    norm.u      = u_c;
    norm.T      = T_b;
    norm.rho    = rho_b;
    norm.mu     = mu_b;
    norm.kappa  = kappa_b;
    norm.E      = norm.u^2;
    norm.e      = norm.u^2;
    norm.P      = norm.rho*norm.u^2;
    norm.c_p    = c_p_b;
    norm.c_v    = c_v_b;
    norm.lambda = norm.mu; %-2/3*norm.mu;
    norm.y      = delta;

%     norm.E      = E_b;
%     norm.e      = e_b;
%     norm.P      = P_b;
%     norm.c_p    = c_p_b;
%     norm.c_v    = c_v_b;
%     norm.lambda = -2/3*norm.mu; %norm.mu
%     norm.y      = delta;




elseif strcmp(bScaling,'wall')
    % Normalization metrsics wall/center(max)
    if isempty(varargin)
        [~,idx] = min(abs(Dy.du_dy-0));
        u_c = u(idx);
    else
        u_c = varargin{1};
    end
    norm.u      = u_c;
    norm.T      = T(1);
    norm.rho    = rho(1);
    norm.mu     = mu(1);
    norm.kappa  = kappa(1);
    norm.E      = norm.u^2;
    norm.e      = norm.u^2;
    norm.P      = norm.rho*norm.u^2;
    norm.c_p    = c_p(1); %norm.kappa/norm.mu; %c_p(1);
    norm.c_v    = c_v(1); %norm.kappa/norm.mu; %c_v(1);
    norm.lambda = mu(1);  %lambda(1); % For LST modal-analysis: lambda(1);
    norm.y      = delta;

end

end
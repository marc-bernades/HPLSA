function dx = CentralDerivative_d2_2ndOrder_1D(u,y)


% Allocate memory for derivatives
dx = zeros(1,length(u));

% Extremes upwind first order (Boundary points not needed)
% dx(1)    = (2*u(3)   - 2*u(2)     + u(1))./(0.5*(y(3) - y(1)))^2;           %Forward
% dx(end)  = (2*u(end) - 2*u(end-1) + u(end-2))./(0.5*(y(end) - y(end-2)))^2; %Backward

dx(1)    = (u(3)   - 2*u(2)     + u(1))./(0.5*(y(3) - y(1))).^2;           %Forward
dx(end)  = (u(end) - 2*u(end-1) + u(end-2))./(0.5*(y(end) - y(end-2))).^2; %Backward


% Sweep internal points
for i = 2:(length(u)-1)
    dx(i)    = (u(i+1) - 2*u(i) + u(i-1))./(0.5*(y(i+1) - y(i-1)))^2;
end


end
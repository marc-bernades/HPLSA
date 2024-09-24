function dx = CentralDerivative_d1_2ndOrder_1D(u,y)


% Allocate memory for derivatives
dx = zeros(1,length(u));

% Extremes upwind first order (Boundary points not needed)
dx(1)    = (u(2) - u(1))./(y(2) - y(1));
dx(end)  = (u(end) - u(end - 1))./(y(end) - y(end - 1));

% dx(1)    = (-u(3)   + 4*u(2)     - 3*u(1))./((y(3) - y(1)));           %Forward
% dx(end)  = (3*u(end) - 4*u(end-1) + u(end-2))./((y(end) - y(end-2))); %Backward


% Sweep internal points
for i = 2:(length(u)-1)
    dx(i)    = (u(i+1) - u(i-1))./(y(i+1) - y(i-1));
    %     dx(i)    = (u(i+1) - u(i))./(y(i+1) - y(i));
end


end
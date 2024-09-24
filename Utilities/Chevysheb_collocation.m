function [y,D] = Chevysheb_collocation(N,delta)

% Reference Spectral methods in MATLAB"
% SIAM Publications LL.N. Trefethen
y = delta*(1 - cos (pi*(0:N)/N)');
c =[2;ones(N-1,1);2].*(-1).^(0:N)';
Y =repmat(y,1,N+1);
dY = Y-Y';
D = (c*(1./c)')./(dY+eye(N+1)); % off-diagonal entries
D = D-diag(sum(D'));


end
function df = BaseRootFindingMinimization(x)

%     /// Globally Convergent Methods for Nonlinear Systems of Equations: fdjac method
%     /// W. H. Press, S. A. Teukolsky, W. T. Vetterling, B. P. Flannery.
%     /// Numerical recipes in C++.
%     /// Cambridge University Press, 2001.

%     //const double EPS = 1.0e-8;
EPS = 1.0e-6;
n = size(x);

for j = 1:n
    temp = x(j);
    h = EPS*abs( temp );
    if( h == 0.0 )
        h = EPS;
    end
    x(j) = temp + h;
    h = x(j) - temp;
    x(j) = temp;
    for i = 1:n
        df(i,j) = ( f(i) - fvec(i) )/h;

    end
end
end


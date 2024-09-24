function P = Calculate_P_PengRobinson(rho,T,Substance)


[ a,b,R,~,~,~ ] = PengRobinson( T, Substance );

v = 1./rho;

P = R.*T./(v - b) - a./(v.^2 + 2*b*v - b.^2);


end
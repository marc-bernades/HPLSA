function [ GG, t_all, GG_max, q_in, q_out, Val_eigen_tg, Vec_eigen_tg] = LST_TransientGrowth( N,y, q_vars, Val_eigen, Vec_eigen, ...
    Imag_max, Imag_min, Real_max, Real_min, N_t, t_max, m_d, m_T , N_modes)


% Reduce eigenvalue and eigenvector
% Val_eigen_tg = zeros(N*q_vars,1);
% Vec_eigen_tg = zeros(N*q_vars,N*q_vars);
Val_eigen_tg = zeros(N_modes,1);
Vec_eigen_tg = zeros(N*q_vars,N_modes);


num     = 0;
num_INF = 0;
for i = 1:N_modes %N*q_vars
    if(isinf(real(Val_eigen(i)))||isinf(imag(Val_eigen(i))))
        num_INF=num_INF+1;
    else
        if(imag(Val_eigen(i))>Imag_min...
                && imag(Val_eigen(i))<Imag_max...
                && real(Val_eigen(i))<Real_max...
                && real(Val_eigen(i))>Real_min)
            num = num+1;
            Val_eigen_tg(num)   = Val_eigen(i);
            Vec_eigen_tg(:,num) = Vec_eigen(:,i);
        end
    end
end

% Val_eigen_tg(num+1:N*q_vars)   = [];
% Vec_eigen_tg(:,num+1:N*q_vars) = [];

Val_eigen_tg(num+1:N_modes)   = [];
Vec_eigen_tg(:,num+1:N_modes) = [];

N_tg = num;
fprintf('Transient growth modes = %3d, ', N_tg);
fprintf('Infinit modes modes    = %3d, ', num_INF);
fprintf('\n')


% Differentiate y-vector
Y=zeros(N,1);
for i=2:N-1
    Y(i)=0.5*(y(i+1)-y(i-1));
end
Y(1)=0.5*y(2);
Y(N)=0.5*(y(N)-y(N-1));

% Weight matrix
if length(m_d) <= 1
    M  = sqrt(diag([m_d, 1, 1, 1, m_T]));
    % X*y all components (rho,u,v,w,T) for each grid point
    MM   = kron(sqrt(diag(Y)),M);

else
    for jj = 1:length(m_d)

        M{jj}  = sqrt(Y(jj))*sqrt(diag([m_d(jj), 1, 1, 1, m_T(jj)]));

    end
    MM  = blkdiag(M{:});
end


% Operator for SVD
work = MM*Vec_eigen_tg;
A_tg = work'*work;

% Transient growth
img   = sqrt(-1);
t_all = linspace(0,t_max,N_t);
GG    = zeros(N_t,1);

[U,S,V] = svd(A_tg);
s       = sqrt(diag(S));
F       = diag(s)*U';
F_inv   = U*diag(ones(size(s))./s);

for i = 1:N_t
    t = t_all(i);

    lambda=zeros(N_tg,1);
    for j = 1:N_tg
        lambda(j) = exp(-img*Val_eigen_tg(j)*t);
    end
    LAMDA = diag(lambda);
    clear lambda

    evol  = F*LAMDA*F_inv;
    G     = norm(evol)^2;
    GG(i) = G;
end

% Optimal Perturbation
[GG_max, idx_max] = max(GG);

t      = t_all(idx_max);

% if GG_max > 1E4 && t < 10
%     GG_max = NaN;
% end

lambda = zeros(N_tg,1);
for j = 1:N_tg
    lambda(j) = exp(-img*Val_eigen_tg(j)*t);
end
LAMDA   = diag(lambda);
clear lambda
evol    = F*LAMDA*F_inv;
[U,S,V] = svd(evol);
q_in    = Vec_eigen_tg*F_inv*V(:,1);
q_out   = Vec_eigen_tg*F_inv*U(:,1);





end
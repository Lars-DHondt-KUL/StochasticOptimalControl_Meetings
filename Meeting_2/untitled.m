
import casadi.*
% mex f_integration.c -largeArrayDims  % Matlab

%%

% x_mean = zeros(2,N+1);
% F = zeros(1,N+1);
% Pvar = zeros(4,N+1);
% K = 0;
% D = 0;
% sigma_w = 10^(-2);

x_mean_MX = MX.sym('x_mean',2,N+1);
F_MX = MX.sym('F',1,N+1);
Pvar_MX = MX.sym('Pvar',4,N+1);
K_MX = MX.sym('K',1,1);
D_MX = MX.sym('D',1,1);
sigma_w_MX = MX.sym('sigma_w',1,1);

%%

C = Importer('f_integration.c','clang');
f_integration = external('map60_f_integration_k',C);

[obj_N,eq_constr_N] = f_integration(x_mean_MX(:,k), x_mean_MX(:,k+1), F_MX(:,k), F_MX(:,k+1),...
    Pvar_MX(:,k), Pvar_MX(:,k+1), K_MX, D_MX, sigma_w_MX);

%%



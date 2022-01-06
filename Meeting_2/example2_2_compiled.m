
% this one does not work




% function [failed] = example2_2c(make_fig)

clear all; clc; close all;
make_fig = 1;
addpath([pwd '\Integrator']);


import casadi.*
dt = 0.167; T = 10.02; t = 0:dt:T;
N = round(T/dt);
k = 1:N;


%- Stochastic Dynamics of the controlled system
% Variables initialisaiton
x_mean = SX.sym('x_mean',2);
w =  SX.sym('w',1); % disturbing force
F =  SX.sym('F',1);
K = SX.sym('X',1,1);
B = SX.sym('X',1,1);
% Dynamic System parameters
g = 9.81; l = 1.2; m = 5;

% Dynamics of stochastic problem ( xdot = f(x,u,...) )
x_mean_dot = [ x_mean(2); (-m*g*l*sin(x_mean(1)) + F + w + K*(x_mean(1)-pi) + B*x_mean(2))/(m*l^2)];
fcn_dyn = Function('fcn_dyn',{x_mean,F,w,K,B},{x_mean_dot});
% Linearization of dynamics around mean state (note that fcn_A is dependent on the feedback gains as well!!)
fcn_A = Function('fcn_A',{x_mean,F,w,K,B},{jacobian(x_mean_dot,x_mean)});
fcn_C = Function('fcn_C',{x_mean,F,w,K,B},{jacobian(x_mean_dot,w)});

clear x_mean F K B w

%%

x_mean_k_SX = SX.sym('x_mean_k',2); % x_mean(:,i);
x_mean_k1_SX = SX.sym('x_mean_k1',2); % x_mean(:,i+1);
F_k_SX = SX.sym('F_k',1); % F(i);
F_k1_SX = SX.sym('F_k1',1); % F(i+1);
Pvar_k_SX = SX.sym('Pvar_k',4); % Pvar(:,i);
Pvar_k1_SX = SX.sym('Pvar_k1',4); % Pvar(:,i+1);
K_SX = SX.sym('K',1,1);
D_SX = SX.sym('D',1,1);
sigma_w_SX = SX.sym('sigma_w',1,1);

P_var_k = reshape(Pvar_k_SX,2,2);
P_var_k1 = reshape(Pvar_k1_SX,2,2);

% integration of mean state trajectory
x_mean_dot = fcn_dyn(x_mean_k_SX,F_k_SX,0,K_SX,D_SX);
x_mean_dot_next = fcn_dyn(x_mean_k1_SX,F_k1_SX,0,K_SX,D_SX);
eq_constr_1_k = G_Trapezoidal(x_mean_k_SX,x_mean_k1_SX,x_mean_dot,x_mean_dot_next,dt);

% evaluate linearization of dynamics
A_k = fcn_A(x_mean_k_SX,F_k_SX,0,K_SX,D_SX);
A_k_next = fcn_A(x_mean_k1_SX,F_k1_SX,0,K_SX,D_SX);
C_k = fcn_C(x_mean_k_SX,F_k_SX,0,K_SX,D_SX);

% integration of covariance matrix
% - positive definiteness preserving Lyapunov discretisation (ignore for now)
DG_DW = DG_DW_Trapezoidal(C_k,dt);
DG_DX = DG_DX_Trapezoidal(A_k,dt);
DG_DZ = DG_DZ_Trapezoidal(A_k_next,dt);
M_k = (DG_DZ)^(-1);
% - integration step
P_k = M_k*(DG_DX*P_var_k*DG_DX' + DG_DW*sigma_w_SX*DG_DW')*M_k';
eq_constr_2_k = reshape((P_var_k1 - P_k),4,1);

% expected effort
obj_k1 =(F_k1_SX + K_SX*(x_mean_k1_SX(1)-pi) + D_SX*x_mean_k1_SX(2))^2 + [K_SX D_SX]*P_k*[K_SX;D_SX];

f_integration_k = Function('f_integration_k',{x_mean_k_SX,x_mean_k1_SX,F_k_SX,F_k1_SX,...
    Pvar_k_SX,Pvar_k1_SX,K_SX,D_SX,sigma_w_SX},{obj_k1,[eq_constr_1_k;eq_constr_2_k]});

f_integration_map = f_integration_k.map(N);


%%
% f_integration_map.generate('f_integration.c');
% 
% 
% C = Importer('f_integration.c','clang');
% f_external = external('map60_f_integration_k',C);

%%
C = CodeGenerator('f_integration.c');
C.add(f_integration_map);
C.add(f_integration_map.jacobian());
C.generate();

%%

CI = Importer('f_integration.c','clang');
% f_external = external('map60_f_integration_k',CI);
% f_integration_N = external('map60_f_integration_k',CI);
f_integration_N = external('f_integration',CI);

%%
% x_mean_MX = MX.sym('x_mean',2,N+1);
% F_MX = MX.sym('F',1,N+1);
% Pvar_MX = MX.sym('Pvar',4,N+1);
% K_MX = MX.sym('K',1,1);
% D_MX = MX.sym('D',1,1);
% sigma_w_MX = MX.sym('sigma_w',1,1);
% 
% [obj_N_MX,eq_constr_N_MX] = f_external(x_mean_MX(:,k), x_mean_MX(:,k+1), F_MX(:,k), F_MX(:,k+1),...
%     Pvar_MX(:,k), Pvar_MX(:,k+1), K_MX, D_MX, sigma_w_MX);
% 
% f_integration_N = Function('f_integration_N',{x_mean_MX,F_MX,Pvar_MX,K_MX,D_MX,sigma_w_MX},...
%     {sum(obj_N_MX),vertcat(eq_constr_N_MX)});
% 
% f_integration_N.generate('f_integration_N.c');
%%

opti = casadi.Opti(); % Create opti instance

noise_param = opti.parameter();
sigma_w = 10^(noise_param-3);

F = opti.variable(1,N+1);      % force/torque
x_mean = opti.variable(2,N+1); % mean trajectory
Pvar = opti.variable(4,N+1);   % covariance matrix
K = opti.variable(1,1);        % stiffness
D = opti.variable(1,1);        % damping
opti.subject_to(K == 0);
opti.subject_to(D == 0);


P_1 = 0.0001*eye(2);
opti.set_initial(Pvar,0.0001*ones(4,N+1));
opti.subject_to(Pvar(:,1) == P_1(:));

% boundary conditions mean state
opti.subject_to(x_mean(:,end) == [pi; 0]);
opti.subject_to(x_mean(:,1) == 0);

[obj_N,eq_constr_N] = f_integration_N(x_mean(:,k), x_mean(:,k+1), F(:,k), F(:,k+1),...
    Pvar(:,k), Pvar(:,k+1), K, D, sigma_w);

% [obj_N,eq_constr_N] = f_integration_N(x_mean, F, Pvar, K, D, sigma_w);

opti.subject_to(eq_constr_N == 0);

% expected effort
obj_1 = (F(:,1) + K*(x_mean(1,1)-pi) + D*x_mean(2,1))^2 + [K D]*P_1*[K;D];
obj = obj_1 + sum(obj_N);

% cost function
% Add minimization of variance, otherwise feedback gains will be zero and feedforward only is optimal solution
opti.minimize( norm(obj + 1e4*(Pvar(1,N+1))) );

% solver options
optionssol.ipopt.linear_solver = 'mumps';
optionssol.ipopt.tol = 1e-5; 
optionssol.ipopt.dual_inf_tol = 1e-5;
optionssol.ipopt.constr_viol_tol = 1e-7;

opti.solver('ipopt',optionssol);

%%
% Plot variables
if make_fig
    figure();
    cs = linspecer(8,'blue');
    cs = cs(2:end,:);
end

% OCP for seven different noise levels
failed = zeros(1,7);

for noise_factor = 1:7

%     try
        opti.set_value(noise_param,noise_factor);
        sol = opti.solve();
    
        x_mean_sol = sol.value(x_mean);
        F_sol = sol.value(F);
        P_sol = sol.value(Pvar);
        K_sol = sol.value(K);
        D_sol = sol.value(D);
        F_sol = F_sol + K_sol*(x_mean_sol(1,:)-pi) + D_sol*x_mean_sol(2,:);
        
        if make_fig
            subplot(3,2,1)
            plot(t,x_mean_sol(1,:),'Color',cs(noise_factor,:));
            title('position');
            hold on;
            subplot(3,2,2)
            plot(t,180/pi*sqrt(P_sol(1,:)),'Color',cs(noise_factor,:));
            title('position SD');
            hold on;
            set(gca, 'YScale', 'log')
            
            subplot(3,2,3)
            plot(t,x_mean_sol(2,:),'Color',cs(noise_factor,:));
            title('velocity');
            hold on;
            
            subplot(3,2,4)
            plot(t,180/pi*sqrt(P_sol(4,:)),'Color',cs(noise_factor,:));
            title('velocity SD');
            hold on;
            set(gca, 'YScale', 'log')
            
            subplot(3,2,5)
            plot(t,F_sol(1,:),'Color',cs(noise_factor,:));
            title('torque');
            hold on;
        
            subplot(3,2,6)
            plot(noise_factor,'d','Color',cs(noise_factor,:));
            title('noise factor');
            hold on;
        end
%     catch
%         failed(noise_factor) = 1;
%     end


end
clear
clc

addpath([pwd '\Integrator']);

%% original code


import casadi.*
dt = 0.167; 
T = 10.02; t = 0:dt:T;
N = round(T/dt);



%- Stochastic Dynamics of the controlled system
% Variables initialisaiton
x_mean = MX.sym('x_mean',2);
w = MX.sym('w',1); % disturbing force
F = MX.sym('F',1);
K = MX.sym('X',1,1);
B = MX.sym('X',1,1);
% Dynamic System parameters
g = 9.81; l = 1.2; m = 5;

% Dynamics of stochastic problem ( xdot = f(x,u,...) )
x_mean_dot = [ x_mean(2); (-m*g*l*sin(x_mean(1)) + F + w + K*(x_mean(1)-pi) + B*x_mean(2))/(m*l^2)];
fcn_dyn = Function('fcn_dyn',{x_mean,F,w,K,B},{x_mean_dot});
% Linearization of dynamics around mean state (note that fcn_A is dependent on the feedback gains as well!!)
fcn_A = Function('fcn_A',{x_mean,F,w,K,B},{jacobian(x_mean_dot,x_mean)});
fcn_C = Function('fcn_C',{x_mean,F,w,K,B},{jacobian(x_mean_dot,w)});

clear x_mean F K B w


% OCP for seven different noise levels
noise_factor = 1;
sigma_w = 10^(noise_factor-3);

F = MX.sym('F',1,N+1);      % force/torque
x_mean = MX.sym('x_mean',2,N+1); % mean trajectory
Pvar = MX.sym('Pvar',4,N+1);   % covariance matrix
K = MX.sym('K',1,1);        % stiffness
D = MX.sym('D',1,1);        % damping

eq_constr = {};
obj = 0;

for i = 1:N
    % integration of mean state trajectory
    x_mean_dot = fcn_dyn(x_mean(:,i),F(i),0,K,D);
    x_mean_dot_next = fcn_dyn(x_mean(:,i+1),F(i+1),0,K,D);
    eq_constr{end+1} = G_Trapezoidal(x_mean(:,i),x_mean(:,i+1),x_mean_dot,x_mean_dot_next,dt);
    
    % evaluate linearization of dynamics
    A_k = fcn_A(x_mean(:,i),F(i),0,K,D);
    A_k_next = fcn_A(x_mean(:,i+1),F(i+1),0,K,D);
    C_k = fcn_C(x_mean(:,i),F(i),0,K,D);

    % integration of covariance matrix
    % - positive definiteness preserving Lyapunov discretisation (ignore for now)
    DG_DW = DG_DW_Trapezoidal(C_k,dt);
    DG_DX = DG_DX_Trapezoidal(A_k,dt);
    DG_DZ = DG_DZ_Trapezoidal(A_k_next,dt);
    M_k = (DG_DZ)^(-1);
    % - integration step
    P_k = M_k*(DG_DX*reshape(Pvar(:,i),2,2)*DG_DX' + DG_DW*sigma_w*DG_DW')*M_k';
    eq_constr{end+1} = Pvar(:,i+1) - P_k(:);
    
    % expected effort
    obj = obj + (F(:,i+1) + K*(x_mean(1,i+1)-pi) + D*x_mean(2,i+1))^2 + [K D]*P_k*[K;D]; % expected effort
end

eq_constr = vertcat(eq_constr{:});

f_OCP_original = Function('f_OCP_original',{x_mean,F,Pvar,K,D,},{obj,eq_constr});

%% adapted code
clearvars -except f_OCP_original N

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

%

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

%
noise_factor = 1;
sigma_w = 10^(noise_factor-3);

F = MX.sym('F',1,N+1);      % force/torque
x_mean = MX.sym('x_mean',2,N+1); % mean trajectory
Pvar = MX.sym('Pvar',4,N+1);   % covariance matrix
K = MX.sym('K',1,1);        % stiffness
D = MX.sym('D',1,1);        % damping

[obj_N,eq_constr_N] = f_integration_map(x_mean(:,k), x_mean(:,k+1), F(:,k), F(:,k+1),...
    Pvar(:,k), Pvar(:,k+1), K, D, sigma_w);

obj = sum(obj_N);
eq_constr = vertcat(eq_constr_N(:));

f_OCP_adapted = Function('f_OCP_adapted',{x_mean,F,Pvar,K,D,},{obj,eq_constr});


%% compare them
clearvars -except f_OCP_original f_OCP_adapted N

x_mean = lhsdesign(2,N+1);
F = lhsdesign(1,N+1);
Pvar = lhsdesign(4,N+1);
K = 0;
D = 0;


[obj1,eq1] = f_OCP_original(x_mean,F,Pvar,K,D);

[obj2,eq2] = f_OCP_adapted(x_mean,F,Pvar,K,D);

diff_obj = obj1 - obj2;
diff_eq = eq1 - eq2;

disp(full(diff_obj))
disp(rms(full(diff_eq)))

%% generate c code
import casadi.*

C = CodeGenerator('f_OCP_original.c');
C.add(f_OCP_original);
% C.add(f_OCP_original.jacobian());
C.generate();
clearvars C

C = CodeGenerator('f_OCP_adapted.c');
C.add(f_OCP_adapted);
% C.add(f_OCP_adapted.jacobian());
C.generate();
clearvars C





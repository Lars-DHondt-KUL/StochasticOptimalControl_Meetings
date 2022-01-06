function [Results] = example2_2c(make_fig,noise_factors,dt,fb)

% clear all; clc; close all;
% make_fig = 0;
% noise_factors = [1,3,5];
% dt = 0.167;
% fb = 0;

import casadi.*

T = 10.02; t = 0:dt:T;
N = round(T/dt);
k = 1:N;


%- Stochastic Dynamics of the controlled system
% Variables initialisaiton
x_mean = SX.sym('x_mean',2);
w = SX.sym('w',1); % disturbing force
F = SX.sym('F',1);
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

%% Generate expression graph for integration over 1 mesh interval

x_mean_k_SX = SX.sym('x_mean_k',2);     % x_mean(:,i);
x_mean_k1_SX = SX.sym('x_mean_k1',2);   % x_mean(:,i+1);
F_k_SX = SX.sym('F_k',1);               % F(i);
F_k1_SX = SX.sym('F_k1',1);             % F(i+1);
Pvar_k_SX = SX.sym('Pvar_k',4);         % Pvar(:,i);
Pvar_k1_SX = SX.sym('Pvar_k1',4);       % Pvar(:,i+1);
K_SX = SX.sym('K',1,1);
D_SX = SX.sym('D',1,1);
sigma_w_SX = SX.sym('sigma_w',1,1);

P_var_k_SX = reshape(Pvar_k_SX,2,2);    % Pvar as a square matrix
P_var_k1_SX = reshape(Pvar_k1_SX,2,2);  % Pvar as a square matrix

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
P_k = M_k*(DG_DX*P_var_k_SX*DG_DX' + DG_DW*sigma_w_SX*DG_DW')*M_k';
eq_constr_2_k = reshape((P_var_k1_SX - P_k),4,1);

% expected effort
obj_k1 =(F_k1_SX + K_SX*(x_mean_k1_SX(1)-pi) + D_SX*x_mean_k1_SX(2))^2 + [K_SX D_SX]*P_k*[K_SX;D_SX];

% casadi function with expression graph
f_integration_k = Function('f_integration_k',{x_mean_k_SX,x_mean_k1_SX,F_k_SX,F_k1_SX,...
    Pvar_k_SX,Pvar_k1_SX,K_SX,D_SX,sigma_w_SX},{obj_k1,[eq_constr_1_k;eq_constr_2_k]});

%% Map the expression graph over N mesh intervals
f_integration_map = f_integration_k.map(N);

% f_integration_N.generate('f_integration_map.c');

%% Create proxy function to reduce the inputs
% proxy variables for opti variables
F_SX = SX.sym('F',1,N+1);           % force/torque
x_mean_SX = SX.sym('x_mean',2,N+1); % mean trajectory
Pvar_SX = SX.sym('Pvar',4,N+1);     % covariance matrix
K_SX = SX.sym('K',1,1);             % stiffness
D_SX = SX.sym('D',1,1);             % damping
sigma_w_SX = SX.sym('sigma_w',1,1); % noise level

P_1_SX = reshape(Pvar_SX(:,1),2,2); % covariance matrix at initial mesh point

% % boundary conditions mean state
% eq_constr_final_SX = x_mean_SX(:,end) - [pi; 0];
% eq_constr_init_SX = x_mean_SX(:,1);
% 
% % P at 1st mesh point        
% eq_constr_init_P_SX = P_1_SX - 0.0001*eye(2);

% expected effort
obj_1_SX = (F_SX(:,1) + K_SX*(x_mean_SX(1,1)-pi) + D_SX*x_mean_SX(2,1))^2 + [K_SX D_SX]*P_1_SX*[K_SX;D_SX];

% Add minimization of variance, otherwise feedback gains will be zero and feedforward only is optimal solution
obj_P_SX = Pvar_SX(1,N+1) *1e4;

% Evaluate expression graph 
[obj_SX,eq_constr_SX] = f_integration_map(x_mean_SX(:,k), x_mean_SX(:,k+1), F_SX(:,k), F_SX(:,k+1),...
    Pvar_SX(:,k), Pvar_SX(:,k+1), K_SX, D_SX, sigma_w_SX);

% compose objective and constraint functions
obj_SX = norm(obj_1_SX + sum(obj_SX) + obj_P_SX);
% eq_constr_SX = vertcat(eq_constr_SX(:),eq_constr_final_SX(:),eq_constr_init_SX(:),eq_constr_init_P_SX(:));
eq_constr_SX = vertcat(eq_constr_SX(:));

% define OCP casadi function
f_OCP = Function('f_OCP',{x_mean_SX,F_SX,Pvar_SX,K_SX,D_SX,sigma_w_SX},{obj_SX,eq_constr_SX});

%%

% C = CodeGenerator('f_OCP_adapted2.c');
% C.add(f_OCP);
% C.generate();

%% Formulate optimisation problem
opti = casadi.Opti();           % Create opti instance

% variables
% force/torque
F = opti.variable(1,N+1);

% mean trajectory
x_mean = opti.variable(2,N+1);

% covariance matrix
Pvar = opti.variable(4,N+1);
opti.set_initial(Pvar,0.0001*ones(4,N+1));

% feedback gains
K = opti.variable(1,1);         % stiffness
D = opti.variable(1,1);         % damping
if ~fb
    opti.subject_to(K == 0);
    opti.subject_to(D == 0);
end


% boundary conditions mean state
opti.subject_to(x_mean(:,end) == [pi; 0]);
opti.subject_to(x_mean(:,1) == 0);
opti.subject_to(Pvar(:,1) == 0.0001);

% parameter
noise_param = opti.parameter(); % noise variable
sigma_w = 10^(noise_param-3);   % noise

% Evaluate expression graph of OCP with opti variables
[obj,eq_constr] = f_OCP(x_mean, F, Pvar, K, D, sigma_w);

% constraints
opti.subject_to(eq_constr == 0);

% cost function
opti.minimize(obj);

% solver options
optionssol.ipopt.linear_solver = 'mumps';
optionssol.ipopt.tol = 1e-5; 
optionssol.ipopt.dual_inf_tol = 1e-5;
optionssol.ipopt.constr_viol_tol = 1e-7;

opti.solver('ipopt',optionssol);

%% Solve optimisation problem
% Plot variables
if make_fig
    figure();
    cs = linspecer(noise_factors(end)+1,'blue');
    cs = cs(2:end,:);
end

% OCP for seven different noise levels
for ii = 1:length(noise_factors)
    noise_factor = noise_factors(ii);

    try % error handling, so the sequence does not stop if solving one OCP fails

        % set parameter value
        opti.set_value(noise_param,noise_factor);

        % solve OCP
        sol = opti.solve();

        % extract results
        x_mean_sol = sol.value(x_mean);
        F_sol = sol.value(F);
        P_sol = sol.value(Pvar);
        K_sol = sol.value(K);
        D_sol = sol.value(D);
        F_sol = F_sol + K_sol*(x_mean_sol(1,:)-pi) + D_sol*x_mean_sol(2,:);
        
        % prepare results to return
        Results.(['R' num2str(noise_factor)]).stats = sol.stats;
        Results.(['R' num2str(noise_factor)]).x_mean = x_mean_sol;
        Results.(['R' num2str(noise_factor)]).F = F_sol;
        Results.(['R' num2str(noise_factor)]).P = P_sol;
        Results.(['R' num2str(noise_factor)]).K = K_sol;
        Results.(['R' num2str(noise_factor)]).D = D_sol;

        % Plot variables
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
    catch err_msg
        % if there was an error, pass it on with the results
        Results.(['R' num2str(noise_factor)]).failed = err_msg;
    end


end
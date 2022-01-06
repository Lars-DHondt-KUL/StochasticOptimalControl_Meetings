function [Results] = example2_2(make_fig,noise_factors,dt,fb)

% clear all; clc; close all;
% make_fig = 1;
% noise_factors = [1,3];


import casadi.*
% dt = 0.167; 
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

% Plot variables
if make_fig
    figure();
    cs = linspecer(noise_factors(end)+1,'blue');
    cs = cs(2:end,:);
end

% OCP for seven different noise levels
for ii = 1:length(noise_factors)
    noise_factor = noise_factors(ii);
    
    sigma_w = 10^(noise_factor-3);

    opti = casadi.Opti(); % Create opti instance
    
    F = opti.variable(1,N+1);      % force/torque
    x_mean = opti.variable(2,N+1); % mean trajectory
    Pvar = opti.variable(4,N+1);   % covariance matrix
    K = opti.variable(1,1);        % stiffness
    D = opti.variable(1,1);        % damping
    if ~fb
        opti.subject_to(K == 0);
        opti.subject_to(D == 0);
    end

    
    P_k = 0.0001*eye(2);        
    opti.subject_to(Pvar(:,1) == P_k(:));

    % expected effort
    obj = (F(:,1) + K*(x_mean(1,1)-pi) + D*x_mean(2,1))^2 + [K D]*P_k*[K;D];
    for i = 1:N
        % integration of mean state trajectory
        x_mean_dot = fcn_dyn(x_mean(:,i),F(i),0,K,D);
        x_mean_dot_next = fcn_dyn(x_mean(:,i+1),F(i+1),0,K,D);
        opti.subject_to(G_Trapezoidal(x_mean(:,i),x_mean(:,i+1),x_mean_dot,x_mean_dot_next,dt) == 0);
        
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
        opti.subject_to(Pvar(:,i+1) == P_k(:));
        
        % expected effort
        obj = obj + (F(:,i+1) + K*(x_mean(1,i+1)-pi) + D*x_mean(2,i+1))^2 + [K D]*P_k*[K;D]; % expected effort
    end
    
    % boundary conditions mean state
    opti.subject_to(x_mean(:,end) == [pi; 0]);
    opti.subject_to(x_mean(:,1) == 0);
    
    % cost function
    opti.minimize( norm(obj + 1e4*(P_k(1,1))) );  % Add minimization of variance, otherwise feedback gains will be zero and feedforward only is optimal solution
    
    optionssol.ipopt.linear_solver = 'mumps';
    optionssol.ipopt.tol = 1e-5; 
    optionssol.ipopt.dual_inf_tol = 1e-5;
    optionssol.ipopt.constr_viol_tol = 1e-7;
    
    opti.solver('ipopt',optionssol);
    
    try
        sol = opti.solve();
        x_mean_sol = sol.value(x_mean);
        F_sol = sol.value(F);
        P_sol = sol.value(Pvar);
        K_sol = sol.value(K);
        D_sol = sol.value(D);
        F_sol = F_sol + K_sol*(x_mean_sol(1,:)-pi) + D_sol*x_mean_sol(2,:);
        
        Results.(['R' num2str(noise_factor)]).stats = sol.stats;
        Results.(['R' num2str(noise_factor)]).x_mean = x_mean_sol;
        Results.(['R' num2str(noise_factor)]).F = F_sol;
        Results.(['R' num2str(noise_factor)]).P = P_sol;
        Results.(['R' num2str(noise_factor)]).K = K_sol;
        Results.(['R' num2str(noise_factor)]).D = D_sol;

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
        Results.(['R' num2str(noise_factor)]).failed = err_msg;
    end

end


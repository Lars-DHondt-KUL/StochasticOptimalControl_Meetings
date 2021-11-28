clear all; clc; close all;

import casadi.*
dt = 0.001; T = 10; t = 0:dt:T;
N = round(T/dt);



%- Stochastic Dynamics of the controlled system
% Variables initialisaiton
x_mean = MX.sym('x_mean',2);
w =  MX.sym('w',1); % disturbing force
F =  MX.sym('F',1);
K = MX.sym('X',1,1);
B = MX.sym('X',1,1);
% Dynamic System parameters
g = 9.81; l = 1.2; m = 5;

% Dynamics of stochastic problem ( xdot = f(x,u,...) )
x_mean_dot = [ x_mean(2); (-m*g*l*sin(x_mean(1)) + F + w)/(m*l^2)];
fcn_dyn = Function('fcn_dyn',{x_mean,F,w,K,B},{x_mean_dot});

% Linearization depending on the state.
fcn_A = Function('fcn_A',{x_mean,F,w,K,B},{jacobian(x_mean_dot,x_mean)});
fcn_C = Function('fcn_C',{x_mean,F,w,K,B},{jacobian(x_mean_dot,w)});

P_0 = 1e-6*eye(2); % small initial uncertainty



% Simulate P forward for upright position and visualize


% Simulate P forward for downward position and visualize




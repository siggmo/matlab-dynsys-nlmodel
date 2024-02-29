% Example demonstrating the use of the nlmodel framework to simulate a
% state feedback controller applied to the linearization and the nonlinear
% model.

clear variables;
opt=odeset('RelTol',1e-6);

% NONLINEAR MODEL
[rocket, par] = rocket2d_fullstate();

% LINEARIZATION
% obtain linear model at upright equilibrium
x_e = [0, 0, 0, 0, 0, 0]';
u_e = [par.m*par.g; 0; 0];  % upward thrust to keep the rocket hovering
ss_rocket = rocket.linearize(x_e, u_e);
[A, B, ~, ~] = ssdata(ss_rocket);

% CONTROLLER
p = -1:-1:-6;
F = place(A, B, p);
K = ss(-F);

% PLANT CONFIGURATION
% Plant with controller, but system missing. States and control input as
% performance output. w as generalized disturbance.
systemnames = 'K';
inputvar = '[yt{6}; w{3}]';
outputvar = '[K+w; yt(1:3); K+w]';
input_to_K = '[yt]';
Pk = sysic;

% INPUT FUNCTION
w = @(t) t*zeros([3 1]);    % zero disturbance input

% CLOSED LOOP
% linear
cl_lin = lft(ss_rocket, Pk);
[A, B, C, D] = ssdata(cl_lin);
f_lin = @(t, x) A*x + B*w(t); % closed loop linear system with w acting on it
h_lin = @(x, u) C*x + D*u; % linear system output map

% nonlinear
cl_nl = nl_upper_lft(rocket, Pk, x_e, u_e);
f_nl = @(t, x) cl_nl.f(x, w(t));
h_nl = @cl_nl.h;

% SIMULATION
% initialization
Ts = 0;
Te = 5;
T = Ts:0.05:Te;

x_0 = 0.1*[1 1 1 0 0 0]'; % deviation from equilibrium

sols = [];

% linear
lin_cl_sol = ode45(f_lin, [Ts Te], x_0, opt);
sols = add_sol(sols, lin_cl_sol, T, w, h_lin, 'linear', 'blue');

% nonlinear
nl_cl_sol = ode45(f_nl, [Ts Te], x_e + x_0, opt);
sols = add_sol(sols, nl_cl_sol, T, w, h_nl, 'nonlinear', 'green');

% PLOTTING
title = '2D rocket state feedback';
plot_sols(title, cl_lin.outputName, sols);
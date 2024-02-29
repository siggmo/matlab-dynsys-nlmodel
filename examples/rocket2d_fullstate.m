function [rocket, par] = rocket2d_fullstate()
%ROCKET2D Returns a nlmodel object of a 2D rocket with full state as
%output.
%
%   Parameter values are returned, too.
%
%   Model taken from: 'A Robust Control Approach for Rocket Landing' by
%   Reuben Ferrante, 2017, University of Edinburgh

    states = {'x'; 'z'; 'theta'; 'x dot'; 'z dot'; 'theta dot'};
    inputs = {'F_E'; 'F_S'; '\phi'};
    outputs = states;
    
    m = 1;
    J = 1;
    l_1 = 1;
    l_2 = 0.8;
    l_n = 0.1;
    g = 9.81;
    
    % model dynamics f
    f_4 = @(x, u) 1/m*(u(1)*sin(u(3)+x(3)) + u(2)*cos(x(3)));
    f_5 = @(x, u) 1/m*(u(1)*cos(u(3)+x(3)) - u(1)*sin(x(1)) - m*g);
    f_6 = @(x, u) 1/J*(-u(1)*sin(u(3))*(l_1+l_n*cos(u(3))) + u(2)*l_2);
    f = @(x, u) [x(4); x(5); x(6); f_4(x, u); f_5(x, u); f_6(x, u)];
    
    % output map
    h = @(x) [x(1); x(2); x(3); x(4); x(5); x(6)];

    % jacobians, as nested function
    function [A, B, C, D] = jacobians(x, u)
        df4dx3 = @(x, u) 1/m*(u(1)*cos(u(3)+x(3)) - u(2)*sin(x(3)));
        df5dx3 = @(x, u) -1/m*(u(1)*sin(u(3)+x(3)) + u(2)*cos(x(3)));
        Dxf = @(x, u) [zeros(3) eye(3);
                       0 0 df4dx3(x, u) 0 0 0;
                       0 0 df5dx3(x, u) 0 0 0;
                       0 0 0 0 0 0];
    
        Duf = @(x, u) [zeros(3);
                       1/m*sin(u(3)+x(3)) 1/m*cos(x(3)) 1/m*u(1)*cos(u(3)+x(3));
                       1/m*cos(u(3)+x(3)) -1/m*sin(x(3)) -1/m*u(1)*sin(u(3)+x(3));
                       -1/J*sin(u(3))*(l_1+l_n*cos(u(3))) 1/J*l_2 1/J*u(1)*(-cos(u(3))*(l_1+l_n*cos(u(3))) + l_n*(sin(u(3))^2))];
        
        Dxh = @(x) eye(6);
        Duh = zeros([6 3]);
    
        A = Dxf(x, u);
        B = Duf(x, u);
        C = Dxh(x);
        D = Duh;
    end
    
    % model of rocket
    rocket = nlmodel(f, h, 6, 3, 6, 'stateName', states, 'inputName', inputs, 'outputName', outputs, 'jacobians', @jacobians);

    % expose parameters
    par.m   = m;
    par.J   = J;
    par.l_1 = l_1;
    par.l_2 = l_2;
    par.l_n = l_n;
    par.g   = g;
end
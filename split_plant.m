function [A, B1, B2, C1, C2, D11, D12, D21, D22] = split_plant(P, nz, ny, nw, nu)
% SPLIT_PLANT extracts the submatrices of a generalized plant P.
%   
%   P is expected to have the following structure:
%   
%   /z\   /P11 P12\/w\
%   \y/ = \P21 P22/\u/
%   
%   where
%   z: performance output
%   y: measurement output
%   w: generalized disturbance
%   u: control input
%   
%   You have to specify the number of components of each channel (nz, ny,
%   nw, nu).

    arguments
        P ss
        nz (1, 1) double {mustBeInteger, mustBeGreaterThanOrEqual(nz, 0)}
        ny (1, 1) double {mustBeInteger, mustBeGreaterThanOrEqual(ny, 0)}
        nw (1, 1) double {mustBeInteger, mustBeGreaterThanOrEqual(nw, 0)}
        nu (1, 1) double {mustBeInteger, mustBeGreaterThanOrEqual(nu, 0)}
    end

    % extract the state space matrices of the plant P
    A = P.A;

    nx = size(A, 1);

    B1  = P.B(   1:nx   ,    1:nw   );  % w->x_dot
    B2  = P.B(   1:nx   , nw+1:nw+nu);  % u->x_dot

    C1  = P.C(   1:nz   ,    1:nx   );  % x->z
    C2  = P.C(nz+1:nz+ny,    1:nx   );  % x->y

    D11 = P.D(   1:nz   ,    1:nw   );  % w->z
    D12 = P.D(   1:nz   , nw+1:nw+nu);  % u->z
    D21 = P.D(nz+1:nz+ny,    1:nw   );  % w->y
    D22 = P.D(nz+1:nz+ny, nw+1:nw+nu);  % u->y
end
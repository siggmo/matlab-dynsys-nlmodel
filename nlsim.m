function [Y, X] = nlsim(M, u, T, x0)
%NLSIM Simulate time response of nonlinear dynamic systems to arbitrary inputs.
%
%   Uses ode45 solver to simulate the response of the dynamic system M to
%   the input function u=u(t). The response is returned at the time samples
%   specified by T. The initial condition can be specified with by x0, and
%   is set to zero if omitted.
%   Y and X have as many rows as there are outputs and states respectively,
%   and as many columns as there are time samples in T.

    arguments
        M nlmodel
        u (1, 1) function_handle
        T (1, :) double
        x0 (:, 1) double = 0
    end

    if x0 == 0
        x0 = zeros([M.nx 1]);
    end

    f = @(t, x) M.f(x, u(t));
    h = @M.h;
    Ts = T(1);
    Te = T(end);

    % static systems don't need to be integrated
    if M.nx ~= 0
        %opt = odeset('RelTol',1e-6);
        sol = ode45(f, [Ts Te], x0);%, opt);
        X = deval(sol, T);
    else
        X = x0 * ones([1 length(T)]);
    end

    Y = zeros([M.ny, length(T)]);
    
    if M.possibleDirectFeedThrough
        for i=1:length(T)
            Y(:, i) = h(X(:, i), u(T(i)));
        end
    else
        for i=1:length(T)
            Y(:, i) = h(X(:, i));
        end
    end
end
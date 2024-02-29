function M = ss2nlmodel(SYS, x_e, u_e, y_e)
%SS2NLMODEL convert state space model to nlmodel.
%
%   Constructs a nlmodel object from a linear state space model,
%   and applies the offsets (if given, otherwise they are assumed to be
%   zero).
%   Be aware that delays, uncertainties etc. are not captured here, only 
%   the dynamics given by the four system matrices A, B, C, D.
%
%   USAGE
%   * M = ss2nlmodel(SYS)
%       where SYS is a ss object (zero offsets assumed). Creates a nlmodel
%       with the following dynamics
%           f = f(x, u) = A*x + B*u
%           h = h(x, u) = C*x + D*u (or if D=0: h = h(x) = C*x)
%   * M = ss2nlmodel(SYS, x_e, u_e, y_e)
%       where SYS is a ss object, and x_e, u_e and y_e are offsets. Creates
%       a nlmodel with these dynamics:
%           f = f(x, u) = A*(x-x_e) + B*(u-u_e)
%           h = h(x, u) = y_e + C*(x-x_e) + D*(u-u_e) (h=h(x) if D=0)

    arguments
        SYS StateSpaceModel
        x_e (:, 1) double = 0
        u_e (:, 1) double = 0
        y_e (:, 1) double = 0
    end

    % extract system matrices
    % TODO: check if model is really supported. How?
    [A, B, C, D] = ssdata(SYS);

    % determine dimensions
    nx = size(A, 1);
    nu = size(B, 2);
    ny = size(C, 1);
    possibleDirectFeedthrough = any(D, "all") & ~isempty(D);

    % if no offset was given, use omit offsets in formulas. Otherwise check
    % dimensions and apply offsets.
    if nargin == 1
        % model dynamics
        f = @(x, u) A*x + B*u;
        
        % output map
        if possibleDirectFeedthrough
            h = @(x, u) C*x + D*u;
        else
            h = @(x) C*x;
        end
    elseif nargin == 4
        % check dimension of offsets
        if ~(length(x_e) == nx)
            error("MATLAB:nl_sys_toolbox:ss2nlmodel:state_offset_dim", "Incompatible dimension of state offset!");
        end
        if ~(length(u_e) == nu)
            error("MATLAB:nl_sys_toolbox:ss2nlmodel:input_offset_dim", "Incompatible dimension of input offset!");
        end
        if ~(length(y_e) == ny)
            error("MATLAB:nl_sys_toolbox:ss2nlmodel:output_offset_dim", "Incompatible dimension of output offset!");
        end

        % model dynamics
        f = @(x, u) A*(x - x_e) + B*(u - u_e);
        
        % output map
        if possibleDirectFeedthrough
            h = @(x, u) C*(x - x_e) + D*(u - u_e) + y_e;
        else
            h = @(x) C*(x - x_e) + y_e;
        end
    else
        error("MATLAB:nl_sys_toolbox:ss2nlmodel:argument_count", "ss2nlmodel requires exactly one or four input arguments!");
    end

    % jacobians: linear model is its own linearization
    function [A_, B_, C_, D_] = jacobians(x, u)
        A_ = A;
        B_ = B;
        C_ = C;
        D_ = D;
    end

    % channel names
    stateName = SYS.stateName;
    inputName = SYS.inputName;
    outputName = SYS.outputName;

    % nlmodel
    M = nlmodel(f, h, nx, nu, ny, 'stateName', stateName, 'inputName', inputName, 'outputName', outputName, 'jacobians', @jacobians);
end
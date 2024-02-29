function M = nl_upper_lft(M_nl, P, x_e, u_e)
    %NL_UPPER_LFT perform upper LFT on a nlmodel and ss object.
    % 
    %   USAGE
    %   M = NL_UPPER_LFT(M_nl, P)             (without offsets / zero offsets)
    %   or
    %   M = NL_UPPER_LFT(M_nl, P, x_e, u_e)   (with offsets)
    %
    %   OVERVIEW
    %   M = NL_UPPER_LFT(M_nl, P, x_e, u_e) forms the following feedback  
    %   interconnection of the nonlinear system M_nl and the linear system P:
    %		
    %                        +-------+
    %                        |       |
    %                  +---->| M_nl  |-----+
    %                  |     |       |     |
    %                  |     +-------+     |
    %          u + u_e |                   | y - y_e, (y_e = h(x_e, u_e))
    %                  |     +-------+     |
    %                  +-----|       |<----+
    %                        |   P   |
    %            z <---------|       |-------- w
    %                        +-------+
    %
    %   The feedback loop connects the upper outputs of P (as many as M_nl has 
    %   inputs) to the inputs of M_nl and all the outputs of M_nl to the upper 
    %   inputs of P. P therefore needs to have at least as many inputs as M_nl 
    %   has outputs and at least as many outputs as M_nl has inputs. The 
    %   remaining inputs and outputs of P are exposed as inputs and outputs of 
    %   M. The resulting system M is a nlmodel object and maps the input vector 
    %   w to the output vector z. This operation is referred to as a upper 
    %   linear fractional transformation or upper LFT.
    %   
    %   OFFSETS
    %   When dealing with linearizations of nonlinear systems, one typically
    %   linearizes around equilibrium points. If one is specified, nl_upper_lft  
    %   takes care of the offsets. Otherwise, they are assumed to be zero.
    %   To the input w and output z no offsets are being applied.
    %
    %   DIRECT FEED THROUGH
    %   Attention: Direct feed through in the nonlinear system is not allowed!
    %   Therefore only nlmodel-objects with output map h=h(x) must be used. The
    %   system P however may very well have a direct feed through term, which  
    %   for example is inevitable for state feedback controllers.
    %
    %   INTERPRETATION AS GENERALIZED PLANT
    %   P can be viewed as a generalized plant. The standard configuration of 
    %   a generalized plant is:
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
    %   This function was intended to interconnect a nonlinear system with 
    %   controllers that were designed for the linearized system. In this
    %   context, P would contain the controller K together with the 
    %   configuration for performance output z_p and generalized disturbance 
    %   w_p. Such a plant would look like this:
    %   
    %   /z_nl\   /P11 P12\/w_nl\      (z = z_nl, w = w_nl)
    %   \z_p / = \P21 P22/\w_p /      (y = z_p,  u = w_p)
    %
    %   The inputs and outputs are then typically defined as follows:
    %   z = z_nl: control input to the linearized system
    %   y = z_p : performance output
    %   w = w_nl: measurement output from the linearized system
    %   u = w_p : generalized disturbance
    %   
    %   Then, when interconnecting with M_nl using nl_upper_lft, z_nl goes to 
    %   the inputs of M_nl, and the output of M_nl goes to w_nl, with offsets
    %   being applied respectively to take the point of linearization into 
    %   account.

    arguments
        M_nl nlmodel
        P ss
        x_e (:, 1) double = 0
        u_e (:, 1) double = 0
    end
    
    % make sure model does not have direct feed through
    if M_nl.possibleDirectFeedThrough
        error("MATLAB:nl_sys_toolbox:nl_upper_lft:direct_feedthrough", "Output map of M_nl must not depend on u, but only on x! Direct feed through is not supported.");
    end

    % determine state dimensions
    n_nl = M_nl.nx;     % number of components of nonlinear state vector
    n_p = size(P.A, 1); % number of components of linear system state vector

    % determine and check input/output dimensions
    n_z_nl = M_nl.nu;                   % z_nl is the input to M_nl
    n_z_p = size(P.C, 1) - n_z_nl;      % z_p contains the remaining output components (performance output)
    if n_z_p < 0
        error("MATLAB:nl_sys_toolbox:nl_upper_lft:too_few_outputs", "Plant has too few output components for interconnection!");
    end

    n_w_nl = M_nl.ny;                   % the output of M_nl enters P via w_nl
    n_w_p = size(P.D, 2) - n_w_nl;      % w_p contains the remaining inputs components (generalized disturbance)
    if n_w_p < 0
        error("MATLAB:nl_sys_toolbox:nl_upper_lft:too_few_inputs", "Plant has too few input components for interconnection!");
    end
    
    % extract partitioned matrices
    [A, B1, B2, C1, C2, D11, D12, D21, D22] = split_plant(P, n_z_nl, n_z_p, n_w_nl, n_w_p);

    % Construct closed-loop dynamics. States stacked, M_nl top, plant bottom.
    % If no offset was given, omit them in the formulas. This reduces
    % numerical errors. Otherwise check dimensions.
    h = @M_nl.h;
    f = @M_nl.f;
    if nargin == 2
        % dynamics without any offsets
        f_cl_ = @(x_nl, x_p, w_p) [f(x_nl, C1*x_p + D11*h(x_nl) + D12*w_p);
                                   A*x_p + B1*h(x_nl) + B2*w_p];
        h_cl_ = @(x_nl, x_p, w_p) C2*x_p + D21*h(x_nl) + D22*w_p;

        % still set offsets to zero vectors of appropriate dimension, since
        % they are needed later for the jacobians function. It's not
        % allowed to define functions inside if-statements, so the
        % jacobians-function has to incorporate offsets, even if they are
        % zero.
        x_e = zeros([M_nl.nx 1]);
        u_e = zeros([M_nl.nu 1]);
        y_e = M_nl.h(x_e);         % compute output offset
    elseif nargin == 4
        % check offset dimensions
        if ~isequal(length(x_e), M_nl.nx)
            error("MATLAB:nl_sys_toolbox:nl_upper_lft:state_offset", "Dimension of state offset does not match!");
        end
        if ~isequal(length(u_e), M_nl.nu)
            error("MATLAB:nl_sys_toolbox:nl_upper_lft:input_offset", "Dimension of input offset does not match!");
        end

        % compute output offset
        y_e = M_nl.h(x_e);

        % dynamics with offsets applied
        f_cl_ = @(x_nl, x_p, w_p) [f(x_nl, C1*x_p + D11*(h(x_nl)-y_e) + D12*w_p + u_e);
                                   A*x_p + B1*(h(x_nl)-y_e) + B2*w_p];
        h_cl_ = @(x_nl, x_p, w_p) C2*x_p + D21*(h(x_nl)-y_e) + D22*w_p;
    else
        error("MATLAB:nl_sys_toolbox:nl_upper_lft:argument_count", "nl_upper_lft requires exactly two or four input arguments!");
    end

    % use stacked state only now, for better readability of above formulas
    f_cl = @(x, u) f_cl_(x(1:n_nl, :), x(n_nl+1:n_nl+n_p, :), u);
    h_cl = @(x, u) h_cl_(x(1:n_nl, :), x(n_nl+1:n_nl+n_p, :), u);

    % channel names of closed loop
    stateName_cl = [M_nl.stateName; P.stateName];
    inputName_cl = P.inputName(n_w_nl+1:n_w_nl+n_w_p, :);
    outputName_cl = P.outputName(n_z_nl+1:n_z_nl+n_z_p, :);

    % dimensions of closed loop
    nx_cl = n_nl+n_p;
    nu_cl = n_w_p;
    ny_cl = n_z_p;

    % jacobians of closed loop, using chain rule.
    function [A_cl, B_cl, C_cl, D_cl] = jacobians_cl_(x, w_p)
        x_nl = x(1:n_nl, :);
        x_p = x(n_nl+1:n_nl+n_p, :);
        u_ = C1*x_p + D11*(h(x_nl)-y_e) + D12*w_p + u_e;

        [A_nl, B_nl, C_nl, ~] = M_nl.jacobians(x_nl, u_); % no direct feedthrough -> D_nl=0

        A_cl = [A_nl + B_nl*D11*C_nl, B_nl*C1;
                B1*C_nl             , A];
        B_cl = [B_nl*D12;
                B2];
        C_cl = [D21*C_nl, C2];
        D_cl = D22;
    end

    % only pass jacobians function if the jacobians of the nonlinear system are
    % available.
    if M_nl.jacobiansPresent
        M = nlmodel(f_cl, h_cl, nx_cl, nu_cl, ny_cl, 'stateName', stateName_cl, 'inputName', inputName_cl, 'outputName', outputName_cl, 'jacobians', @jacobians_cl_);
    else
        M = nlmodel(f_cl, h_cl, nx_cl, nu_cl, ny_cl, 'stateName', stateName_cl, 'inputName', inputName_cl, 'outputName', outputName_cl);
    end
end
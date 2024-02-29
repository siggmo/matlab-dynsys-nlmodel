function [testResult, samplingInput] = compareToLinearSS(M_nl, M_lin)
    %COMPARETOLINEARSS  Compares a nlmodel object with an ss model object
    %
    %   Three different comparison methods are being applied:
    %    - compare linearization of nlmodel with ss model if possible. This is
    %       done by extracting the system matrices and taking the 2-norm of the
    %       difference.
    %    - compare nlmodel objects by converting ss model to an nlmodel. Then
    %       comparison is done by randomly sampling the difference of both
    %       vector fields f and h of both models and taking the maximum.
    %    - simulate both models using nlsim and compare the system responses to
    %       some arbitrarily chosen input signal and zero initial value.
    
        arguments
            M_nl nlmodel
            M_lin ss
        end
    
        %%% CHECK DIMENSIONS
        if not (size(M_nl) == size(M_lin) & M_nl.nx == length(M_lin.stateName))
            error("MATLAB:nl_sys_toolbox:nlmodel:test:compareModel", "Models must have the same number of states, inputs and outputs to be able to be compared.");
        end

        %%% CHANNEL NAMES
        % check if channel names are the same
        if isequal(M_nl.stateName, M_lin.stateName) & isequal(M_nl.inputName, M_lin.inputName) & isequal(M_nl.outputName, M_lin.outputName)
            channelNamesOk = true;
        else
            channelNamesOk = false;
        end
    
        %%% MODEL COMPARISON: LINEAR SS MODEL OBJECTS
        % linearize nlmodel object for comparison, if possible
        % this only validates the implementation of the jacobians in the nlmodel!
        if M_nl.jacobiansPresent
            M_nl_ss = M_nl.linearize();
            [A_lin, B_lin, C_lin, D_lin] = ssdata(M_lin);
            [A_nl, B_nl, C_nl, D_nl] = ssdata(M_nl_ss);
            sysMatrixError = norm([A_lin - A_nl, B_lin - B_nl; C_lin - C_nl, D_lin - D_nl], 2);
        end
    
        %%% MODEL COMPARISON: NLMODEL OBJECTS
        % here the vector fields are compared by evaluating them at random points.
        n = 1000;
    
        % initialize variables
        diff_f = zeros([n 1]);
        diff_h = zeros([n 1]);
        X = zeros([M_nl.nx n]);
        U = zeros([M_nl.nu n]);
    
        % evaluate at n random points
        for i = 1:n
            x = (rand(M_nl.nx, 1) - 0.5) * 1000;
            u = (rand(M_nl.nu, 1) - 0.5) * 1000;
    
            diff_f(i) = norm(M_nl.f(x, u) - (M_lin.A*x + M_lin.B*u), 2);
            diff_h(i) = norm(M_nl.h_xu(x, u) - (M_lin.C*x + M_lin.D*u), 2);
            X(:, i) = x;
            U(:, i) = u;
        end
    
        % store results and also those points where the difference is maximal
        [samplingErrorF, idx] = max(diff_f);    % no abs() needed since all values are positive
        samplingInput.maxErrF_x = X(:, idx);
        samplingInput.maxErrF_u = U(:, idx);
        [samplingErrorH, idx] = max(diff_h);
        samplingInput.maxErrH_x = X(:, idx);
        samplingInput.maxErrH_u = U(:, idx);
    
        %%% MODEL COMPARISON: SIMULATION
        % Input function of appropriate dimension and time vector
        u = @(t) [sin(t)]*ones([M_nl.nu 1]);
        T = 0:0.01:5;
    
        % simulate ss model. Use own function nlsim, since lsim didn't work properly.
        % TODO: check why lsim doesn't work!!!
        M_lin_nl = nlmodel(M_lin);
        [Y_lin, X_lin] = nlsim(M_lin_nl, u, T);
    
        % simulate nlmodel itself
        [Y_nl, X_nl] = nlsim(M_nl, u, T);
    
        % compare
        % TODO: this might not be a good choice of norms. It would be much better to
        % compute the energy of the difference, i.e. taking the integral over
        % the 2-norm of the difference of the trajectory points.
        simulationErrorX = norm(X_lin - X_nl, "inf");
        simulationErrorY = norm(Y_lin - Y_nl, "inf");
    
        % combine results into struct
        testResult.channelNamesOk = channelNamesOk;
        testResult.sysMatrixError = sysMatrixError;
        testResult.samplingErrorF = samplingErrorF;
        testResult.samplingErrorH = samplingErrorH;
        testResult.simulationErrorX = simulationErrorX;
        testResult.simulationErrorY = simulationErrorY;
    end
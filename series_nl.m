function M = series_nl(M1, M2, n)
    %SERIES_NL Series connection of two nlmodel objects.
    %
    %                                  +------+
    %                           v2 --->|      |
    %                  +------+        |  M2  |-----> y2
    %                  |      |------->|      |
    %         u1 ----->|      |y1   u2 +------+
    %                  |  M1  |
    %                  |      |---> z1
    %                  +------+
    %
    %   M = SERIES_NL(M1,M2,n) connects the nlmodel objects M1 and M2 in 
    %   series. The n upper outputs of M1 are connected with the n lower
    %   inputs of M2. If n is omitted or equals -1, SERIES_NL connects M1 and 
    %   M2 in cascade and returns M = M2 * M1.
    %   Setting n=0 naturally results in M = [M2; M1].
    
        arguments
            M1 nlmodel
            M2 nlmodel
            n (1, 1) double {mustBeInteger, mustBeGreaterThanOrEqual(n, -1)} = -1
        end
    
        % if n was omitted, set it to the maximum possible value
        if n == -1
            n = min(M1.ny, M2.nu);
        end
    
        % check if input/output dimensions allow for a series connection
        if ~((n <= M1.ny) && (n <= M2.nu))
            error("MATLAB:nl_sys_toolbox:series:illegal_n", "The specified number of to be connected inputs/outputs exceeds the available number of inputs/outputs for the given systems!");
        end
    
        % gather data. Use h_xu instead of h, because then the most general
        % case can be implemented.
        nx1 = M1.nx;
        nu1 = M1.nu;
        ny1 = n;
        nz1 = M1.ny - ny1;
    
        nx2 = M2.nx;
        nu2 = n;
        nv2 = M2.nu - nu2;
        ny2 = M2.ny;
    
        f1 = @M1.f;
        h1 = @M1.h_xu;
        f2 = @M2.f;
        h2 = @M2.h_xu;
    
        % split output function of M1
        function y1 = h11(x1, u1)
            h_eval = h1(x1, u1);
            y1 = h_eval(1:ny1);
        end
        function z1 = h12(x1, u1)
            h_eval = h1(x1, u1);
            z1 = h_eval(ny1+1:ny1+nz1);
        end
    
        % construct dynamics of series connection
        f_ = @(x1, x2, v2, u1) [f1(x1, u1);
                                f2(x2, [v2; h11(x1, u1)])];
        f = @(x, u) f_(x(1:nx1, :), x(nx1+1:nx1+nx2, :), u(1:nv2, :), u(nv2+1:nv2+nu1, :));
    
        % construct output map of series connection
        h_ = @(x1, x2, v2, u1) [h2(x2, [v2; h11(x1, u1)]);
                                h12(x1, u1)];
        % At this point we can only guarantee that no direct feed through is
        % possible for M if both for M1 and M2 a direct feed through is 
        % impossible, too.
        if ~M1.possibleDirectFeedThrough && ~M2.possibleDirectFeedThrough
            h = @(x) h_(x(1:nx1, :), x(nx1+1:nx1+nx2, :), zeros([nv2 1]), zeros([nu1 1]));
        else
            h = @(x, u) h_(x(1:nx1, :), x(nx1+1:nx1+nx2, :), u(1:nv2, :), u(nv2+1:nv2+nu1, :));
        end    
    
        % dimensions of series connection
        nx = nx1 + nx2;
        nu = nv2 + nu1;
        ny = ny2 + nz1;
    
        % channel names of series connection
        stateName = [M1.stateName; M2.stateName];
        inputName = [M2.inputName(1:nv2); M1.inputName];
        outputName = [M2.outputName; M1.outputName(ny1+1:ny1+nz1)];
    
        M = nlmodel(f, h, nx, nu, ny, "stateName", stateName, "inputName", inputName, "outputName", outputName);
    end
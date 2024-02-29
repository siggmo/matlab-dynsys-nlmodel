classdef nlmodelTest < matlab.unittest.TestCase
    properties (Constant)
        Tol = eps*10^6;
        int_state = [0 10];
        int_io = [1 10];
    end

    methods(Test)
        %%% Test if operations on linear models yield the same result as
        %%% when using the well tested matlab functions for dynamic system
        %%% models.
        function AdditionLinear(testCase)
            % generate test input data
            nx1 = randi(testCase.int_state);
            nx2 = randi(testCase.int_state);
            nu = randi(testCase.int_io);
            ny = randi(testCase.int_io);

            M1 = rss(nx1, ny, nu);
            M2 = rss(nx2, ny, nu);
            %disp(M1.A);
            %disp(M2.A);

            M1_nl = nlmodel(M1);
            M2_nl = nlmodel(M2);

            % perform operation both on nlmodel and ss objects
            M_nl = M1_nl + M2_nl; % add nlmodel objects
            M_lin = M1 + M2;      % add ss objects

            compare(testCase, M_nl, M_lin);
        end

        function CompositionLinear(testCase)
            M1 = rss(6, 2, 4);
            M2 = rss(3, 4, 2);

            M1_nl = nlmodel(M1);
            M2_nl = nlmodel(M2);

            M_nl = M1_nl * M2_nl; % multiply nlmodel objects
            M_lin = M1 * M2;      % multiply ss objects

            compare(testCase, M_nl, M_lin);
        end

        function VertcatLinear(testCase)
            M1 = rss(6, 4, 2);
            M2 = rss(3, 4, 2);

            M1_nl = nlmodel(M1);
            M2_nl = nlmodel(M2);

            M_nl = [M1_nl; M2_nl]; % vertcat nlmodel objects
            M_lin = [M1; M2];      % vertcat ss objects

            compare(testCase, M_nl, M_lin);
        end

        function HorzcatLinear(testCase)
            M1 = rss(6, 4, 2);
            M2 = rss(3, 4, 2);

            M1_nl = nlmodel(M1);
            M2_nl = nlmodel(M2);

            M_nl = [M1_nl M2_nl]; % horzcat nlmodel objects
            M_lin = [M1 M2];      % horzcat ss objects

            compare(testCase, M_nl, M_lin);
        end

        function additionLinear_oneEmptyStateSystem(testCase)
            % this particular configuration threw an error earlier.
            M1 = rss(1, 1, 4);
            M2 = rss(0, 1, 4);

            % create nlmodels from random linear state space models
            M1_nl = nlmodel(M1);
            M2_nl = nlmodel(M2);
            
            % perform operation on both ss and nlmodel objects
            M_nl = M1_nl + M2_nl; % add nlmodel objects
            M_lin = M1 + M2;      % add ss objects
    
            compare(testCase, M_nl, M_lin);
        end

        function additionLinear_twoEmptyStateSystems(testCase)
            % this particular configuration threw an error earlier:
            %   Error using odearguments
            %   When the first argument to ode45 is a function handle, the tspan and y0 arguments must be supplied.
            %   
            %   Error in ode45 (line 107)
            %     odearguments(odeIsFuncHandle,odeTreatAsMFile, solver_name, ode, tspan, y0, options, varargin);
            %   
            %   Error in nlsim (line 28)
            %       sol = ode45(f, [Ts Te], x0);%, opt);
            %   
            %   Error in compareToLinearSS (line 47)
            %       Y_lin = nlsim(M_lin_nl, u, T);
            %   
            %   Error in nlmodelTestAdditionLinear (line 50)
            %       [diffL2, diffHinf, samplingError, simulationError] = compareToLinearSS(M_nl, M_lin);

            M1 = rss(0, 7, 3);
            M2 = rss(0, 7, 3);

            % create nlmodels from random linear state space models
            M1_nl = nlmodel(M1);
            M2_nl = nlmodel(M2);
            
            % perform operation on both ss and nlmodel objects
            M_nl = M1_nl + M2_nl; % add nlmodel objects
            M_lin = M1 + M2;      % add ss objects
    
            compare(testCase, M_nl, M_lin);
        end

        function testPossibleDirectFeedThroughDetection(testCase)
            % for the following system the direct feed through was not
            % detected, although the D matrix obviously is unequal to zero.
            A = ...
              [-1.4253334647493592 1.2514503443878728 1.8805882300735153;
               -2.0684085917398694 -3.3102941372746963 2.0557921646435946;
               0.90797904886415148 -2.6340971056748814 -2.4080266827425891];
            
            B = ...
              [-1.5470965138061403 0 0.1518074558867922 -0;
               0.021011792198648668 0 0 -1.5632119583496094;
               0 -1.5176353831941849 0.1514644302194455 -0.40237248489200927];
            
            C = ...
              [0.17488024622359219 0.57607358600852177 0.65974079998493518;
               0.91180619197093427 -0.37203037983946863 -0.457321769085304];
            
            D = ...
              [-0 0 0 -1.2907106948694238;
               -0 0.22161743457856298 -0 0.20065188402208484];

            M = nlmodel(ss(A, B, C, D));

            testCase.assertTrue(M.possibleDirectFeedThrough);
        end

        function multiplicationLinear_matMulError(testCase)
            % this particular configuration threw an error earlier:
            % Percent done: 23.0Error using  * 
            % Incorrect dimensions for matrix multiplication. Check that the number of
            % columns in the first matrix matches the number of rows in the second
            % matrix. To operate on each element of the matrix individually, use TIMES
            % (.*) for elementwise multiplication.
            % 
            % Error in ss2nlmodel>@(x,u)A*x+B*u (line 43)
            %         f = @(x, u) A*x + B*u;
            % 
            % Error in nlmodel>@(varargin)M1.f(varargin{:}) (line 336)
            %             f1 = @M1.f;
            % 
            % Error in nlmodel>@(x,u)[f1(x(1:nx1),h2(x(nx1+1:nx1+nx2),u));f2(x(nx1+1:nx1+nx2),u)] (line 345)
            %                 f_prod = @(x, u) [f1(x(1:nx1), h2(x(nx1+1:nx1+nx2), u));
            % 
            % Error in nlmodel/checkDynamics (line 620)
            %                 f_eval = M.f(testState, testInput);
            % 
            % Error in nlmodel/checkAll (line 688)
            %             M.checkDynamics();
            % 
            % Error in nlmodel/nlmodel_ (line 143)
            %                 M.checkAll();
            % 
            % Error in nlmodel (line 164)
            %                 M = M.nlmodel_(varargin{:});
            % 
            % Error in  *  (line 378)
            %                 M = nlmodel(f_prod, h_prod, nx_prod, nu_prod, ny_prod, 'stateName', stateName_prod, 'inputName', inputName_prod, 'outputName', outputName_prod, 'jacobians', jacobians_prod);
            % 
            % Error in nlmodelTestAdditionLinear (line 111)
            %             M_nl = M1_nl * M2_nl;
            % 
            % Caused by:
            %     Provided function handle f does not meet the requirements. It
            %     doesn't accept column vectors of dimensions nx and nu, or contains
            %     some other errors.
            %     Failed to create nlmodel object.
            % 
            % Related documentation        

            %%% REASON FOR ERROR: slicing was not correctly implemented in
            %%% the mtimes function. changed all occurences of x(i:j) to
            %%% x(i:j, :). This fixed the problem.
            
            nx1 = 0;
            nx2 = 1;
            nu1 = 4;
            nu2 = 4;
            ny1 = 7;
            ny2 = 4;

            M1 = rss(nx1, ny1, nu1);
            M2 = rss(nx2, ny2, nu2);

            % create nlmodels from random linear state space models
            M1_nl = nlmodel(M1);
            M2_nl = nlmodel(M2);
            
            % perform operation on both ss and nlmodel objects
            M_nl = M1_nl * M2_nl;
            M_lin = M1 * M2;
    
            compare(testCase, M_nl, M_lin);
        end
    end

    methods
        function compare(testCase, M_nl, M_lin)
            % compare both models
            [testResult, ~] = compareToLinearSS(M_nl, M_lin);
            channelNamesOk = testResult.channelNamesOk;
            sysMatrixError = testResult.sysMatrixError;
            samplingErrorF = testResult.samplingErrorF;
            samplingErrorH = testResult.samplingErrorH;
            simulationErrorX = testResult.simulationErrorX;
            simulationErrorY = testResult.simulationErrorY;
            
            % check comparison results
            testCase.assertTrue(channelNamesOk);
            testCase.assertLessThan(sysMatrixError, testCase.Tol);
            testCase.assertLessThan(samplingErrorF, testCase.Tol);
            testCase.assertLessThan(samplingErrorH, testCase.Tol);
            testCase.assertLessThan(simulationErrorX, testCase.Tol);
            testCase.assertLessThan(simulationErrorY, testCase.Tol);
        end
    end
end
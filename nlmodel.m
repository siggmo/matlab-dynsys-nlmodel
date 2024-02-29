classdef nlmodel
% NLMODEL objects representing a nonlinear state space I/O model.
%
%   GENERAL DESCRIPTION
%   Class representing nonlinear continuous time state space models 
%   consisting of a right side f (vector field) and an output map h:
%       dx/dt = f(x(t), u(t))
%           y = h(x(t), u(t)) or y = h(x(t))
%   Description cell arrays as well as the jacobians can additionally be 
%   added. Property validation and consistency checks (mostly dimensions) 
%   are performed automatically.
%
%   USAGE
%   Objects can be instatiated in four ways:
%   * M = nlmodel(f, h, nx, nu, ny)
%       where f=f(x, u) and h=h(x, u) or h=h(x) are function handles of
%       dynamics and output map, and nx, nu and ny specify the dimension of
%       the model.
%       Only in this syntax optional name-value arguments are supported to
%       set the properties inputName, stateName, outputName and jacobians. 
%       See the respective property descriptions for details.
%       For example:
%           M = nlmodel(@(x, u) x + u, @(x) x, 1, 1, 1, "inputName",
%           "controlinput", "stateName", "x", "outputName", "y")
%   * M = nlmodel(SYS)
%       where SYS is a ss object (zero offsets assumed). Creates a nlmodel
%       with the following dynamics
%           f = f(x, u) = A*x + B*u
%           h = h(x, u) = C*x + D*u (or if D=0: h = h(x) = C*x
%   * M = nlmodel(SYS, x_e, u_e, y_e)
%       where SYS is a ss object, and x_e, u_e and y_e are offsets. Creates
%       a nlmodel with these dynamics:
%           f = f(x, u) = A*(x-x_e) + B*(u-u_e)
%           h = h(x, u) = y_e + C*(x-x_e) + D*(u-u_e) (h=h(x) if D=0)
%   * M = nlmodel(M) where M is already a nlmodel object. Then the
%       constructor just returns M without altering it.
%   For detailed information on the input arguments look at the respective
%   property descriptions.
%
%   DIRECT FEED THROUGH
%   Some functions require that no direct feed through is present in the
%   model. Since that is difficult to check for arbitrary functions
%   h=h(x, u), one needs to explicitly make the output map independent of
%   u, i.e. h=h(x), if one wants to use the above mentioned functions.
%   The dependent property possibleDirectFeedThrough therefore checks the
%   number of arguments of h.

    %% Class Properties
    % MATHEMATICAL MODEL and DIMENSIONS
    properties (SetAccess = protected)
        % those default values don't mean anything since they are being
        % overwritten in the constructor. The function_handle default value,
        % set to @deal, is actually necessary since matlab cannot construct a
        % default value of class function_handle.
        f   (1, 1) function_handle        = @deal           % model dynamics, f=f(x, u)
        h   (1, 1) function_handle        = @deal           % output map, h=h(x, u) or h=h(x)

        nx  (1, 1) double {mustBeInteger} = 0               % number of state components
        nu  (1, 1) double {mustBeInteger} = 0               % number of input components
        ny  (1, 1) double {mustBeInteger} = 0               % number of output components

        jacobiansPresent    (1, 1) logical         = false  % tells if jacobians are available
    end
    
    % DIRECT FEED THROUGH
    properties(Dependent)
        possibleDirectFeedThrough       % true if h=h(x, u), i.e. h possibly dependent on input u.
        h_xu                            % returns a function handle h=h(x, u), even if actually h=h(x)
    end

    % JACOBIANS
    properties (Dependent)
        jacobians                       % returns jacobians of f and h at point (x, u): [A, B, C, D]. Throws an error if jacobians function handle is not present.
    end
    properties (Access=private)
        jacobians_          (1, 1) function_handle = @deal  % jacobians of f and h: [A, B, C, D] = jac(x, u)
    end

    % CHANNEL NAMES
    % this doubled structure for the name cell arrays is needed because,
    % for validating the size of these arrays in the set method, one has to
    % access the nx, nu, ny properties. This can cause problems with
    % matlabs loading procedure.
    properties (Dependent)
        stateName   (:, 1) cell {mustBeText}    % State names. Default empty
        inputName   (:, 1) cell {mustBeText}    % Input channel names. Default empty
        outputName  (:, 1) cell {mustBeText}    % Output channel names. Default empty
    end
    properties (Access=private)
        stateName_   (:, 1) cell {mustBeText} = {}
        inputName_   (:, 1) cell {mustBeText} = {}
        outputName_  (:, 1) cell {mustBeText} = {}
    end

    %% Class Methods
    methods (Hidden)
        % actual constructor
        function M = nlmodel_(M, f, h, nx, nu, ny, NameValueArgs)
        % create nlmodel object. For information on the arguments, see
        % property descriptions.

            arguments
                M
                f
                h
                nx
                nu
                ny
                NameValueArgs.stateName
                NameValueArgs.inputName
                NameValueArgs.outputName
                NameValueArgs.jacobians
            end

            %%% default description cell arrays
            if ~isfield(NameValueArgs, 'stateName')
                NameValueArgs.stateName = repmat({''}, nx, 1);  % cellstr("x" + string(1:nx))';
            end
            if ~isfield(NameValueArgs, 'inputName')             
                NameValueArgs.inputName = repmat({''}, nu, 1);  % cellstr("u" + string(1:nu))';
            end
            if ~isfield(NameValueArgs, 'outputName')
                NameValueArgs.outputName = repmat({''}, ny, 1);  % cellstr("y" + string(1:ny))';
            end

            %%% store everything (executes property validation)
            M.nx = nx;
            M.nu = nu;
            M.ny = ny;
            M.f = f;
            M.h = h;
            M.stateName = NameValueArgs.stateName;
            M.inputName = NameValueArgs.inputName;
            M.outputName = NameValueArgs.outputName;
            if isfield(NameValueArgs, 'jacobians')
                M.jacobiansPresent = true;
                M.jacobians_ = NameValueArgs.jacobians;
            else
                M.jacobiansPresent = false;
            end

            %%% check consistency
            try
                M.checkAll();
            catch ME
                msg = "Failed to create nlmodel object.";
                causeException = MException("MATLAB:nl_sys_toolbox:nlmodel:creation_failed", msg);
                ME = addCause(ME, causeException);
                rethrow(ME);
            end
        end
    end

    methods
        %% CONSTRUCTOR
        % wrapper constructor, in order to allow for models to be created
        % from ss objects. Offsets may be specified (x_e, u_e, y_e).
        % If a nlmodel object is given, then it is just returned unaltered.
        function M = nlmodel(varargin)
            if class(varargin{1}) == "ss" && (nargin == 1 || nargin == 4)
                M = ss2nlmodel(varargin{:});
            elseif class(varargin{1}) == "nlmodel" && nargin == 1
                M = varargin{1};
            else
                M = M.nlmodel_(varargin{:});
            end
        end


        %% LINEARIZATION
        function SYS = linearize(M, x_e, u_e)
        %LINEARIZE Obtains linearization of nlmodel object at equilibrium.
        %
        %   Returns a ss object with system matrices A, B, C, D obtained
        %   from the jacobians function handle of the nlmodel object.
        %
        %   USAGE
        %   SYS = linearize(M)
        %       gives linearization around zero equilibrum.
        %   SYS = linearize(M, x_e, u_e)
        %       gives linearization around specified equilibrium.
        %   
        %   It is checked if the given (or zero) equilirium actually is an
        %   equilibrium. Otherwise linearize throws an error.
        %
        %   Jacobians must be present!

            arguments
                M {mustBeOfClass(M,'nlmodel')}
                x_e (:, 1) double = 0
                u_e (:, 1) double = 0
            end

            % to create a linearization of the model, jacobians must be
            % available.
            if M.jacobiansPresent
                % if no equilibrium was given, set it to zero vectors of appropriate
                % dimensions. Otherwise check dimensions.
                if nargin == 1
                    x_e = zeros([M.nx 1]);
                    u_e = zeros([M.nu 1]);
    
                    % check if zero equilibrium actually is an equilibrium
                    if ~(abs(M.f(x_e, u_e)) < 1e4*eps)
                        error("MATLAB:nl_sys_toolbox:nlmodel:linearize:no_equilibrium", "Zero is no equilibrium of the given system. Specify an equilibrium in order to obtain a linearization!");
                    end
                elseif nargin == 3
                    if ~(length(x_e) == M.nx)
                        error("MATLAB:nl_sys_toolbox:linearize:state_offset_dim", "Specified equilibrium state x_e has incompatible length!");
                    end
                    if ~(length(u_e) == M.nu)
                        error("MATLAB:nl_sys_toolbox:linearize:input_offset_dim", "Specified equilibrium input u_e has incompatible length!");
                    end
    
                    % check if given (x_e, u_e) actually is an equilibrium
                    if ~(abs(M.f(x_e, u_e)) < 1e4*eps)
                        error("MATLAB:nl_sys_toolbox:linearize:no_equilibrium", "(x_e, u_e) must be an equilibrium point!");
                    end
                else
                    error("MATLAB:nl_sys_toolbox:linearize:argument_count", "nlmodel.linearize requires exactly none or two input arguments!");
                end

                % everything okay. Obtain system matrices and create
                % ss-object.
                [A, B, C, D] = M.jacobians(x_e, u_e);
                SYS = ss(A, B, C, D, "inputName", M.inputName, "stateName", M.stateName, "outputName", M.outputName);
            else
                error("MATLAB:nl_sys_toolbox:linearize:jac_missing", "Model does not include jacobians!");
            end
        end


        %% OPERATOR overloading functions
        function M = plus(M1, M2)
        %PLUS  Adds two nlmodel objects together.
        %
        %   M = PLUS(M1,M2) performs M = M1 + M2. This is equivalent 
        %   to connecting M1 and M2 in parallel.
        %   States are being stacked, state of first summand top, state
        %   of second summand bottom.

            arguments
                M1 {mustBeOfClass(M1,'nlmodel')}
                M2 {mustBeOfClass(M2,'nlmodel')}
            end

            % dimensions of addition
            nx_add = M1.nx + M2.nx;
            nu_add = M1.nu;
            ny_add = M1.ny;

            % Channel names of addition
            stateName_add = [M1.stateName; M2.stateName];

            % check input dimensions
            if ~(M1.nu == M2.nu)
                error("MATLAB:nl_sys_toolbox:nlmodel:plus:input_dim", "When adding systems, they must have the same input dimension!");
            end

            % check input names
            if isEmptyCharCell(M1.inputName)
                inputName_add = M2.inputName;
            elseif isEmptyCharCell(M2.inputName)
                inputName_add = M1.inputName;
            elseif ~isequal(M1.inputName, M2.inputName)
                warning("MATLAB:nl_sys_toolbox:nlmodel:plus:input_names", "Ignoring all input names because of name conflicts.");
                inputName_add = repmat({''}, nu_add, 1);
            else
                inputName_add = M1.inputName;
            end
            
            % check output dimensions
            if ~(M1.ny == M2.ny)
                error("MATLAB:nl_sys_toolbox:nlmodel:plus:output_dim", "When adding systems, they must have the same output dimension!");
            end

            % check output names
            if isEmptyCharCell(M1.outputName)
                outputName_add = M2.outputName;
            elseif isEmptyCharCell(M2.outputName)
                outputName_add = M1.outputName;
            elseif ~isequal(M1.outputName, M2.outputName)
                warning("MATLAB:nl_sys_toolbox:nlmodel:plus:output_names", "Ignoring all output names because of name conflicts.");
                outputName_add = repmat({''}, ny_add, 1);
            else
                outputName_add = M1.outputName;
            end

            % dynamics of addition
            f_add = @(x, u) [M1.f(x(1:M1.nx, :), u);
                             M2.f(x(M1.nx+1:M1.nx+M2.nx, :), u)];

            % output map of addition. Handle possible direct feed through.
            if M1.possibleDirectFeedThrough
                if M2.possibleDirectFeedThrough
                    h_add = @(x, u) M1.h(x(1:M1.nx, :), u) + M2.h(x(M1.nx+1:M1.nx+M2.nx, :), u);
                else
                    h_add = @(x, u) M1.h(x(1:M1.nx, :), u) + M2.h(x(M1.nx+1:M1.nx+M2.nx, :));
                end
            else
                if M2.possibleDirectFeedThrough
                    h_add = @(x, u) M1.h(x(1:M1.nx, :)) + M2.h(x(M1.nx+1:M1.nx+M2.nx, :), u);
                else
                    h_add = @(x) M1.h(x(1:M1.nx, :)) + M2.h(x(M1.nx+1:M1.nx+M2.nx, :));
                end
            end

            % only pass jacobians function if jacobians of both summands are
            % available. Function definitions are not allowed inside if
            % statements.
            if M1.jacobiansPresent && M2.jacobiansPresent
                jacobians_add = @(x, u) jacobians_plus_(M1, M2, x, u);
                M = nlmodel(f_add, h_add, nx_add, nu_add, ny_add, 'stateName', stateName_add, 'inputName', inputName_add, 'outputName', outputName_add, 'jacobians', jacobians_add);
            else
                M = nlmodel(f_add, h_add, nx_add, nu_add, ny_add, 'stateName', stateName_add, 'inputName', inputName_add, 'outputName', outputName_add);
            end
        end

        function M = mtimes(M1, M2)
        %MTIMES  Multiplies two nlmodel objects together.
        %
        %   M = MTIMES(M1,M2) performs the multiplication M = M1 * M2. 
        %   This is equivalent to connecting M1 and M2 in series as 
        %   follows:
        %
        %      u ----> M2 ----> M1 ----> y
        %   
        %   If one of them is scalar, it is expanded to appropriate
        %   dimension first. Because of that, order of multiplication
        %   matters, even if one of the arguments is scalar!

            arguments
                M1 {mustBeOfClass(M1,'nlmodel')}
                M2 {mustBeOfClass(M2,'nlmodel')}
            end

            % Expand scalars to static gains of appropriate dimension
            if isscalar(M1) && isnumeric(M1)
                M1 = ss2nlmodel(ss(M1*eye(M2.ny), 'inputName', M2.outputName, 'outputName', M2.outputName));
            elseif isscalar(M2) && isnumeric(M2)
                M2 = ss2nlmodel(ss(M2*eye(M1.nu), 'inputName', M1.inputName, 'outputName', M1.inputName));
            end

            % check if systems are compatible
            if ~(M2.ny == M1.nu)
                error("MATLAB:nl_sys_toolbox:nlmodel:mtimes:incompatible_dimensions", "Number of outputs of right system must equal number of inputs of the left system!");
            end
%             if ~isequal(M2.outputName, M1.inputName)
%                 warning("MATLAB:nl_sys_toolbox:nlmodel:mtimes:names", "Output names of right system do not match input names of left system. Proceeding...");
%             end

            % extract system data for better readability
            nx1 = M1.nx;
            f1 = @M1.f;
            h1 = @M1.h;
            
            nx2 = M2.nx;
            f2 = @M2.f;
            h2 = @M2.h;

            % dynamics of product. Handle possible direct feed through
            if M2.possibleDirectFeedThrough
                f_prod = @(x, u) [f1(x(1:nx1, :), h2(x(nx1+1:nx1+nx2, :), u));
                             f2(x(nx1+1:nx1+nx2, :), u)];
            else
                f_prod = @(x, u) [f1(x(1:nx1, :), h2(x(nx1+1:nx1+nx2, :)));
                             f2(x(nx1+1:nx1+nx2, :), u)];
            end

            % output map of addition. Handle possible direct feed through.
            if M1.possibleDirectFeedThrough
                if M2.possibleDirectFeedThrough
                    h_prod = @(x, u) h1(x(1:nx1, :), h2(x(nx1+1:nx1+nx2, :), u));
                else
                    h_prod = @(x) h1(x(1:nx1, :), h2(x(nx1+1:nx1+nx2, :)));
                end
            else
                h_prod = @(x) h1(x(1:nx1, :));
            end

            % channel names of product.
            stateName_prod = [M1.stateName; M2.stateName];
            inputName_prod = M2.inputName;
            outputName_prod = M1.outputName;
        
            % dimensions of product
            nx_prod = nx1 + nx2;
            nu_prod = M2.nu;
            ny_prod = M1.ny;
                    
            % only pass jacobians function if jacobians of both operands are
            % available. Function definitions are not allowed inside if
            % statements.
            if M1.jacobiansPresent && M2.jacobiansPresent
                jacobians_prod = @(x, u) jacobians_mtimes_(M1, M2, x, u);
                M = nlmodel(f_prod, h_prod, nx_prod, nu_prod, ny_prod, 'stateName', stateName_prod, 'inputName', inputName_prod, 'outputName', outputName_prod, 'jacobians', jacobians_prod);
            else
                M = nlmodel(f_prod, h_prod, nx_prod, nu_prod, ny_prod, 'stateName', stateName_prod, 'inputName', inputName_prod, 'outputName', outputName_prod);
            end
        end

        function M = vertcat(varargin)
        %VERTCAT  Vertical concatenation of nlmodel objects.
        %
        %   M = VERTCAT(M1,M2,...) performs the concatenation operation
        %   [M1 ; M2 ; ...].

            arguments (Repeating)
                varargin {mustBeOfClass(varargin,'nlmodel')}
            end

            ni = nargin;

            M = varargin{1};

            ignoringInputName = false;
            for j=2:ni
                Mj = varargin{j};

                % Check input dimensions
                if ~(M.nu == Mj.nu)
                    error("MATLAB:nl_sys_toolbox:nlmodel:vertcat:input_dimension", "Input dimension of all systems must be equal!");
                end

                % Check input names
                if ~ignoringInputName
                    if isEmptyCharCell(M.inputName)
                        inputName_vertcat = Mj.inputName;
                    elseif isEmptyCharCell(Mj.inputName)
                        inputName_vertcat = M.inputName;
                    elseif ~isequal(M.inputName, Mj.inputName)
                        warning("MATLAB:nl_sys_toolbox:nlmodel:vertcat:input_names", "Ignoring all input names because of name conflicts.");
                        inputName_vertcat = repmat({''}, M.nu, 1);
                        ignoringInputName = true;
                    else
                        inputName_vertcat = M.inputName;
                    end
                end

                % Dynamics of concatenation
                f_vertcat = @(x, u) [M.f(x(1:M.nx, :), u);
                                     Mj.f(x(M.nx+1:M.nx+Mj.nx, :), u)];

                % Output map of concatenation. Handle possible direct feed
                % through.
                if M.possibleDirectFeedThrough && Mj.possibleDirectFeedThrough
                    h_vertcat = @(x, u) [M.h(x(1:M.nx, :), u);
                                         Mj.h(x(M.nx+1:M.nx+Mj.nx, :), u)];
                elseif M.possibleDirectFeedThrough && ~Mj.possibleDirectFeedThrough
                    h_vertcat = @(x, u) [M.h(x(1:M.nx, :), u);
                                         Mj.h(x(M.nx+1:M.nx+Mj.nx, :))];                    
                elseif ~M.possibleDirectFeedThrough && Mj.possibleDirectFeedThrough
                    h_vertcat = @(x, u) [M.h(x(1:M.nx, :));
                                         Mj.h(x(M.nx+1:M.nx+Mj.nx, :), u)];                    
                else
                    h_vertcat = @(x) [M.h(x(1:M.nx, :));
                                      Mj.h(x(M.nx+1:M.nx+Mj.nx, :))];                    
                end

                % channel names
                stateName_vertcat = [M.stateName; Mj.stateName];
                outputName_vertcat = [M.outputName; Mj.outputName];

                % dimensions
                nx_vertcat = M.nx + Mj.nx;
                nu_vertcat = M.nu;
                ny_vertcat = M.ny + Mj.ny;

                % only pass jacobians function if jacobians of both operands are
                % available. Function definitions are not allowed inside if
                % statements.
                if M.jacobiansPresent && Mj.jacobiansPresent
                    jacobians_vertcat = @(x, u) jacobians_vertcat_(M, Mj, x, u);
                    M = nlmodel(f_vertcat, h_vertcat, nx_vertcat, nu_vertcat, ny_vertcat, 'stateName', stateName_vertcat, 'inputName', inputName_vertcat, 'outputName', outputName_vertcat, 'jacobians', jacobians_vertcat);
                else
                    M = nlmodel(f_vertcat, h_vertcat, nx_vertcat, nu_vertcat, ny_vertcat, 'stateName', stateName_vertcat, 'inputName', inputName_vertcat, 'outputName', outputName_vertcat);
                end

                % alternative, but leads to more complicated functions,
                % which might be less performant... and i don't have
                % jacobians incorporated in series_nl yet.
                % M = series_nl(Mj, M, 0);
            end
        end

        function M = horzcat(varargin)
        %HORZCAT  Horizontal concatenation of nlmodel objects.
        %
        %   M = HORZCAT(M1,M2,...) performs the concatenation operation
        %   [M1 , M2 , ...].

            arguments (Repeating)
                varargin {mustBeOfClass(varargin,'nlmodel')}
            end

            ni = nargin;

            M = varargin{1};

            ignoringOutputName = false;
            for j=2:ni
                Mj = varargin{j};

                % Check output dimensions
                if ~(M.ny == Mj.ny)
                    error("MATLAB:nl_sys_toolbox:nlmodel:horzcat:output_dimension", "Output dimension of all systems must be equal!");
                end

                % Check output names
                if ~ignoringOutputName
                    if isEmptyCharCell(M.outputName)
                        outputName_horzcat = Mj.outputName;
                    elseif isEmptyCharCell(Mj.outputName)
                        outputName_horzcat = M.outputName;
                    elseif ~isequal(M.outputName, Mj.outputName)
                        warning("MATLAB:nl_sys_toolbox:nlmodel:horzcat:output_names", "Ignoring all output names because of name conflicts.");
                        outputName_horzcat = repmat({''}, M.ny, 1);
                        ignoringOutputName = true;
                    else
                        outputName_horzcat = M.outputName;
                    end
                end

                % Dynamics of concatenation
                f_horzcat = @(x, u) [M.f(x(1:M.nx, :), u(1:M.nu, :));
                                     Mj.f(x(M.nx+1:M.nx+Mj.nx, :), u(M.nu+1:M.nu+Mj.nu, :))];

                % Output map of concatenation. Handle possible direct feed
                % through.
                if M.possibleDirectFeedThrough && Mj.possibleDirectFeedThrough
                    h_horzcat = @(x, u) M.h(x(1:M.nx, :), u(1:M.nu, :)) + Mj.h(x(M.nx+1:M.nx+Mj.nx, :), u(M.nu+1:M.nu+Mj.nu, :));
                elseif M.possibleDirectFeedThrough && ~Mj.possibleDirectFeedThrough
                    h_horzcat = @(x, u) M.h(x(1:M.nx, :), u(1:M.nu, :)) + Mj.h(x(M.nx+1:M.nx+Mj.nx, :));
                elseif ~M.possibleDirectFeedThrough && Mj.possibleDirectFeedThrough
                    h_horzcat = @(x, u) M.h(x(1:M.nx, :)) + Mj.h(x(M.nx+1:M.nx+Mj.nx, :), u(M.nu+1:M.nu+Mj.nu, :));
                else
                    h_horzcat = @(x) M.h(x(1:M.nx, :)) + Mj.h(x(M.nx+1:M.nx+Mj.nx, :));
                end

                % channel names
                stateName_horzcat = [M.stateName; Mj.stateName];
                inputName_horzcat = [M.inputName; Mj.inputName];

                % dimensions
                nx_horzcat = M.nx + Mj.nx;
                nu_horzcat = M.nu + Mj.nu;
                ny_horzcat = M.ny;

                % only pass jacobians function if jacobians of both operands are
                % available. Function definitions are not allowed inside if
                % statements.
                if M.jacobiansPresent && Mj.jacobiansPresent
                    jacobians_horzcat = @(x, u) jacobians_horzcat_(M, Mj, x, u);
                    M = nlmodel(f_horzcat, h_horzcat, nx_horzcat, nu_horzcat, ny_horzcat, 'stateName', stateName_horzcat, 'inputName', inputName_horzcat, 'outputName', outputName_horzcat, 'jacobians', jacobians_horzcat);
                else
                    M = nlmodel(f_horzcat, h_horzcat, nx_horzcat, nu_horzcat, ny_horzcat, 'stateName', stateName_horzcat, 'inputName', inputName_horzcat, 'outputName', outputName_horzcat);
                end
            end
        end

        function [NY, NU] = size(M)
        %SIZE  Size of input/output models.
        %
        %   S = SIZE(M) returns
        %      * S = [NY NU] for a single model M with NY outputs and NU inputs
            NY = M.ny;
            NU = M.nu;
        end

        %% PROPERTY Getters and Setters
        function M = set.f(M, f)
            % SET function for f property

            % check number of input arguments
            if nargin(f) ~= 2
                error("MATLAB:nl_sys_toolbox:nlmodel:dynamics:nargin", "f (model dynamics) is required to have exactly two arguments!");
            end
            M.f = f;
        end
        
        function M = set.h(M, h)
            % SET function for h property

            % check number of input arguments
            h_narg = nargin(h);
            if ~(h_narg == 1 || h_narg == 2)
                error("MATLAB:nl_sys_toolbox:nlmodel:output_map:nargin", "h (output map) is required to have one or two arguments!");
            end
            M.h = h;
        end

        function ret = get.possibleDirectFeedThrough(M)
            % GET function for possibleDirectFeedThrough property

            % if the output map h only has one argument (x), then no direct
            % feedthrough is possible
            if nargin(M.h) == 1
                ret = false;
            else
                ret = true;
            end
        end

        function ret = get.h_xu(M)
            % GET function for h_xu property
            if M.possibleDirectFeedThrough
                ret = @M.h;
            else
                ret = @(x, u) M.h(x);
            end
        end

        function ret = get.jacobians(M)
            % GET function for jacobians property
            if M.jacobiansPresent
                ret = @M.jacobians_;
            else
                error("MATLAB:nl_sys_toolbox:nlmodel:jacobians:not_available", "Model does not contain a jacobians function handle!");
            end
        end

        % stateName
        function M = set.stateName(M, stateName)
            % SET function for stateName property
            if ~length(stateName) == M.nx
                error("MATLAB:nl_sys_toolbox:nlmodel:stateName:size", "Number of provided stateNames does not equal number of state components!");
            else
                M.stateName_ = stateName;
            end
        end
        function stateName = get.stateName(M)
            % GET function for stateName property
            stateName = M.stateName_;
        end

        % inputName
        function M = set.inputName(M, inputName)
            % SET function for inputName property
            if ~length(inputName) == M.nu
                error("MATLAB:nl_sys_toolbox:nlmodel:inputName:size", "Number of provided inputNames does not equal number of input components!");
            else
                M.inputName_ = inputName;
            end
        end
        function inputName = get.inputName(M)
            % GET function for inputName property
            inputName = M.inputName_;
        end
        
        % outputName
        function M = set.outputName(M, outputName)
            % SET function for outputName property
            if ~length(outputName) == M.ny
                error("MATLAB:nl_sys_toolbox:nlmodel:outputName:size", "Number of provided outputNames does not equal number of output components!");
            else
                M.outputName_ = outputName;
            end
        end
        function outputName = get.outputName(M)
            % GET function for outputName property
            outputName = M.outputName_;
        end
    end

    %% CONSISTENCY checking methods (private)
    methods (Access=private)
        function checkDynamics(M)
            testState = zeros([M.nx, 1]);
            testInput = zeros([M.nu, 1]);

            % check dimension of function input arguments
            try
                f_eval = M.f(testState, testInput);
            catch ME
                msg = "Provided function handle f does not meet the requirements. It doesn't accept column vectors of dimensions nx and nu, or contains some other errors.";
                causeException = MException("MATLAB:nl_sys_toolbox:nlmodel:dynamics:argin", msg);
                ME = addCause(ME, causeException);
                rethrow(ME);
            end

            % check dimension of function output
            if ~isequal(size(f_eval), [M.nx, 1])
                error("MATLAB:nl_sys_toolbox:nlmodel:dynamics:argout", "Provided function handle f does not meet the requirements. It doesn't return a column vector of dimension nx.");
            end
        end

        function checkOutputMap(M)
            testState = zeros([M.nx, 1]);
            testInput = zeros([M.nu, 1]);

            % check dimension of function input arguments
            try
                h_narg = nargin(M.h);
                if h_narg == 1
                    h_eval = M.h(testState);
                elseif h_narg == 2
                    h_eval = M.h(testState, testInput);
                end
            catch ME
                msg = "Provided function handle h does not meet the requirements. It doesn't accept column vectors of dimensions nx (and nu), or contains some other errors.";
                causeException = MException("MATLAB:nl_sys_toolbox:nlmodel:output_map:argin", msg);
                ME = addCause(ME, causeException);
                rethrow(ME);
            end

            % check dimension of function output
            if ~isequal(size(h_eval), [M.ny, 1])
                error("MATLAB:nl_sys_toolbox:nlmodel:output_map:argout", "Provided function handle h does not meet the requirements. It doesn't return a column vector of dimension ny.");
            end
        end
        
        function checkJacobians(M)
            if M.jacobiansPresent
                testState = zeros([M.nx, 1]);
                testInput = zeros([M.nu, 1]);
    
                % try to evaluate function with test input arguments
                try
                    [A, B, C, D] = M.jacobians(testState, testInput);
                catch ME
                    msg = "Failed to evaluate provided jacobians function handle!";
                    causeException = MException("MATLAB:nl_sys_toolbox:nlmodel:jacobians:eval", msg);
                    ME = addCause(ME, causeException);
                    rethrow(ME);
                end
    
                % check dimension of jacobi matrices
                if ~isequal(size(A), [M.nx M.nx])
                    error("MATLAB:nl_sys_toolbox:nlmodel:jacobians:dim_A", "Dimension of jacobi-matrix A does not match [nx nx]!");
                elseif ~isequal(size(B), [M.nx M.nu])
                    error("MATLAB:nl_sys_toolbox:nlmodel:jacobians:dim_B", "Dimension of jacobi-matrix B does not match [nx nu]!");
                elseif ~isequal(size(C), [M.ny M.nx])
                    error("MATLAB:nl_sys_toolbox:nlmodel:jacobians:dim_C", "Dimension of jacobi-matrix C does not match [ny nx]!")
                elseif ~isequal(size(D), [M.ny M.nu])
                    error("MATLAB:nl_sys_toolbox:nlmodel:jacobians:dim_D", "Dimension of jacobi-matrix D does not match [ny nu]!");
                end
            end
        end

        function checkAll(M)
            M.checkDynamics();
            M.checkOutputMap();
            M.checkJacobians();
        end
    end
end

%% JACOBIANS (local functions)
% jacobians of addition
function [A, B, C, D] = jacobians_plus_(M1, M2, x, u)
    [A1, B1, C1, D1] = M1.jacobians(x(1:M1.nx, :), u);
    [A2, B2, C2, D2] = M2.jacobians(x(M1.nx+1:M1.nx+M2.nx, :), u);

    A = [A1 zeros([M1.nx M2.nx]); zeros([M2.nx M1.nx]) A2];
    B = [B1; B2];
    C = [C1 C2];
    D = D1 + D2;
end

% jacobians of multiplication
function [A, B, C, D] = jacobians_mtimes_(M1, M2, x, u)
    y2 = M2.h_xu(x(M1.nx+1:M1.nx+M2.nx, :), u);
    [A1, B1, C1, D1] = M1.jacobians(x(1:M1.nx, :), y2);
    [A2, B2, C2, D2] = M2.jacobians(x(M1.nx+1:M1.nx+M2.nx, :), u);

    A = [A1 B1*C2; zeros([M2.nx M1.nx]) A2];
    B = [B1*D2; B2];
    C = [C1 D1*C2];
    D = D1*D2;
end

% jacobians of vertical concatenation
function [A, B, C, D] = jacobians_vertcat_(M1, M2, x, u)
    [A1, B1, C1, D1] = M1.jacobians(x(1:M1.nx, :), u);
    [A2, B2, C2, D2] = M2.jacobians(x(M1.nx+1:M1.nx+M2.nx, :), u);

    A = [A1 zeros([M1.nx M2.nx]); zeros([M2.nx M1.nx]) A2];
    B = [B1; B2];
    C = [C1 zeros([M1.ny M2.nx]); zeros([M2.ny M1.nx]) C2];
    D = [D1; D2];
end

% jacobians of horizontal concatenation
function [A, B, C, D] = jacobians_horzcat_(M1, M2, x, u)
    [A1, B1, C1, D1] = M1.jacobians(x(1:M1.nx, :), u(1:M1.nu, :));
    [A2, B2, C2, D2] = M2.jacobians(x(M1.nx+1:M1.nx+M2.nx, :), u(M1.nu+1:M1.nu+M2.nu, :));

    A = [A1 zeros([M1.nx M2.nx]); zeros([M2.nx M1.nx]) A2];
    B = [B1 zeros([M1.nx M2.nu]); zeros([M2.nx M1.nu]) B2];
    C = [C1 C2];
    D = [D1 D2];
end

%% Custom validator functions
function mustBeOfClass(input,className)
    % Test for specific class name
    cname = class(input);
    if ~strcmp(cname,className)
        eid = 'Class:notCorrectClass';
        msg = ['Input must be of class ',className,'.'];
        throwAsCaller(MException(eid,msg))
    end
end

function ret = isEmptyCharCell(input)
    arguments
        input (:,1) cell
    end

    % Test if cell array contains only empty chars
    n = numel(input);
    ret = isequal(input, repmat({''}, n, 1));
end
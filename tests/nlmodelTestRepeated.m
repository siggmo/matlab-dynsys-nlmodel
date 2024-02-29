% automated test script for nlmodel class
% you can specify the to be tested operation, choose from
% ['addition' 'ss2nlmodelOnly' 'multiplication' 'vertcat' 'horzcat' 'upper_lft']
% as well as the number of test runs by setting the workspace variables
% testtype and n.
% After execution the test results can be found in the variables
% testResultNum and testResultBool. To examine the corresponding input
% data, see the variables testInputs and samplingInputs.

clearvars -except testtype n;

% number of runs
if ~exist('n', 'var')
    n = 300;
end
if ~exist('testtype', 'var')
    testtype = "addition";
end

% constants
modelTol = eps*10^5;
simTol = eps*10^8;
int_state = [0 30];
int_io = [0 30];

% store input and result data
testResultNum = cell(n, 6); % channelNamesOk, sysMatrixError, samplingErrorF, samplingErrorH, simulationErrorX, simulationErrorY
testResultBool = cell(n, 6); % are above values below tolerance?

% for progress display
reverseStr = '';

% perform test runs
for i=1:n    
    % perform operation on both ss and nlmodel objects
    switch testtype
        case "ss2nlmodelOnly"
            % generate test input data
            nx1 = randi(int_state);
            nu = randi(int_io);
            ny = randi(int_io);
            
            M1 = rss(nx1, ny, nu);

            M1.StateName = cellstr("M1_x" + string(1:nx1))';
            M1.InputName = cellstr("M1_u" + string(1:nu))';
            M1.OutputName = cellstr("M1_y" + string(1:ny))';

            testInput.nx1 = nx1;
            testInput.nu = nu;
            testInput.ny = ny;
            testInput.M1 = M1;
        
            testInputs(i) = testInput;

            try
                % create nlmodel object from ss object
                M1_nl = nlmodel(M1);

                % perform operation on both ss and nlmodel objects
                M_nl = M1_nl;
                M_lin = M1;
            catch E
                disp("Aktuelle Iteration: " + i);
                rethrow(E);
            end

        case "addition"
            % generate test input data
            nx1 = randi(int_state);
            nx2 = randi(int_state);
            nu = randi(int_io);
            ny = randi(int_io);
            
            M1 = rss(nx1, ny, nu);
            M2 = rss(nx2, ny, nu);

            M1.StateName = cellstr("M1_x" + string(1:nx1))';
            M1.InputName = cellstr("M1_u" + string(1:nu))';
            M1.OutputName = cellstr("M1_y" + string(1:ny))';
            M2.StateName = cellstr("M2_x" + string(1:nx2))';
            M2.InputName = cellstr("M2_u" + string(1:nu))';
            M2.OutputName = cellstr("M2_y" + string(1:ny))';

            testInput.nx1 = nx1;
            testInput.nx2 = nx2;
            testInput.nu = nu;
            testInput.ny = ny;
            testInput.M1 = M1;
            testInput.M2 = M2;
        
            testInputs(i) = testInput;

            try
                % create nlmodel objects from ss objects
                M1_nl = nlmodel(M1);
                M2_nl = nlmodel(M2);

                % perform operation on both ss and nlmodel objects
                M_nl = M1_nl + M2_nl;
                M_lin = M1 + M2;
            catch E
                disp("Aktuelle Iteration: " + i);
                rethrow(E);
            end

        case "multiplication"
            % generate test input data
            nx1 = randi(int_state);
            nx2 = randi(int_state);
            nu2 = randi(int_io);
            ny2 = randi(int_io);
            nu1 = ny2;
            ny1 = randi(int_io);
            
            M1 = rss(nx1, ny1, nu1);
            M2 = rss(nx2, ny2, nu2);

            M1.StateName = cellstr("M1_x" + string(1:nx1))';
            M1.InputName = cellstr("M1_u" + string(1:nu1))';
            M1.OutputName = cellstr("M1_y" + string(1:ny1))';
            M2.StateName = cellstr("M2_x" + string(1:nx2))';
            M2.InputName = cellstr("M2_u" + string(1:nu2))';
            M2.OutputName = cellstr("M2_y" + string(1:ny2))';

            testInput.nx1 = nx1;
            testInput.nx2 = nx2;
            testInput.nu1 = nu1;
            testInput.nu2 = nu2;
            testInput.ny1 = ny1;
            testInput.ny2 = ny2;
            testInput.M1 = M1;
            testInput.M2 = M2;
        
            testInputs(i) = testInput;

            try
                % create nlmodel objects from ss objects
                M1_nl = nlmodel(M1);
                M2_nl = nlmodel(M2);

                % perform operation on both ss and nlmodel objects
                M_nl = M1_nl * M2_nl;
                M_lin = M1 * M2;
            catch E
                disp("Aktuelle Iteration: " + i);
                rethrow(E);
            end

        case "vertcat"
            % generate test input data
            nx1 = randi(int_state);
            nx2 = randi(int_state);
            nu = randi(int_io);
            ny1 = randi(int_io);
            ny2 = randi(int_io);
            
            M1 = rss(nx1, ny1, nu);
            M2 = rss(nx2, ny2, nu);

            M1.StateName = cellstr("M1_x" + string(1:nx1))';
            M1.InputName = cellstr("M1_u" + string(1:nu))';
            M1.OutputName = cellstr("M1_y" + string(1:ny1))';
            M2.StateName = cellstr("M2_x" + string(1:nx2))';
            M2.InputName = cellstr("M2_u" + string(1:nu))';
            M2.OutputName = cellstr("M2_y" + string(1:ny2))';

            testInput.nx1 = nx1;
            testInput.nx2 = nx2;
            testInput.nu = nu;
            testInput.ny1 = ny1;
            testInput.ny2 = ny2;
            testInput.M1 = M1;
            testInput.M2 = M2;
        
            testInputs(i) = testInput;

            try
                % create nlmodel objects from ss objects
                M1_nl = nlmodel(M1);
                M2_nl = nlmodel(M2);

                % perform operation on both ss and nlmodel objects
                M_nl = [M1_nl; M2_nl];
                M_lin = [M1; M2];
            catch E
                disp("Aktuelle Iteration: " + i);
                rethrow(E);
            end

        case "horzcat"
            % generate test input data
            nx1 = randi(int_state);
            nx2 = randi(int_state);
            nu1 = randi(int_io);
            nu2 = randi(int_io);
            ny = randi(int_io);
            
            M1 = rss(nx1, ny, nu1);
            M2 = rss(nx2, ny, nu2);

            M1.StateName = cellstr("M1_x" + string(1:nx1))';
            M1.InputName = cellstr("M1_u" + string(1:nu1))';
            M1.OutputName = cellstr("M1_y" + string(1:ny))';
            M2.StateName = cellstr("M2_x" + string(1:nx2))';
            M2.InputName = cellstr("M2_u" + string(1:nu2))';
            M2.OutputName = cellstr("M2_y" + string(1:ny))';

            testInput.nx1 = nx1;
            testInput.nx2 = nx2;
            testInput.nu1 = nu1;
            testInput.nu2 = nu2;
            testInput.ny = ny;
            testInput.M1 = M1;
            testInput.M2 = M2;
        
            testInputs(i) = testInput;

            try
                % create nlmodel objects from ss objects
                M1_nl = nlmodel(M1);
                M2_nl = nlmodel(M2);

                % perform operation on both ss and nlmodel objects
                M_nl = [M1_nl, M2_nl];
                M_lin = [M1, M2];
            catch E
                disp("Aktuelle Iteration: " + i);
                rethrow(E);
            end

        case "upper_lft"
            % generate test input data
            nx1 = randi(int_state);
            nx2 = randi(int_state);
            nu = randi(int_io);
            ny = randi(int_io);
            nz = randi(int_io);
            nw = randi(int_io);

            M1 = rss(nx1, ny, nu);
            M1.D = M1.D*0;
            M2 = rss(nx2, nu+nz, ny+nw);

            M1.StateName = cellstr("M1_x" + string(1:nx1))';
            M1.InputName = cellstr("M1_u" + string(1:nu))';
            M1.OutputName = cellstr("M1_y" + string(1:ny))';
            M2.StateName = cellstr("M2_x" + string(1:nx2))';
            M2.InputName = cellstr("M2_u" + string(1:ny+nw))';
            M2.OutputName = cellstr("M2_y" + string(1:nu+nz))';

            testInput.nx1 = nx1;
            testInput.nx2 = nx2;
            testInput.nu = nu;
            testInput.ny = ny;
            testInput.nz = nz;
            testInput.nw = nw;
            testInput.M1 = M1;
            testInput.M2 = M2;

            testInputs(i) = testInput;

            try
                % create nlmodel objects from ss objects
                M1_nl = nlmodel(M1);

                % perform operation on both ss and nlmodel objects
                M_nl = nl_upper_lft(M1_nl, M2);
                M_lin = lft(M1, M2);
            catch E
                disp("Aktuelle Iteration: " + i);
                rethrow(E);
            end
    end

    % compare
    [testResult, samplingInput] = compareToLinearSS(M_nl, M_lin);
    channelNamesOk = testResult.channelNamesOk;
    sysMatrixError = testResult.sysMatrixError;
    samplingErrorF = testResult.samplingErrorF;
    samplingErrorH = testResult.samplingErrorH;
    simulationErrorX = testResult.simulationErrorX;
    simulationErrorY = testResult.simulationErrorY;
            
    testResultSys(i).M_nl = M_nl;
    testResultSys(i).M_lin = M_lin;
    testResultNum(i, :) = {double(channelNamesOk), sysMatrixError, samplingErrorF, samplingErrorH, simulationErrorX, simulationErrorY};
    testResultBool(i, :) = {channelNamesOk, sysMatrixError < modelTol, samplingErrorF < modelTol, samplingErrorH < modelTol, simulationErrorX < simTol, simulationErrorY < simTol};

    samplingInputs(i) = samplingInput;

    % Display the progress
    percentDone = 100 * i / n;
    msg = sprintf('Percent done: %3.1f', percentDone); %Don't forget this semicolon
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
end

testResultBool = cell2mat(testResultBool);
testResultNum = cell2mat(testResultNum);
% Dynamic State Estimation for Synchronous Generator Using Various Filters
clear;
close all; 
clc;

%% Start capturing all output to a text file (Ganti nama file)
diaryFileName = 'pfVariations.txt';

% Check if the diary file already exists and delete it if necessary
if exist(diaryFileName, 'file')
    delete(diaryFileName);
end

diary(diaryFileName);
diary on;

%% Parameter Initialization
parameter_ekf_ukf_pf; % Call to load parameters

%% Simulation Settings
iteration = 100;
x_range = [-0.1, 0.1];
num_of_state = 4;
num_of_measurement = 1;
threshold = 0.1; % Example threshold

% Prompt user to select filter methods
disp('Select filter methods:');
disp('1. Extended Kalman Filter (EKF)');
disp('2. Unscented Kalman Filter (UKF)');
disp('3. Particle Filter (PF) 100');
disp('4. Particle Filter (PF) 500');
disp('5. Particle Filter (PF) 1000');
selectedFiltersInput = input('Enter your choice (e.g., "1", "2 3", "1 and 2"): ', 's');

% Process input
selectedFilters = [];
if contains(selectedFiltersInput, '1')
    selectedFilters = [selectedFilters, 1];
end
if contains(selectedFiltersInput, '2')
    selectedFilters = [selectedFilters, 2];
end
if contains(selectedFiltersInput, '3')
    selectedFilters = [selectedFilters, 3];
end
if contains(selectedFiltersInput, '4')
    selectedFilters = [selectedFilters, 4];
end
if contains(selectedFiltersInput, '5')
    selectedFilters = [selectedFilters, 5];
end

% Check if any valid filter is selected
if isempty(selectedFilters)
    error('No valid filter selection. Please enter valid choices.');
end

% Generate deviations using systematic uniform random sampling
segment_size = diff(x_range) / iteration;
deviation = x_range(1) + (0:(iteration-1)) * segment_size + rand(1, iteration) * segment_size;

%% Simulation Iterations
for filterIndex = selectedFilters
    % Set filter method and initialize state based on user selection
    switch filterIndex
        case 1
            filterMethod = 'ekf';
            simModel = "project_simulation_ekf";
        case 2
            filterMethod = 'ukf';
            simModel = "project_simulation_ukf";
        case 3
            filterMethod = 'pf100';
            simModel = "project_simulation_pf100";
        case 4
            filterMethod = 'pf500';
            simModel = "project_simulation_pf500";
        case 5
            filterMethod = 'pf1000';
            simModel = "project_simulation_pf1000";
    end
    % Atur waktu simulasi
    load_system(simModel)
    set_param(simModel, 'StartTime', '0', 'StopTime', '6');

    % Preallocate 'out' structure for each filter method
    out = repmat(struct('itr', []), iteration, 1);

    % Run simulations for selected filter
    fprintf("\nRunning simulations using %s\n", filterMethod);
    for i = 1:iteration
        fprintf("Iteration %i of %i\n", i, iteration);
        fprintf('Deviation: %f\n', deviation(i));

        % Update state with deviation for simulation
        eval(['init_' filterMethod ' = init_state + deviation(i);']);

        % Run simulation
        tic;
        out(i).itr = sim(simModel);
        toc;
    end

    %% Convergence Analysis
    % Create a convergence analysis object with the collected data
    analysisTool = SignalAnalysisTool('Data', out);

    % Store deviations in the analysis tool
    analysisTool.Deviations = deviation;

    %% Perform MSE Error Analysis
    analysisTool.calculateMeanSquaredErrors();
    
    %% Perform MAPE Error Analysis
    analysisTool.calculateMeanAbsolutePercentageError();

    %% Check Convergence
    [isConvergentX, isConvergentY] = analysisTool.checkErrorConvergence(threshold);
    isConvergentX_total = sum(isConvergentX);
    fprintf("Number of convergent states for %s: %i\n", filterMethod, isConvergentX_total);
    isConvergentY_total = sum(isConvergentY);
    fprintf("Number of convergent measurements for %s: %i\n", filterMethod, isConvergentY_total);
    isConvergentXY_total = sum(isConvergentX & isConvergentY);
    fprintf("Number of convergent states and measurements for %s: %i\n", filterMethod, isConvergentXY_total);

    %% Results Plotting
    % Plot error under threshold.
    %analysisTool.plotErrorsBelowThreshold(threshold);
    
    % Plot deviation versus MSE.
    %analysisTool.plotDeviationVersusError(deviation);

    % Plot Monte Carlo results using methods in SignalAnalysisTool
    %analysisTool.plotMonteCarloSimulationResults('Threshold', false);
    %analysisTool.plotMonteCarloSimulationResults('Threshold', true);

    % Plot all MSE timesteps for all Monte Carlo results
    %analysisTool.plotMeanSquaredErrorTimesteps('Threshold', false);
    %analysisTool.plotMeanSquaredErrorTimesteps('Threshold', true);

    %% Store analysis results for each filter (Ganti nama file)
    eval([filterMethod '_pfVariations = analysisTool;'])
    name = [filterMethod '_pfVariations'];
    eval(['save(name, "' filterMethod '_pfVariations");']);
    eval(['clear ' filterMethod '_pfVariations;']);
end

% Stop capturing output and close the diary file
diary off;

% Optionally, display a message that the output has been saved
disp(['Simulation results have been saved to "', diaryFileName, '"']);

% main_plot;

function name = getVarName(~)
    name = inputname(1);
end
% Dynamic State Estimation for Synchronous Generator Using Various Filters
clear; close all; clc;

% Aktifkan beep
beep on;

%% Start capturing all output to a text file (Ganti nama file)
simulationFileName = 'wonoise';     % Nama file log (tanpa ekstensi, sesuaikan nama berdasarkan simulasi
simulationFileName = [simulationFileName '.txt'];

% Check if the diary file already exists and delete it if necessary
if exist(simulationFileName, 'file')
    delete(simulationFileName);
end

diary(simulationFileName);          % Mulai catat semua output ke file
diary on;

% -------------------------------
% Inisialisasi Parameter Sistem
% -------------------------------
parameter_ekf_ukf_pf;              % Panggil skrip parameter (init_state, Q, R, dsb)

% -------------------------------
% Pengaturan Simulasi
% -------------------------------
iteration = 10;
x_range = [-0.1, 0.1];
num_of_state = 4;
num_of_measurement = 1;
threshold = 0.1; % Example threshold

% -------------------------------
% Input Metode Filter dari Pengguna
% -------------------------------
disp('Select filter methods:');
disp('1. Extended Kalman Filter (EKF)');
disp('2. Unscented Kalman Filter (UKF)');
disp('3. Particle Filter (PF)');
selectedFiltersInput = input('Enter your choice (e.g., "1", "2 3", "1 and 2"): ', 's');

% Process input
selectedFilters = [];
if contains(selectedFiltersInput, '1') || contains(selectedFiltersInput, 'EKF')
    selectedFilters = [selectedFilters, 1];
end
if contains(selectedFiltersInput, '2') || contains(selectedFiltersInput, 'UKF')
    selectedFilters = [selectedFilters, 2];
end
if contains(selectedFiltersInput, '3') || contains(selectedFiltersInput, 'PF')
    selectedFilters = [selectedFilters, 3];
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
            simModel = "proj_sim_ekf_wonoise";
        case 2
            filterMethod = 'ukf';
            simModel = "proj_sim_ukf_wonoise";
        case 3
            filterMethod = 'pf';
            simModel = "proj_sim_pf_experiment";
    end
    % Atur waktu simulasi
    load_system(simModel)
    set_param(simModel, 'StartTime', '0', 'StopTime', '6');

    notConvergent = 1;
    p = 4;
    mean_errors_q4 = zeros(100,4);

    minus = 1;
    maxim = 20;
    range = maxim - minus;
    n = 100;
    delta = range/n;
    q4 = minus + delta *  (0:n);

    for j = 1:length(q4)
        % Preallocate 'out' structure for each filter method
        out = repmat(struct('itr', []), iteration, 1);

        % Run simulations for selected filter
        fprintf("\nRunning simulations using %s for j = %i\n", filterMethod, j);
        % disp(p)
        Q0 = diag([ 20 * 0.000158477870294882 ,  250 * 2.448576952494e-06 ,  6.5111462744571e-05 , q4(j) *  8.9480570589606e-05]);
        tic;
        for i = 1:iteration
            %fprintf("Iteration %i of %i\n", i, iteration);
            %fprintf('Deviation: %f\n', deviation(i));

            % Update state with deviation for simulation
            eval(['init_' filterMethod ' = init_state + deviation(i);']);
            % eval(['init_' filterMethod ' = init_state + 0.05;']);

            % Run simulation
            out(i).itr = sim(simModel);
        end
        toc;

        %% Convergence Analysis
        % Create a convergence analysis object with the collected data
        analysisTool = SignalAnalysisTool('Data', out);

        % Store deviations in the analysis tool
        analysisTool.Deviations = deviation;

        %% Perform MSE Error Analysis
        analysisTool.calculateMeanSquaredErrors();

        %% Perform MAPE Error Analysis
        analysisTool.calculateMeanAbsolutePercentageError();

        %% Calculate Statistic of Monte Carlo
        analysisTool.calculateStatOfMonteCarlo();

        %% Caluclate Mean Erros Statistics of Monte Carlo
        mean_errors_q4(j+1,:) = calculate_mean_errors(analysisTool);

        %% Check Convergence
        [isConvergentX, isConvergentY] = analysisTool.checkErrorConvergence(threshold);
        isConvergentX_total = sum(isConvergentX);
        fprintf("Number of convergent states for %s: %i\n", filterMethod, isConvergentX_total);
        isConvergentY_total = sum(isConvergentY);
        fprintf("Number of convergent measurements for %s: %i\n", filterMethod, isConvergentY_total);
        isConvergentXY_total = sum(isConvergentX & isConvergentY);
        fprintf("Number of convergent states and measurements for %s: %i\n", filterMethod, isConvergentXY_total);

        %if mean_errors  == 100
        %    notConvergent = 0;
        %end
        % p = p + 1;
    end
    beep;

    %% Results Plotting

    % Plot error under threshold.
    % analysisTool.plotErrorsBelowThreshold(threshold);

    % Plot deviation versus MSE.
    % analysisTool.plotDeviationVersusError(deviation);

    % Plot Monte Carlo results using methods in SignalAnalysisTool
    % analysisTool.plotMonteCarloSimulationResults('Threshold', false);
    % Buat calculate statistic of monte carlo. Jangan digabung sama
    % plotMontecarloSimulationResults.
    % analysisTool.plotMonteCarloSimulationResults('Threshold', true);

    % Plot all MSE timesteps for all Monte Carlo results
    % analysisTool.plotMeanSquaredErrorTimesteps('Threshold', false);
    % analysisTool.plotMeanSquaredErrorTimesteps('Threshold', true);

    % -------------------------------
    % Simpan Hasil Analisis
    % -------------------------------
    name_var = extractBefore(simulationFileName, '.txt');
    varName = [filterMethod '_' name_var];
    eval([varName ' = analysisTool;']);
    save([varName '.mat'], varName);
    eval(['clear' varName ';']);
end

% Stop capturing output and close the diary file
diary off;

% Optionally, display a message that the output has been saved
disp(['Simulation results have been saved to "', simulationFileName, '"']);

% generate_MC_paper_figures(name_var);
% compare_MCmean;

% Fungsi bantu (tidak dipakai di skrip utama)
function name = getVarName(~)
    name = inputname(1);
end
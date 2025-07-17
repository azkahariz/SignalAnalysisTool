classdef SignalAnalysisTool < handle
    % SIGNALANALYSISTOOL Performs dynamic state and measurement analysis.
    % This class provides functionalities for evaluating the accuracy and convergence 
    % of dynamic state and measurement estimation results from Monte Carlo simulations, 
    % especially in power systems. It supports statistical error analysis and various 
    % plotting capabilities.
    %
    % The analysis includes calculation of Mean Squared Error (MSE), Mean Absolute 
    % Percentage Error (MAPE), convergence checking based on user-defined thresholds, 
    % and statistical summaries (mean and standard deviation) over multiple iterations. 
    % The class also provides visualizations to support interpretation of state and 
    % measurement estimation performance.
    %
    % Properties:
    %   SimulationData - Structure to hold simulation data for states and measurements.
    %   ErrorStruct - Structure to hold error data for states and measurements.
    %   MSEStruct - Structure to hold Mean Squared Error (MSE) data for states and measurements.
    %   TimeToConvergence - Numeric array indicating the time it takes for each state to converge.
    %   NumOfState - Integer representing the number of states in the simulation data.
    %   NumOfMeasurement - Integer representing the number of measurements in the simulation data.
    %   MonteCarloIterations - Integer representing the number of iterations in the simulation data.
    %   Threshold - Numeric value used as a threshold for error to determine convergence.
    %   MAPE - Structure to hold Mean Absolute Percentage Error for states and measurements.
    %   Deviations - Array to hold deviation values corresponding to each iteration.
    %
    % Methods:
    %   SignalAnalysisTool - Constructor for creating an instance of the class.
    %   getArg - Method to retrieve argument values by name or return a default value.
    %   calculateMeanSquaredErrors - Method to calculate MSE for each state and measurement.
    %   checkErrorConvergence - Method to check if the simulation data has converged based on a specified threshold.
    %   plotErrorsBelowThreshold - Method to plot errors for states and measurements that are under a specified threshold.
    %   plotErrorsPerIteration - Helper method for plotting errors per iteration for a given set of data.
    %   plotMonteCarloSimulationResults - Method to plot results for each iteration of Monte Carlo simulations.
    %
    % Private Methods:
    %   storeSimulationData - Method to store simulation data into the SimulationData property.
    %   deriveFieldName - Helper method to determine the field name based on a given pattern.
    %   initializeErrorStructures - Method to initialize the structures used for error and convergence analysis.
    %   formatPlot - Helper method to format the plots with appropriate titles and labels.
    %   formatPlotAxes - Helper method to format the plot axes.
    %
    % Example:
    %   % Create an instance of the class with simulation data
    %   simData = load('simulationData.mat');
    %   analysisTool = SignalAnalysisTool(simData);
    %   analysisTool.calculateMeanSquaredErrors();
    %   [isConvergentX, isConvergentY] = analysisTool.checkErrorConvergence(0.01);
    %   analysisTool.plotErrorsBelowThreshold(0.05);
    %   analysisTool.plotMonteCarloSimulationResults('Threshold', true);
    %
    % See also HANDLE.

    properties
        SimulationData
        ErrorStruct
        MSEStruct
        TimeToConvergence
        NumOfState
        NumOfMeasurement
        MonteCarloIterations
        Threshold
        MAPE
        Deviations
        StatOfMonteCarlo
    end

    %% Constructor and initialization methods
    methods
        function obj = SignalAnalysisTool(varargin)
            % SIGNALANALYSISTOOL Constructs an instance of this class.
            % This constructor initializes an object of the SignalAnalysisTool class
            % with simulation data provided as input arguments. The data is processed
            % and stored within the object for further analysis.
            %
            % Syntax:
            %   obj = SignalAnalysisTool(data)
            %   obj = SignalAnalysisTool('data', data)
            %
            % Description:
            %   The constructor checks the input arguments for the simulation data.
            %   If the correct format and type of data are provided, it stores the data
            %   into the object's properties and initializes other necessary structures.
            %   It accepts data either directly or as a key-value pair with 'data' as the key.
            %
            % Inputs:
            %   varargin - A variable-length input argument list. It can contain the simulation
            %              data structure directly or as a key-value pair where 'data' is the key.
            %
            % Outputs:
            %   obj - An instance of the SignalAnalysisTool class with initialized properties.
            %
            % Errors:
            %   If no data is provided, or an incorrect number or format of inputs is given,
            %   the constructor will throw an error.
            %
            % Example:
            %   simData = load('simulationData.mat');
            %   analysisTool = SignalAnalysisTool(simData);
            %   % or
            %   analysisTool = SignalAnalysisTool('data', simData);

            % Cek jumlah argumen
            data = '';
            if isempty(varargin)
                error('ConvergenceAnalysis:NoData', 'Tidak ada data yang dimasukkan!');
            elseif length(varargin) > 2
                error('ConvergenceAnalysis:TooManyInputs', 'Terlalu banyak input! Harap masukkan hanya satu set data.');
            elseif length(varargin) == 2
                % Cek apakah format input adalah pasangan kunci-nilai
                if strcmpi(varargin{1}, 'data')
                    data = varargin{2};
                else
                    error('ConvergenceAnalysis:InvalidSyntax', 'Format input tidak valid. Gunakan ''data'', data.');
                end
            elseif length(varargin) == 1
                % Cek apakah input adalah struktur data
                if isstruct(varargin{1})
                    data = varargin{1};
                else
                    error('ConvergenceAnalysis:InvalidData', 'Input data harus berupa struktur.');
                end
            end

            %% Store data
            obj.storeSimulationData(data);
            obj.MonteCarloIterations = length(data);
            obj.NumOfState = size(obj.SimulationData.x.values{1},2);
            obj.NumOfMeasurement = size(obj.SimulationData.y.values{1},2);
            obj.initializeErrorStructures();
        end

        function value = getArg(obj, varargin, argName, defaultValue)
            % GETARG Retrieves the value of a named argument from a list, or returns a default value if the name is not found.
            % This method is useful for parsing function inputs that may or may not be present.
            %
            % Syntax:
            %   value = obj.getArg(varargin, argName, defaultValue)
            %
            % Description:
            %   It searches through 'varargin', which is a cell array of input arguments,
            %   looking for the argument named 'argName'. If 'argName' is found, it returns
            %   the subsequent value in the array. If 'argName' is not found, it returns
            %   'defaultValue'.
            %
            % Inputs:
            %   varargin - A cell array containing a list of arguments, typically passed to
            %              a function or method.
            %   argName - A string representing the name of the argument to search for.
            %   defaultValue - The value to return if 'argName' is not found in 'varargin'.
            %
            % Outputs:
            %   value - The value associated with 'argName' if found; otherwise, 'defaultValue'.
            %
            % Example:
            %   function myFunction(varargin)
            %       % Use the getArg method to retrieve an optional argument value
            %       paramValue = obj.getArg(varargin, 'paramName', defaultValue);
            %   end

            value = defaultValue;
            for i = 1:2:length(varargin)
                if strcmpi(varargin{i}, argName)
                    value = varargin{i + 1};
                    break;
                end
            end
        end
    end

    %% Data management and utility methods
    methods(Access = private)
        function obj = storeSimulationData(obj, data)
            % Menyimpan data simulasi ke dalam struktur SimulationData.
            % Parameter:
            %   data - Struktur data dari simulasi Simulink.

            fieldnames = data(1).itr.get; % Mendapatkan nama field dari struktur 'data'

            % Loop untuk memeriksa dan menyimpan setiap field
            for i = 1:length(fieldnames)
                fieldname = fieldnames{i};
                newFieldName = obj.deriveFieldName(fieldname);

                if ~isempty(newFieldName)
                    for j = 1:length(data)
                        obj.SimulationData.(newFieldName).values{j} = data(j).itr.(fieldname).signals.values;
                        obj.SimulationData.(newFieldName).time{j} = data(j).itr.(fieldname).time;
                    end
                end
            end
        end

        function newFieldName = deriveFieldName(obj, fieldname)
            % Menentukan newFieldName berdasarkan pola fieldname
            if startsWith(fieldname, 'x_hat') || startsWith(fieldname, 'y_hat')
                newFieldName = extractBefore(fieldname, "_hat") + "_hat";
            elseif startsWith(fieldname, 'x') || startsWith(fieldname, 'y')
                newFieldName = fieldname;
            else
                newFieldName = '';
            end
        end

        function obj = initializeErrorStructures(obj)
            % Menginisialisasi struktur untuk analisis kesalahan dan konvergensi.
            obj.MSEStruct.x = zeros(obj.MonteCarloIterations, obj.NumOfState);
            obj.MSEStruct.y = zeros(obj.MonteCarloIterations, obj.NumOfMeasurement);
            obj.ErrorStruct.x = cell(obj.MonteCarloIterations, 1);
            obj.ErrorStruct.y = cell(obj.MonteCarloIterations, 1);
            obj.TimeToConvergence.x = nan(obj.NumOfState, 1);
        end
    end

    %% Error analysis methods
    methods
        function obj = calculateMeanSquaredErrors(obj)
            % CALCULATEMEANSQUAREDERRORS Computes the Mean Squared Error
            % (MSE) for each state and measurement. This method calculates
            % the error for each iteration of the simulation data by
            % comparing the estimated states and measurements with the
            % actual values. The MSE is computed and stored within the
            % ErrorStruct property of the object.
            %
            % Syntax:
            %   obj.calculateMeanSquaredErrors()
            %
            % Description:
            %   For each state and measurement in the simulation data, this
            %   method iterates through all iterations and calculates the
            %   mean squared error. The results are aggregated in the
            %   ErrorStruct, which can be used for further analysis such as
            %   convergence checking and plotting error trends.
            %
            % Inputs:
            %   None - Operates on the properties of the SignalAnalysisTool
            %   instance.
            %
            % Outputs:
            %   obj - The current instance of the SignalAnalysisTool class
            %   with updated ErrorStruct
            %         properties containing the calculated errors.
            %
            % Example:
            %   % Assuming an instance 'analysisTool' of SignalAnalysisTool
            %   class has been created
            %   analysisTool.calculateMeanSquaredErrors(); % The
            %   ErrorStruct property of 'analysisTool' is now populated
            %   with MSE values.
            
            % Computes the MSE for each state and measurement.
            for i = 1:obj.MonteCarloIterations
                % Initializing temporary arrays for storing absolute errors
                tempErrorX = zeros(size(obj.SimulationData.x.values{i}));
                tempErrorY = zeros(size(obj.SimulationData.y.values{i}));

                for j = 1:obj.NumOfState
                    % Calculating absolute error for each state
                    tempErrorX(:, j) = abs(obj.SimulationData.x.values{i}(:, j) ...
                        - obj.SimulationData.x_hat.values{i}(:, j));
                    obj.MSEStruct.x(i, j) = mean(tempErrorX(:, j).^2); % Updated indexing
                end
                
                % Iterasi untuk setiap measurement
                for j = 1:obj.NumOfMeasurement
                    % Calculating absolute error for each measurement
                    tempErrorY(:, j) = abs(obj.SimulationData.y.values{i}(:, j) ...
                        - obj.SimulationData.y_hat.values{i}(:, j));
                    obj.MSEStruct.y(i, j) = mean(tempErrorY(:, j).^2); % Updated indexing
                end

                % Storing the absolute errors
                obj.ErrorStruct.x{i} = tempErrorX;
                obj.ErrorStruct.y{i} = tempErrorY;
            end
        end

        function calculateMeanAbsolutePercentageError(obj)
            % CALCULATEMEANABSOLUTEPERCENTAGEERROR Menghitung Mean Absolute Percentage Error (MAPE) untuk setiap state dan measurement.
            % MAPE dihitung sebagai rata-rata persentase kesalahan absolut antara nilai aktual (SimulationData)
            % dan nilai estimasi (ErrorStruct).
            %
            % Syntax:
            %   obj.calculateMeanAbsolutePercentageError()
            %
            % Description:
            %   Fungsi ini menghitung MAPE berdasarkan data yang tersimpan dalam properti SimulationData dan ErrorStruct.
            %   Hasil perhitungan disimpan dalam properti MAPE.

            % Inisialisasi MAPE untuk states dan measurements
            mapeStates = zeros(obj.NumOfState, obj.MonteCarloIterations);
            mapeMeasurements = zeros(obj.NumOfMeasurement, obj.MonteCarloIterations);

            % Hitung MAPE untuk setiap state pada setiap iterasi
            for j = 1:obj.MonteCarloIterations
                for i = 1:obj.NumOfState
                    actualState = obj.SimulationData.x.values{j}(:, i);
                    errorState = obj.ErrorStruct.x{j}(:, i);
                    validIndices = actualState ~= 0; % Menghindari pembagian dengan nol
                    mapeStates(i, j) = mean(abs(errorState(validIndices) ./ actualState(validIndices))) * 100;
                end

                for i = 1:obj.NumOfMeasurement
                    actualMeasurement = obj.SimulationData.y.values{j}(:, i);
                    errorMeasurement = obj.ErrorStruct.y{j}(:, i); % Gunakan error dari ErrorStruct
                    validIndices = actualMeasurement ~= 0; % Menghindari pembagian dengan nol
                    mapeMeasurements(i, j) = mean(abs(errorMeasurement(validIndices) ./ actualMeasurement(validIndices))) * 100;
                end
            end

            % Rata-rata MAPE untuk setiap state dan measurement
            avgMapeStates = mean(mapeStates, 2); % Rata-rata berdasarkan panjang sampel tiap iterasi
            avgMapeMeasurements = mean(mapeMeasurements, 2); % Rata-rata berdasarkan panjang sampel tiap iterasi

            % Simpan hasil MAPE ke dalam properti MAPE
            obj.MAPE.states = mapeStates;
            obj.MAPE.statesAvg = avgMapeStates;
            obj.MAPE.measurements = mapeMeasurements;
            obj.MAPE.measurementsAvg = avgMapeMeasurements;
        end
        
        function calculateStatOfMonteCarlo(obj, varargin)
            % CALCULATESTATOFMONTECARLO Menghitung nilai rata-rata dan
            % standar deviasi tiap state variable dan pengukuran dari
            % sejumlah simulasi Monte Carlo yang dilakukan. Dia dapat 
            % menghitung kondisi saat simulasi tersebut mengambil nilai
            % konvergen saja, atau semua nilai simulasi Monte Carlo.
            %
            % Syntax:
            %   obj.calculateteStatOfMonteCarlo()
            
            % Mengatur opsi default untuk 'Threshold'
            plotOnlyThreshold = false;

            % Mengecek dan mengatur opsi 'Threshold' dari varargin
            if ~isempty(varargin)
                if length(varargin) == 2 %% strcmpi(varargin{1}, 'Threshold')
                    plotOnlyThreshold = varargin{2};
                else
                    error('Invalid input arguments. Use ''Threshold'', true/false.');
                end
            end
            
            % Initialization statistic of monte carlo simulation
            obj.StatOfMonteCarlo.xMean = [];
            obj.StatOfMonteCarlo.yMean = [];
            obj.StatOfMonteCarlo.xStd = [];
            obj.StatOfMonteCarlo.yStd = [];

            % Cek konvergensi jika diperlukan
            if plotOnlyThreshold
                [isConvergentX, isConvergentY] = obj.checkErrorConvergence(obj.Threshold);
                indicesToPlot = find(isConvergentX & isConvergentY);
            else
                indicesToPlot = 1:obj.MonteCarloIterations;
            end
            
            %             numTimeSteps = size(obj.SimulationData.x_hat.values{1}, 1);
            %             numStateVars = obj.NumOfState;
            %             numMonteCarlo = length(indicesToPlot);
            % 
            %             data = zeros(numTimeSteps, numStateVars, numMonteCarlo); % Preallocation
            
            selectedCells = obj.SimulationData.x_hat.values(indicesToPlot);
            data = cat(3, selectedCells{:});
            obj.StatOfMonteCarlo.xMean = mean(data,3);
            obj.StatOfMonteCarlo.xStd = std(data,0,3);

            selectedCells = obj.SimulationData.y_hat.values(indicesToPlot);
            data = cat(3, selectedCells{:});
            obj.StatOfMonteCarlo.yMean = mean(data,3);
            obj.StatOfMonteCarlo.yStd = std(data,0,3);
        end

        function [isConvergentX, isConvergentY] = checkErrorConvergence(obj, threshold)
            % CHECKERRORCONVERGENCE Checks if the state and measurement errors have converged.
            % This method evaluates the convergence of each state and measurement by
            % comparing their error values against a specified threshold.
            %
            % Syntax:
            %   [isConvergentX, isConvergentY] = obj.checkErrorConvergence(threshold)
            %
            % Description:
            %   The method iterates through the error values stored in the ErrorStruct
            %   and compares each error value against the provided 'threshold'. If all
            %   error values for a particular state or measurement are below the threshold,
            %   it is considered to have converged. This method provides a boolean array
            %   indicating the convergence status for each state and measurement.
            %
            % Inputs:
            %   threshold - A numerical value specifying the maximum allowable error for
            %               convergence to be true.
            %
            % Outputs:
            %   isConvergentX - A boolean array indicating the convergence status of each state.
            %   isConvergentY - A boolean array indicating the convergence status of each measurement.
            %
            % Example:
            %   % Assuming an instance 'analysisTool' of SignalAnalysisTool class has been created
            %   % and 'calculateMeanSquaredErrors' has been called
            %   [isConvergentX, isConvergentY] = analysisTool.checkErrorConvergence(0.01);
            %   % The variables 'isConvergentX' and 'isConvergentY' contain boolean arrays
            %   % indicating the convergence of each state and measurement, respectively.

            obj.Threshold = threshold;
            isConvergentX = all(obj.MSEStruct.x < threshold, 2); % Periksa konvergensi state
            isConvergentY = all(obj.MSEStruct.y < threshold, 2); % Periksa konvergensi measurement
        end
    
    end

    %% Plotting methods
    methods
        function obj = plotErrorsBelowThreshold(obj, threshold)
            % PLOTERRORSBELOWTHRESHOLD Plots the errors for states and measurements that are under a specified threshold.
            % This method visualizes the errors that are below the convergence threshold across all iterations,
            % allowing for a quick assessment of how many and which iterations are within acceptable error limits.
            %
            % Syntax:
            %   obj.plotErrorsBelowThreshold(threshold)
            %
            % Description:
            %   Using the threshold provided, this method checks the convergence and plots the errors
            %   for each state and measurement that are below this threshold. It generates a subplot
            %   for states and another for measurements, providing a graphical representation of
            %   the error distribution and convergence behavior over iterations.
            %
            % Inputs:
            %   threshold - A numerical value that defines the upper limit for error values to be included in the plot.
            %
            % Outputs:
            %   None - The method produces plots as its output and does not modify the object or return values.
            %
            % Example:
            %   % Assuming an instance 'analysisTool' of SignalAnalysisTool class has been created
            %   % and 'calculateMeanSquaredErrors' has been called
            %   analysisTool.plotErrorsBelowThreshold(0.05);
            %   % This will generate plots for errors under the threshold value of 0.05.

            [isConvergentX, isConvergentY] = obj.checkErrorConvergence(threshold);

            figure;
            % Subplot untuk variabel x
            subplot(2, 1, 1);  % Dua baris, satu kolom, posisi pertama
            obj.plotErrorsPerIteration(isConvergentX, obj.MSEStruct.x, obj.NumOfState, 'State', 'MSE di bawah threshold unutk tiap State');

            % Subplot untuk variabel y
            subplot(2, 1, 2);  % Dua baris, satu kolom, posisi kedua
            obj.plotErrorsPerIteration(isConvergentY, obj.MSEStruct.y, obj.NumOfMeasurement, 'Measurement', 'MSE di bawah threshold unutk tiap Measurement');
        end

        function plotDeviationVersusError(obj, deviations)
            % PLOTDEVIATIONVERSUSERROR - Plots the errors against deviations for each state and measurement.
            % This method visualizes how deviations affect the errors for each state and measurement.
            %
            % Syntax:
            %   obj.plotDeviationVersusError(deviations)
            %
            % Description:
            %   This method takes an array of deviations and plots them against the calculated errors
            %   for each state and measurement. This helps in understanding the impact of deviations
            %   on the estimation errors.
            %
            % Inputs:
            %   deviations - An array of deviation values corresponding to each iteration.
            %
            % Outputs:
            %   None - The method produces plots as its output.
            %
            % Example:
            %   % Assuming an instance 'analysisTool' of SignalAnalysisTool class has been created
            %   % and 'calculateMeanSquaredErrors' has been called
            %   deviations = [0.1, 0.2, ...]; % Array of deviations
            %   analysisTool.plotDeviationVersusError(deviations);

            % Check if deviations array matches the number of iterations
            if length(deviations) ~= obj.MonteCarloIterations
                error('The length of deviations array must match the number of iterations.');
            end
            
            % Plotting errors against deviations for each state
            figure;
            for i = 1:obj.NumOfState
                subplot(3, 2, i);
                plot(deviations, obj.MSEStruct.x(:,i), '-o');
                title(['State ', num2str(i), ' Error vs Deviation']);
                xlabel('Deviation');
                ylabel('Error');
            end

            % Plotting errors against deviations for measurement
            subplot(3, 2, obj.NumOfState + 1);
            plot(deviations, obj.MSEStruct.y(:, :), '-o');
            title('Measurement Error vs Deviation');
            xlabel('Deviation');
            ylabel('Error');
        end

        function obj = plotErrorsPerIteration(obj, isConvergent, mseData, NumItems, itemName, plotTitle)
            % PLOTERROSPERITERATION Helper method to plot errors per iteration for a given set of data.
            % This method provides detailed plots for the error of each state or measurement
            % across different iterations, given the convergence status.
            %
            % Syntax:
            %   obj.plotErrorsPerIteration(isConvergent, ErrorStruct, NumItems, itemName, plotTitle)
            %
            % Description:
            %   This method takes in the convergence status array, the error structure, and details
            %   about the items to be plotted. It then generates a plot for each item, overlaying
            %   the error values of the iterations that have converged according
            %   to the provided convergence status. It creates a subplot for each item with a title and labels indicating the type of item and the nature of the error being plotted.
            %
            % Inputs:
            %   isConvergent - A boolean array indicating which iterations have converged.
            %   ErrorStruct - A structure containing the error data for each iteration.
            %   NumItems - The number of items (states or measurements) to plot.
            %   itemName - A string representing the name of the item type being plotted (e.g., 'State' or 'Measurement').
            %   plotTitle - A string representing the title to be used for the plots.
            %
            % Outputs:
            %   None - The method produces plots as its output and does not return any value.
            %
            % Example:
            %   % Assuming 'analysisTool' is an instance of SignalAnalysisTool with populated ErrorStruct
            %   isConvergent = [true false true]; % Example convergence status for three iterations
            %   analysisTool.plotErrorsPerIteration(isConvergent, analysisTool.ErrorStruct, ...
            %   analysisTool.NumOfState, 'State', 'Error per State across Iterations');
            %   % This will create plots for each state, showing errors for iterations that have converged.

            convergentIndices = find(isConvergent);

            if isempty(convergentIndices)
                disp(['Tidak ada iterasi yang konvergen di bawah threshold yang diberikan untuk ', lower(itemName), '.']);
                return;
            end

            for j = 1:NumItems
                errorsForItem = mseData(convergentIndices, j);
                % plot(convergentIndices, errorsForItem, '-o', 'DisplayName', [itemName ' ' num2str(j)]);
                bar(convergentIndices, errorsForItem, 'grouped', 'DisplayName',[itemName ' ' num2str(j)]); % Membaut diagram batang grup
                hold on; % Pertahankan plot saat iterasi melalui setiap item
            end
            
            hold off; % Lepaskan plot
            legend show;
            xlabel('Iterasi');
            ylabel('Error');
            title(plotTitle);
            grid on;
            grid minor;
        end

        function plotMonteCarloSimulationResults(obj, varargin)
            % PLOTMONTECARLOSIMULATIONRESULTS Plots results for each iteration of Monte Carlo simulations.
            % This method visualizes the state and measurement estimates for each Monte Carlo
            % iteration and compares them to the actual values. It allows for an assessment
            % of the estimation accuracy over multiple simulation runs.
            %
            % Syntax:
            %   obj.plotMonteCarloSimulationResults()
            %   obj.plotMonteCarloSimulationResults('Threshold', true)
            %
            % Description:
            %   The method can plot all Monte Carlo results or, if a threshold is provided and set to true,
            %   it will only plot the results for iterations where the measurement error is below the threshold.
            %   This can be useful to analyze only the successful iterations according to the set criteria of convergence.
            %   The plots include all state and measurement estimates and can be configured to show only those that meet the convergence criteria.
            %
            % Inputs:
            %   varargin - (Optional) A sequence of parameter/value pairs.
            %              Currently, the method checks for a 'Threshold' parameter which, if true,
            %              enables plotting only the iterations where convergence based on
            %			   the specified error threshold is achieved.
            %              If 'Threshold' is not provided or false, all iterations are plotted.
            %
            % Outputs:
            %   None - The method produces figures containing the plots as its output and does not return any value.
            %
            % Example:
            %   % Assuming 'analysisTool' is an instance of SignalAnalysisTool and you want to plot all results
            %   analysisTool.plotMonteCarloSimulationResults();
            %
            %   % To plot only the results that are below a certain error threshold
            %   analysisTool.plotMonteCarloSimulationResults('Threshold', true);
            %
            % Note:
            %   The method uses several internal properties of the object to access the necessary data for plotting.
            %   Make sure that the simulation data has been properly stored and the errors have been calculated
            %   before calling this method.

            % Mengatur opsi default
            plotOnlyThreshold = false;
            % Mengecek input varargin
            if ~isempty(varargin)
                if length(varargin) == 2 %% strcmpi(varargin{1}, 'Threshold')
                    plotOnlyThreshold = varargin{2};
                else
                    error('Invalid input arguments. Use ''Threshold'', true/false.');
                end
            end

            % Extract the data from the first iteration for reference
            % signals
            x_sig = obj.SimulationData.x.values{1};
            x_time = obj.SimulationData.x.time{1};
            y_sig = obj.SimulationData.y.values{1};
            y_time = obj.SimulationData.y.time{1};

            % Define titles, labels, and display name for plots
            titles = {'x1: Rotor Angle Estimation ($\delta$)', 'x2: Rotor Speed Estimation ($\Delta\omega$)', ...
                'x3: Transient voltage -q axis Estimation ($e''_{q}$)', 'x4: Transient voltage -d axis Estimation ($e`_{d}$)', ...
                'y: Output Power (Pt)'};
            yLabels = {'$Elec. Rad (pu)$', '$Elec. Rad/s (pu)$', '$Voltage (pu)$', '$Voltage (pu)$', '$Active Power (pu)$'};
            dispNames = {'$\hat{x}_1$', '$x_1$', '$\hat{x}_2$', '$x_2$', '$\hat{x}_3$', '$x_3$', '$\hat{x}_4$', '$x_4$', '$\hat{y}$'  , '$y$'};

            % Cek konvergensi jika diperlukan
            if plotOnlyThreshold
                [isConvergentX, isConvergentY] = obj.checkErrorConvergence(obj.Threshold);
                indicesToPlot = find(isConvergentX == 1 & isConvergentY == 1);
            else
                indicesToPlot = 1:obj.MonteCarloIterations;
            end
            
            % Plotting figure
            figure;
            for i = 1:5
                subplot(3,2,i);
                handles = []; % Initialize handles

                %Plot sinyal aktual setelah memplot hasil simulasi Monte
                %Carlo
                % Plot sinyal aktual
                if i <= 4 % For states x1 to x4
                    signalTrue = plot(x_time, x_sig(:,i), 'b', 'linewidth', 3, 'DisplayName', dispNames{2*i});
                    hold on
                else % For measurement y
                    signalTrue = plot(y_time, y_sig, 'b', 'linewidth', 3, 'DisplayName', dispNames{2*i});
                    hold on
                end
                handles = [handles; signalTrue];
                if i <=4
                    widthAxes = max(x_sig(:,i)) - min(x_sig(:,i)) ;
                    axis([-Inf, inf, min(x_sig(:,i)) - widthAxes, max(x_sig(:,i)) + widthAxes ])
                else 
                    widthAxes = max(y_sig) - min(y_sig) ;
                    axis([-Inf, Inf, min(y_sig) - widthAxes, max(y_sig) + widthAxes ])
                end
                
                % Plot hasil simulasi Monte Carlo
                % [rows, cols] = size(x_sig);
                % x_mc_mean = [];
                for jIndex = 1:length(indicesToPlot)
                    j = indicesToPlot(jIndex); % Use a scalar index
                    % Pilih data yang sesuai
                    if i <= 4 % Untutk state x1 sampai x4
                        data = obj.SimulationData.x_hat.values{j}(:,i);
                        plotTime = obj.SimulationData.x_hat.time{j};
                    else % Untuk measurement y
                        data = obj.SimulationData.y_hat.values{j};
                        plotTime = obj.SimulationData.y_hat.time{j};
                    end
                    % x_mc_mean = [x_mc_mean; data'];

                    % Plot data for all iterations
                    plotMC = plot(plotTime, data, 'Color', [0, 1, 0, 0.1], 'linewidth', 0.1);
                    hold on;
                    if j == indicesToPlot(end)
                        set(plotMC, 'DisplayName', dispNames{2*i-1});
                        handles = [handles; plotMC];
                    end
                end

                obj.calculateStatOfMonteCarlo('Threshold', plotOnlyThreshold)

                % Rata-rata dan std dari tiap step pd simulasi monte carlo
                % xMean = mean(x_mc_mean,1);
                % xStd  = std(x_mc_mean,0,1);
                %
                % Mencari 3*xstd atas dan bawah
                % upperBound = (xMean + 3*xStd);
                % lowerBound = (xMean - 3*xStd);

                if i <=4 
                    upperBound = obj.StatOfMonteCarlo.xMean(:,i)' + 3*obj.StatOfMonteCarlo.xStd(:,i)';
                    lowerBound = obj.StatOfMonteCarlo.xMean(:,i)' - 3*obj.StatOfMonteCarlo.xStd(:,i)';
                else
                    upperBound = obj.StatOfMonteCarlo.yMean' + 3*obj.StatOfMonteCarlo.yStd';
                    lowerBound = obj.StatOfMonteCarlo.yMean' - 3*obj.StatOfMonteCarlo.yStd';
                end
                if i <= 4
                    plotMean = plot(plotTime, obj.StatOfMonteCarlo.xMean(:,i), 'LineWidth', 2, 'Color', 'Red','LineStyle','--');
                    namePlotMean = strcat( '$\mu_{\hat{x}_', '{', num2str(i) ,'}}$');
                    set(plotMean, 'DisplayName', namePlotMean);
                    handles = [handles; plotMean];
    
                    plotStd = fill([plotTime' fliplr(plotTime')], [upperBound fliplr(lowerBound)], 'yellow', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
                    namePlotStd = strcat('$\bar{\hat{x}}_',num2str(i),' \pm 3 \, \sigma_{\hat{x}_',num2str(i),'}$');
                    set(plotStd, 'DisplayName', namePlotStd);
                    handles = [handles; plotStd];
                else
                    plotMean = plot(plotTime, obj.StatOfMonteCarlo.yMean, 'LineWidth', 2, 'Color', 'Red','LineStyle','--');
                    namePlotMean = strcat('$\mu_{\hat{y}}$');
                    set(plotMean, 'DisplayName', namePlotMean);
                    handles = [handles; plotMean];
    
                    plotStd = fill([plotTime' fliplr(plotTime')], [upperBound(1,:) fliplr(lowerBound(1,:))], 'yellow', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
                    namePlotStd = strcat('$\bar{\hat{y}} \pm 3 \, \sigma_{\hat{y}}$');
                    set(plotStd, 'DisplayName', namePlotStd);
                    handles = [handles; plotStd];
                end


                % Plot sinyal aktual setelah memplot hasil simulasi Monte
                % Carlo
                % % Plot sinyal aktual
                % if i <= 4 % For states x1 to x4
                %     signalHandle = plot(x_time, x_sig(:,i), 'b', 'linewidth', 2, 'DisplayName', dispNames{2*i});
                % else % For measurement y
                %     signalHandle = plot(y_time, y_sig, 'b', 'linewidth', 2, 'DisplayName', dispNames{2*i});
                % end
                % handles(end+1) = signalHandle;


                uistack(signalTrue,'top')
                uistack(plotMean,'top')
                % Hanya panggil legend jika ada handle yang valid
                if ~isempty(handles)
                    legend(handles, 'location', 'southeast', 'interpreter', 'latex');
                end

                obj.formatPlot(titles{i}, yLabels{i});
                hold off;
            end
            set(gcf, 'PaperPositionMOde', 'auto');
        end

        function plotMeanSquaredErrorTimesteps(obj, varargin)
            % PLOTMEANSQUAREDERRORTIMESTEPS - Plots the Mean Squared Error (MSE) for each state and measurement
            % over time for each Monte Carlo simulation, with an option to plot only those simulations
            % that meet a convergence threshold.
            %
            % This method visualizes the MSE for each state and measurement across all timesteps
            % for each Monte Carlo iteration. It allows for an optional argument to plot only
            % the iterations that are below a specified convergence threshold.
            %
            % Syntax:
            %   obj.plotMeanSquaredErrorTimesteps()
            %   obj.plotMeanSquaredErrorTimesteps('Threshold', true/false)
            %
            % Inputs:
            %   varargin - (Optional) A sequence of parameter/value pairs.
            %              'Threshold' - A boolean indicating whether to plot only the simulations
            %                            that meet the convergence threshold. If true, only iterations
            %                            with errors below the object's Threshold property are plotted.
            %                            If false or omitted, all iterations are plotted.
            %
            % Outputs:
            %   None - The method produces a figure with subplots as its output.
            %
            % Example:
            %   % Example 1: Plot MSE for all Monte Carlo iterations
            %   obj.plotMeanSquaredErrorTimesteps();
            %
            %   % Example 2: Plot MSE only for iterations that are below the convergence threshold
            %   obj.plotMeanSquaredErrorTimesteps('Threshold', true);
            %
            % Notes:
            %   - The method uses several internal properties of the object to access the necessary data for plotting.
            %   - Ensure that the simulation data has been properly stored and the errors have been calculated
            %     before calling this method.
            %   - The plots include all state and measurement estimates and can be configured to show only those
            %     that meet the convergence criteria.


            % Mengatur opsi default untuk 'Threshold'
            plotOnlyThreshold = false;

            % Mengecek dan mengatur opsi 'Threshold' dari varargin
            if ~isempty(varargin)
                if strcmpi(varargin{1}, 'Threshold') && length(varargin) == 2
                    plotOnlyThreshold = varargin{2};
                else
                    error('Invalid input arguments. Use ''Threshold'', true/false.');
                end
            end

            % Mendefinisikan judul dan label untuk plot
            titles = {'x1: Rotor Angle Estimation ($\delta$)', 'x2: Rotor Speed Estimation ($\Delta\omega$)', ...
                'x3: Transient voltage -q axis Estimation ($e''_{q}$)', 'x4: Transient voltage -d axis Estimation ($e`_{d}$)', ...
                'y: Output Power (Pt)'};
            yLabels = {'$Elec. Rad (pu)$', '$Elec. Rad/s (pu)$', '$Voltage (pu)$', '$Voltage (pu)$', '$Active Power (pu)$'};

            % Cek konvergensi jika opsi 'Threshold' diaktifkan
            if plotOnlyThreshold
                [~, isConvergentY] = obj.checkErrorConvergence(obj.Threshold);
                indicesToPlot = find(isConvergentY);
            else
                indicesToPlot = 1:obj.MonteCarloIterations;
            end

            % Mempersiapkan figure
            figure;

            % Plot MSE untuk setiap state pada setiap iterasi Monte Carlo
            for i = 1:obj.NumOfState
                subplot(3, 2, i);
                hold on;
                for jIndex = 1:length(indicesToPlot)
                    j = indicesToPlot(jIndex);  
                    plot(obj.SimulationData.x_hat.time{j}, obj.ErrorStruct.x{j}(:, i), 'LineWidth', 0.1);
                end
                hold off;
                obj.formatPlot(titles{i}, yLabels{i});
            end

            % Plot MSE untuk measurement pada setiap iterasi Monte Carlo
            subplot(3, 2, 5);
            hold on;
            for jIndex = 1:length(indicesToPlot)
                j = indicesToPlot(jIndex);
                plot(obj.SimulationData.y_hat.time{j}, obj.ErrorStruct.y{j}, 'LineWidth', 0.1);
            end
            hold off;
            obj.formatPlot(titles{5}, yLabels{5});
        end

        function plotMeanAbsolutePercentageErrorTimesteps(obj)
            % PLOTMEANABSOLUTEPERCENTAGEERRORTIMESTEPS - Plots the Mean Absolute Percentage Error (MAPE) for each state and measurement
            % over time for each Monte Carlo simulation.
            %
            % Syntax:
            %   obj.plotMeanAbsolutePercentageErrorTimesteps()
            %
            % Description:
            %   This method visualizes the MAPE for each state and measurement across all timesteps
            %   for each Monte Carlo iteration.

            % Mendefinisikan judul dan label untuk plot
            titles = {'x1: Rotor Angle MAPE ($\delta$)', 'x2: Rotor Speed MAPE ($\Delta\omega$)', ...
                'x3: Transient voltage -q axis MAPE ($e''_{q}$)', 'x4: Transient voltage -d axis MAPE ($e`_{d}$)', ...
                'y: Output Power MAPE (Pt)'};
            yLabels = {'% Error', '% Error', '% Error', '% Error', '% Error'};

            % Mempersiapkan figure
            figure;

            % Plot MAPE untuk setiap state pada setiap iterasi Monte Carlo
            for i = 1:obj.NumOfState
                subplot(3, 2, i);
                hold on;
                for j = 1:obj.MonteCarloIterations
                    plot(obj.SimulationData.x.time{j}, obj.MAPE.states(i, j), 'LineWidth', 0.1);
                end
                hold off;
                obj.formatPlot(titles{i}, yLabels{i});
            end

            % Plot MAPE untuk measurement pada setiap iterasi Monte Carlo
            subplot(3, 2, obj.NumOfState + 1);
            hold on;
            for j = 1:obj.MonteCarloIterations
                plot(obj.SimulationData.y.time{j}, obj.MAPE.measurements(:, j), 'LineWidth', 0.1);
            end
            hold off;
            obj.formatPlot(titles{obj.NumOfState + 1}, yLabels{obj.NumOfState + 1});
        end
    end
    
    %% Additional helper methods for plotting
    methods(Access = private)
        function formatPlot(obj, plot_title, ylabelText)
            % FORMATPLOT Helper method to format the plots with appropriate titles and labels.
            % This method applies consistent formatting to the plots created by other methods
            % in the SignalAnalysisTool class. It sets the title and y-axis label with LaTeX
            % interpreter for better visual appearance.
            %
            % Syntax:
            %   obj.formatPlot(plot_title, ylabelText)
            %
            % Description:
            %   This private method is called internally to set the titles and y-axis labels
            %   of the plots generated by the class's public methods. It uses LaTeX formatting
            %   for the title and y-label to enhance the readability and presentation of plots.
            %
            % Inputs:
            %   plot_title - A string containing the title for the plot, which will be interpreted in LaTeX.
            %   ylabelText - A string containing the y-axis label for the plot, which will be interpreted in LaTeX.
            %
            % Outputs:
            %   None - The method does not return any value.
            %
            % Example Usage:
            %   % This method is used internally and is not intended to be called directly by users.
            %   % However, if needed to be called explicitly, here is an example:
            %   obj.formatPlot('State Estimation Error', '$\text{Error (units)}$');
            %
            % Note:
            %   This method is part of the private API of the SignalAnalysisTool class and is not intended
            %   for public use. It is documented here for completeness and for understanding the internals
            %   of the class's operation.

            title(plot_title, 'Interpreter', 'latex');
            xlabel('$Time(s)$', 'Interpreter', 'latex');
            ylabel(ylabelText, 'Interpreter', 'latex');
            obj.formatPlotAxes();
        end

        function formatPlotAxes(~)
            % FORMATPLOTAXES Helper method to format the axes of the plot.
            % This method sets consistent styling for the plot axes, enhancing the
            % visual appearance and readability. It is called internally to apply
            % standard formatting across various plots generated by the class.
            %
            % Syntax:
            %   obj.formatPlotAxes()
            %
            % Description:
            %   The method adjusts various properties of the axes, such as box style,
            %   tick direction, minor ticks, grid lines, and color settings. This ensures
            %   that all plots generated by the class have a uniform and professional look.
            %
            % Inputs:
            %   None - This method does not take any input arguments.
            %
            % Outputs:
            %   None - The method does not return any value.
            %
            % Example Usage:
            %   % This method is used internally and is not intended to be called directly by users.
            %   % However, if needed to be called explicitly for custom plots, here is an example:
            %   obj.formatPlotAxes();
            %
            % Note:
            %   This method is a part of the private API of the SignalAnalysisTool class and is not intended
            %   for public use. It is documented here for completeness and understanding the class's functionality.

            set(gca, 'Box', 'on', 'TickDir', 'out', 'TickLength', [.007 .007], ...
                'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
                'XGrid', 'on', 'XColor', [0.06 0.06 0.03], 'YColor', [0.06 0.06 0.03]);
            hold off;
        end
    end
end
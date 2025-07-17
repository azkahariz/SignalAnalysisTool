function compare_MCmean(name_var)
% Clear workspace, command window, and figures
clc;
close all;

% Load data from .mat files
load(['ekf_' name_var '.mat']);
load(['ukf_' name_var '.mat']);
load(['pf_' name_var '.mat']);

% Retrieve and display detailed information about variables in the workspace
variableInfo = whos; % Get details of all variables currently in the workspace
variableNames = {variableInfo.name}; % Extract the names of the variables into a cell array

disp('Variables in the loaded .mat files:');
disp(variableNames); % Print the variable names for verification

% Dynamically assign loaded variables to specific names (ekf, ukf, pf)
ekfVarName = variableNames{contains(variableNames, "ekf", 'IgnoreCase', true)};
ukfVarName = variableNames{contains(variableNames, "ukf", 'IgnoreCase', true)};
pfVarName  = variableNames{contains(variableNames, "pf", 'IgnoreCase', true)};

% Assign variables to meaningful names
ekf = eval(ekfVarName);
ukf = eval(ukfVarName);
pf  = eval(pfVarName);

% Extract true state and calculate errors for each filter
x_true = pf.SimulationData.x.values{1};

error_ekf = abs(ekf.StatOfMonteCarlo.xMean - x_true);
error_ukf = abs(ukf.StatOfMonteCarlo.xMean - x_true);
error_pf  = abs(pf.StatOfMonteCarlo.xMean - x_true);

% Compute mean error for each filter
mean_ekf = mean(error_ekf);
mean_ukf = mean(error_ukf);
mean_pf  = mean(error_pf);

std_ekf = std(error_ekf);
std_ukf = std(error_ukf);
std_pf  = std(error_pf);

% Combine mean errors into a matrix for plotting
mean_errors = [mean_ekf; mean_ukf; mean_pf];

% Number of states to compare
num_states = size(mean_errors, 2);

% Plot comparison for each state
for i = 1:num_states
    figure; % Create a new figure for each state
    bar(mean_errors(:, i), 'grouped'); % Bar chart for the current state
    grid on; % Add grid
    xlabel('Filter Method'); % Label x-axis
    ylabel('Mean Error'); % Label y-axis
    title(['Comparison for State ', num2str(i)]); % Add title
    xticks(1:3); % Set x-axis tick positions
    xticklabels({'EKF', 'UKF', 'PF'}); % Set x-axis tick labels
end

end
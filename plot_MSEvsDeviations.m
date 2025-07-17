clc;
%clear;
%close all;

%% Load the required .mat files for EKF, UKF, and PF results
% load('ResultMatFile\ts100\ekf_ts100.mat'); % Load EKF results
% load('ResultMatFile\ts100\ukf_ts100.mat'); % Load UKF results
% load('ResultMatFile\ts100\pf_ts100.mat');  % Load PF results

%% Retrieve detailed information about variables in the workspace
variableInfo = whos; % Get details of all variables currently in the workspace
variableNames = {variableInfo.name}; % Extract the names of the variables into a cell array

%% Display the names of the variables in the workspace
disp('Variables in the loaded .mat files:');
disp(variableNames); % Print the variable names for verification

%% Dynamically assign the loaded variables to specific names (ekf, ukf, pf)
% Assuming the variables are in the same order as the files loaded
variableNames{contains(variableNames, "ekf", 'IgnoreCase', true)}
ekf = eval(variableNames{contains(variableNames, "ekf", 'IgnoreCase', true)}); % Assign the first variable to 'ekf'
variableNames{contains(variableNames, "ukf", 'IgnoreCase', true)}
ukf = eval(variableNames{contains(variableNames, "ukf", 'IgnoreCase', true)}); % Assign the second variable to 'ukf'
variableNames{contains(variableNames, "pf", 'IgnoreCase', true)}
pf  = eval(variableNames{contains(variableNames, "pf", 'IgnoreCase', true)}); % Assign the third variable to 'pf'

% Optional: Clear unnecessary variables if needed
% clear variableInfo variableNames

%% Titles state space
titles = {'x1: Rotor Angle Estimation ($\delta$)',...
          'x2: Rotor Speed Estimation ($\Delta\omega$)', ...
          'x3: Transient voltage -q axis Estimation ($e''_{q}$)',...
          'x4: Transient voltage -d axis Estimation ($e`_{d}$)', ...
          'y: Output Power (Pt)'};
%% Save data to logarithmic function
logMSEStructEKF = [];
logMSEStructUKF = [];
logMSEStructPF  = [];

for i = 1:4
    logMSEStructEKF = [logMSEStructEKF log10(ekf.MSEStruct.x(:,i))];
    logMSEStructUKF = [logMSEStructUKF log10(ukf.MSEStruct.x(:,i))];
    logMSEStructPF  = [logMSEStructPF  log10(pf.MSEStruct.x(:,i))];
end

%% Plotting the data
fig = figure('Position', [100, 100, 1200, 600]); % High-quality figure size
t = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact'); % Compact spacing

for j = 1:4
    nexttile; % Move to the next subplot
    hold on;
    grid on
    grid minor;

    % Plot data
    plotEKF = plot(ekf.Deviations, logMSEStructEKF(:,j),'Color','red','Marker','o',...
        'DisplayName', 'EKF Method');
    plotUKF = plot(ukf.Deviations, logMSEStructUKF(:,j),'Color','blue','Marker','square',...
        'DisplayName','UKF Method');
    plotPF = plot(pf.Deviations,  logMSEStructPF(:,j) ,'Color','green','Marker','v',...
        'DisplayName','PF Method');
    title(titles{j}, 'Interpreter', 'latex');
    if j ==2
        plotLegend = [plotEKF, plotUKF, plotPF];
    end
    % Axis labels
    xlabel('Initial Deviation', 'Interpreter', 'latex');
    ylabel('$\log_{10}(MSE)$', 'Interpreter', 'latex');
end

% Add a global title
% sgtitle('Comparison of MSE for EKF, UKF, and PF Methods', 'FontSize', 13, 'Interpreter', 'latex');

% Add a single legend in the last tile
legend(plotLegend, 'Location', 'northeastoutside', 'Interpreter', 'latex', ...
    'FontSize', 10);

Print = input('Mau di print atau ngga? (Y/N): ', 's');

if 'y' == lower(Print)
    figureHandle = gcf; % Get current figure handle
    set(figureHandle, 'Units', 'normalized', 'OuterPosition', [0, 0, 1, 1]); % Fullscreen

    % Ekspor figure ke file PDF
    try
        exportFileName = input("Nama file: ",'s');
        exportFileName = [exportFileName '.pdf'];
        exportgraphics(figureHandle, exportFileName, 'ContentType', 'vector');
        fprintf('Figure telah diekspor ke: %s\n', exportFileName);
    catch ME
        warning('Gagal mengekspor grafik ke file PDF: %s', ME.message);
    end
else
    fprintf("Nothing happen!\n")
end
clc;
clear;
close all;

% Load the required .mat files for EKF, UKF, and PF results
load('ekf_pfVariations.mat'); % Load EKF results
load('ukf_pfVariations.mat'); % Load UKF results
load('pf100_pfVariations.mat');  % Load PF results
load('pf500_pfVariations.mat');  % Load PF results
load('pf1000_pfVariations.mat');  % Load PF results

%% Input Data
logMSEStructEKF = [];
logMSEStructUKF = [];
logMSEStructPF100  = [];
logMSEStructPF500  = [];
logMSEStructPF1000  = [];

for i = 1:4
    logMSEStructEKF = [logMSEStructEKF log10(ekf_pfVariations.MSEStruct.x(:,i))];
    logMSEStructUKF = [logMSEStructUKF log10(ukf_pfVariations.MSEStruct.x(:,i))];
    logMSEStructPF100  = [logMSEStructPF100  log10(pf100_pfVariations.MSEStruct.x(:,i))];
    logMSEStructPF500  = [logMSEStructPF500  log10(pf500_pfVariations.MSEStruct.x(:,i))];
    logMSEStructPF1000  = [logMSEStructPF1000  log10(pf1000_pfVariations.MSEStruct.x(:,i))];
end

%% Plotting the data
titles = {'x1: Rotor Angle Estimation ($\delta$)',...
          'x2: Rotor Speed Estimation ($\Delta\omega$)', ...
          'x3: Transient voltage -q axis Estimation ($e''_{q}$)',...
          'x4: Transient voltage -d axis Estimation ($e`_{d}$)', ...
          'y: Output Power (Pt)'};
fig = figure('Position', [100, 100, 1200, 600]); % High-quality figure size
t = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact'); % Compact spacing

for j = 1:4
    nexttile; % Move to the next subplot
    hold on;
    grid on
    grid minor;

    % Plot data
    plotEKF = plot(ekf_pfVariations.Deviations, logMSEStructEKF(:,j), ...
        'Color','red','Marker','o',...
        'DisplayName', 'EKF');
    plotUKF = plot(ukf_pfVariations.Deviations, logMSEStructUKF(:,j), ...
        'Color','blue','Marker','square',...
        'DisplayName','UKF');
    plotPF100 = plot(pf100_pfVariations.Deviations, logMSEStructPF100(:,j), ...
        'Color','green','Marker','v',...
        'DisplayName','PF 100');
    plotPF500 = plot(pf500_pfVariations.Deviations, logMSEStructPF500(:,j), ...
        'Color','#EDB120','Marker','hexagram',...
        'DisplayName','PF 500');
    plotPF1000 = plot(pf1000_pfVariations.Deviations, logMSEStructPF1000(:,j), ...
        'Color','#A2142F','Marker','<',...
        'DisplayName','PF 1000');
    title(titles{j}, 'Interpreter', 'latex');
    if j == 2
        plotLegend = [plotEKF, plotUKF, plotPF100, plotPF500, plotPF1000];
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
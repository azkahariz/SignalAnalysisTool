clear;
clc;
close all

load("ekf_wonoise.mat")
%%
plotOnlyThreshold = input("Masukan plot : ");
x = ekf_wonoise.SimulationData.x.values{1};

% Cek konvergensi jika diperlukan
threshold = findConvergenceThreshold(ekf_wonoise,1e-5,1e-6);

if plotOnlyThreshold
    [isConvergentX, isConvergentY] = ekf_wonoise.checkErrorConvergence(threshold);
    indicesToPlot = find(isConvergentX & isConvergentY);
else
    indicesToPlot = 1:ekf_wonoise.MonteCarloIterations;
end

numTimeSteps = size(ekf_wonoise.SimulationData.x_hat.values{1}, 1); % 12001
numStateVars = size(ekf_wonoise.SimulationData.x.values{1}, 2); % 4
numMonteCarlo = length(indicesToPlot); % 18
data = zeros(numTimeSteps, numStateVars, numMonteCarlo); % Preallocation

selectedCells = ekf_wonoise.SimulationData.x_hat.values(indicesToPlot);
data = cat(3, selectedCells{:});

meanEkf = mean(data,3);
stdEkf = std(data,0,3);

time = ekf_wonoise.SimulationData.x.time{1};

%% Plotting
figure;
for i = 1:size(meanEkf, 2)
    subplot(2, 2, i)
    
    % Plot true state
    plotTrue = plot(time, x(:, i), 'Color', 'blue', 'LineStyle', '-', 'LineWidth', 3, ...
         'DisplayName', '$x_{true}$');
    hold on
    
    % Plot estimated state
    plotEst = plot(time, meanEkf(:, i), 'Color', 'red', 'LineStyle', '--', 'LineWidth', 2, ...
         'DisplayName', '$x_{est}$');
    
    % Enable grid
    grid on
    grid minor
    
    % Show legend and set LaTeX interpreter
    legend([plotTrue, plotEst], 'interpreter', 'latex'); % Set the interpreter to 'latex'
end

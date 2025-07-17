clc;

% Retrieve detailed information about variables in the workspace
variableInfo = whos; % Get details of all variables currently in the workspace
variableNames = {variableInfo.name}; % Extract the names of the variables into a cell array

% Display the names of the variables in the workspace
disp('Variables in the loaded .mat files:');
disp(variableNames); % Print the variable names for verification

% Dynamically assign the loaded variables to specific names (ekf, ukf, pf)
% Assuming the variables are in the same order as the files loaded
variableNames{contains(variableNames, "ekf", 'IgnoreCase', true)}
ekf = eval(variableNames{contains(variableNames, "ekf", 'IgnoreCase', true)}); % Assign the first variable to 'ekf'
variableNames{contains(variableNames, "ukf", 'IgnoreCase', true)}
ukf = eval(variableNames{contains(variableNames, "ukf", 'IgnoreCase', true)}); % Assign the second variable to 'ukf'
variableNames{contains(variableNames, "pf", 'IgnoreCase', true)}
pf  = eval(variableNames{contains(variableNames, "pf", 'IgnoreCase', true)}); % Assign the third variable to 'pf'

%%
close all;
ekf.plotMonteCarloSimulationResults('Threshold', true);
ukf.plotMonteCarloSimulationResults('Threshold', true);
pf.plotMonteCarloSimulationResults('Threshold', true);

%%
ekf.plotErrorsBelowThreshold(0.01);
ukf.plotErrorsBelowThreshold(0.01);
pf.plotErrorsBelowThreshold(0.01);

%%
[ekfX, ekfY] = ekf.checkErrorConvergence(ekf.Threshold);
ekfXTotal = sum(ekfX)
ekfYTotal = sum(ekfY)

[ukfX, ukfY] = ukf.checkErrorConvergence(ukf.Threshold);
ukfXTotal = sum(ukfX)
ukfYTotal = sum(ukfY)

[pfX, pfY] = pf.checkErrorConvergence(pf.Threshold);
pXTotal = sum(pfX)
pYTotal = sum(pfY)
clc;
clear;
close all;

% Load the required .mat files for EKF, UKF, and PF results
load('ResultMatFile\meas0.1\ekf_with_measurement_noise0.1.mat'); % Load EKF results
variableInfo = whos; % Get details of all variables currently in the workspace
variableNames = {variableInfo.name}; % Extract the names of the variables into a cell array
load('ResultMatFile\meas0.055\ekf_with_measurement_noise0.055.mat'); % Load EKF results
load('ResultMatFile\meas0.01\ekf_with_measurement_noise0.01.mat'); % Load EKF results
% load('ResultMatFile\meas0.1\ukf_with_measurement_noise0.1.mat'); % Load UKF results
% load('ResultMatFile\meas0.1\pf_with_measurement_noise0.1.mat');  % Load PF results
clc;
clear;
close all;

%% Load the required .mat files for EKF, UKF, and PF results
ekf = LoadMatFileOneVarible('ResultMatFile\process0.001\ekf_with_process_noise0.001.mat');
ukf = LoadMatFileOneVarible('ResultMatFile\process0.001\ukf_with_process_noise0.001.mat');
pf  = LoadMatFileOneVarible('ResultMatFile\process0.001\pf_with_process_noise0.001.mat');

%%
% Mengambil error dari Monte Carlo ke-1
error = ekf.ErrorStruct.x{1};
errorLog = log(error);
time = ekf.SimulationData.x.time{1};

for i = 1:4
    figure(i)
    plot(time, error(:,i))
    grid  on
    grid minor
end

function MatFile = LoadMatFileOneVarible(name)
% Asumsi
loadFile = load(name);
fieldNames = fieldnames(loadFile);
MatFile = loadFile.(fieldNames{1});
end

%% Cek MSE nya untuk monte Carlo ke-1
MSE = ekf.MSEStruct;
% ========================================================================
% File       : main.m
% Deskripsi  : Simulasi Monte Carlo untuk Estimasi Status Dinamis Generator
%              Sinkron menggunakan metode EKF, UKF, dan PF.
%              Termasuk proses pemilihan metode, deviasi inisialisasi,
%              perulangan simulasi, analisis MSE/MAPE/konvergensi, dan
%              visualisasi hasil simulasi.
% Dibuat Oleh: Azka Hariz
% ========================================================================
clear; close all; clc;

% Aktifkan bunyi sistem
beep on;

% -------------------------------
% Setup pencatatan hasil (log)
% -------------------------------
simulationFileName = 'sim_pf_experiment';         % Nama file log (tanpa ekstensi, sesuaikan nama berdasarkan simulasi
simulationFileName = [simulationFileName '.txt']; % Tambahkan ekstensi

% Hapus file lama jika sudah ada
if exist(simulationFileName, 'file')
    delete(simulationFileName);
end

diary(simulationFileName);                        % Mulai catat semua output ke file
diary on;

% -------------------------------
% Inisialisasi Parameter Sistem
% -------------------------------
parameter_ekf_ukf_pf;                   % Panggil skrip parameter (init_state, Q, R, dsb)

% -------------------------------
% Pengaturan Simulasi
% -------------------------------
iteration = 100;                        % Jumlah iterasi Monte Carlo
x_range = [-0.1, 0.1];                  % Range deviasi untuk inisialisasi
num_of_state = 4;
num_of_measurement = 1;
threshold = 0.1;                        % Threshold konvergensi MSE/MAPE

% -------------------------------
% Input Metode Filter dari Pengguna
% -------------------------------
disp('Select filter methods:');
disp('1. Extended Kalman Filter (EKF)');
disp('2. Unscented Kalman Filter (UKF)');
disp('3. Particle Filter (PF)');
selectedFiltersInput = input('Enter your choice (e.g., "1", "2 3", "1 and 2"): ', 's');

% Proses input string untuk identifikasi filter yang dipilih
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

if isempty(selectedFilters)
    error('No valid filter selection. Please enter valid choices.');
end

% -------------------------------
% Pembuatan Deviasi Inisialisasi
% -------------------------------
segment_size = diff(x_range) / iteration;
deviation = x_range(1) + (0:(iteration-1)) * segment_size + rand(1, iteration) * segment_size;

% ========================================================================
% ITERASI SIMULASI MONTE CARLO UNTUK MASING-MASING FILTER
% ========================================================================
for filterIndex = selectedFilters
    % Konfigurasi berdasarkan pilihan filter
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

    % Load dan set waktu simulasi pada model Simulink
    load_system(simModel)
    set_param(simModel, 'StartTime', '0', 'StopTime', '6');

    % Inisialisasi struktur hasil simulasi
    out = repmat(struct('itr', []), iteration, 1);

    % Jalankan simulasi untuk filter yang dipilih
    fprintf("\nRunning simulations using %s \n", filterMethod);

    tic;
    for i = 1:iteration
        fprintf("Iteration %i of %i\n", i, iteration);
        fprintf('Deviation: %f\n', deviation(i));

        % Modifikasi kondisi awal dengan deviasi
        eval(['init_' filterMethod ' = init_state + deviation(i);']);

        % Jalankan simulasi
        out(i).itr = sim(simModel);
    end
    toc;

    % -------------------------------
    % Analisis Hasil Simulasi
    % -------------------------------
    analysisTool = SignalAnalysisTool('Data', out);
    analysisTool.Deviations = deviation;

    % Hitung MSE, MAPE, statistik
    analysisTool.calculateMeanSquaredErrors();
    analysisTool.calculateMeanAbsolutePercentageError();
    analysisTool.calculateStatOfMonteCarlo();

    % Cek konvergensi berdasarkan threshold MSE/MAPE
    [isConvergentX, isConvergentY] = analysisTool.checkErrorConvergence(threshold);
    fprintf("Number of convergent states for %s: %i\n", filterMethod, sum(isConvergentX));
    fprintf("Number of convergent measurements for %s: %i\n", filterMethod, sum(isConvergentY));
    fprintf("Number of convergent states and measurements for %s: %i\n", ...
    filterMethod, sum(isConvergentX & isConvergentY));

    % -------------------------------
    % Visualisasi dan Plot Hasil
    % -------------------------------
    analysisTool.plotErrorsBelowThreshold(threshold);
    analysisTool.plotDeviationVersusError(deviation);
    analysisTool.plotMonteCarloSimulationResults('Threshold', false);
    analysisTool.plotMonteCarloSimulationResults('Threshold', true);
    analysisTool.plotMeanSquaredErrorTimesteps('Threshold', false);
    analysisTool.plotMeanSquaredErrorTimesteps('Threshold', true);
    
    % -------------------------------
    % Simpan Hasil Analisis
    % -------------------------------
    name_var = extractBefore(simulationFileName, '.txt');
    varName = [filterMethod '_' name_var];
    eval([varName ' = analysisTool;']);
    save([varName '.mat'], varName);
    eval(['clear ' varName ';']);
end

% -------------------------------
% Selesai: Tutup file log
% -------------------------------
diary off;
disp(['Simulation results have been saved to "', simulationFileName, '"']);

% Fungsi bantu (tidak dipakai di skrip utama)
function name = getVarName(~)
    name = inputname(1);
end
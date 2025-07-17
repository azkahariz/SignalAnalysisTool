function plot_figure(fileName, exportFileName, export)
% PLOT_FIGURE Memuat file data, menampilkan variabel, memplot hasil, dan
% mengekspor grafik ke PDF.
%
%   plot_figure(fileName, exportFileName) memuat file MATLAB (.mat),
%   menampilkan nama-nama variabel di dalam file beserta nilainya, memplot
%   hasil simulasi Monte Carlo, dan menyimpan grafik ke file PDF.
%
%   Input:
%       fileName (string)       : Nama file .mat yang akan dimuat.
%       exportFileName (string) : Nama file output PDF untuk menyimpan
%       grafik.
%
%   Contoh penggunaan:
%       files = {...
%           {'analysisResults_ekf_gaussian.mat', 'MC_grafik_EKF.pdf'}, ...
%           {'analysisResults_ukf_gaussian.mat', 'MC_grafik_UKF.pdf'}, ...
%           {'analysisResults_pf_gaussian.mat', 'MC_grafik_PF.pdf'}};
%
%       for i = 1:length(files)
%           inputFile = files{i}{1}; outputFile = files{i}{2};
%           fprintf('\n--- Memproses file: %s --x`-\n', inputFile);
%           plot_figure(inputFile, outputFile);
%       end
%
%   Fungsi ini mengasumsikan bahwa data yang dimuat memiliki metode
%   `plotMonteCarloSimulationResults` untuk plotting dengan opsi
%   'Threshold'.
%
%   Penulis: Azka Hariz Sartono 
%   Tanggal: 21 November 2024

    % Memuat file ke dalam struct
    out = load(fileName);

    % Mendapatkan nama semua variabel (field) dalam struct
    fields = fieldnames(out);

    % Menampilkan nama dan nilai setiap variabel
    fprintf('\n--- Daftar Variabel di File: %s ---\n', fileName);
    for i = 1:length(fields)
        varName = fields{i};
        varValue = out.(varName); % Mengakses nilai variabel
        fprintf('Variabel: %s\n', varName);
        disp(varValue); % Menampilkan nilai variabel
    end

    % Memastikan salah satu field memiliki metode
    % plotMonteCarloSimulationResults
    for i = 1:length(fields)
        varValue = out.(fields{i});
        if ismethod(varValue, 'plotMonteCarloSimulationResults')
            % Memplot hasil simulasi Monte Carlo
            fprintf('Memplot hasil simulasi Monte Carlo untuk variabel: %s\n', fields{i});
            varValue.plotMonteCarloSimulationResults('Threshold', false);
            break;
        end
    end

    % Jika tidak ada field dengan metode tersebut
    if ~exist('varValue', 'var') || ~ismethod(varValue, 'plotMonteCarloSimulationResults')
        error('Tidak ada variabel dengan metode plotMonteCarloSimulationResults dalam file.');
    end

    % Set figure yang aktif menjadi fullscreen
    figureHandle = gcf; % Get current figure handle
    set(figureHandle, 'Units', 'normalized', 'OuterPosition', [0, 0, 1, 1]); % Fullscreen

    % Ekspor figure ke file PDF
    if export
        try
            exportgraphics(figureHandle, exportFileName, 'ContentType', 'vector');
            fprintf('Figure telah diekspor ke: %s\n', exportFileName);
        catch ME
            warning('Gagal mengekspor grafik ke file PDF: %s', ME.message);
        end
    end
end

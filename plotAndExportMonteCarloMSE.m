function plotAndExportMonteCarloMSE(name_var)
    % plotAndExportMonteCarloMSE
    %
    % Deskripsi:
    %   Fungsi ini memvisualisasikan hasil analisis Monte Carlo dari metode 
    %   Extended Kalman Filter (EKF), Unscented Kalman Filter (UKF), dan 
    %   Particle Filter (PF) berdasarkan Mean Squared Error (MSE) terhadap 
    %   initial value simulasi Monte Carlo. Grafik yang dihasilkan menampilkan 
    %   sumbu x sebagai initial value simulasi Monte Carlo, dan sumbu y 
    %   menampilkan nilai log10 dari MSE.
    %
    % Parameter:
    %   name_var (string): Nama skenario simulasi Monte Carlo yang telah 
    %   dihasilkan sebelumnya. Fungsi ini akan memuat file .mat terkait 
    %   skenario tersebut, termasuk file untuk metode EKF, UKF, dan PF, untuk 
    %   dibandingkan.
    %
    %   Contoh penggunaan `name_var`: 
    %       Jika skenario simulasi memiliki nama 'ts100' untuk simulasi dengan 
    %       time sampling 100 Hz, fungsi akan memuat file-file yang relevan 
    %       (misalnya `ekf_ts100.mat`, `ukf_ts100.mat`, `pf_ts100.mat`) untuk 
    %       melakukan perbandingan.
    %
    % Output:
    %   Fungsi ini menghasilkan:
    %       - Grafik dengan sumbu x: initial value simulasi Monte Carlo, dan
    %         sumbu y: nilai log10(MSE) dari hasil metode EKF, UKF, dan PF.
    %       - File grafik diekspor ke format gambar .pdf dengan 
    %         nama file yang relevan berdasarkan skenario simulasi.
    %
    % Catatan:
    %   - Pastikan file .mat dengan nama sesuai skenario tersedia di direktori 
    %     kerja. File .mat tersebut harus berisi data MSE dari metode EKF, UKF, 
    %     dan PF untuk dapat diproses oleh fungsi ini.
    %
    % Contoh Pemanggilan:
    %   plotAndExportMonteCarloMSE('ts100');
    %   Fungsi akan memuat file `ekf_ts100.mat`, `ukf_ts100.mat`, dan 
    %   `pf_ts100.mat`, kemudian membandingkan nilai MSE dari ketiga metode 
    %   tersebut dalam grafik logaritmik.
    %
    % Penulis:
    %   Azka Hariz Sartono
    %   Tanggal Terakhir Diperbarui: 16-12-2024

    clc;
    
    ekf = load(['ekf_' name_var '.mat']);
    ukf = load(['ukf_' name_var '.mat']);
    pf  = load(['pf_' name_var '.mat']);

    ekf = ekf.(['ekf_' name_var]);
    ukf = ukf.(['ukf_' name_var]);
    pf  = pf.(['pf_' name_var]);

    %% Parameter Sistem Tetap
    Tm = 0.8;       % Mechanical Torque Input
    Efd = 2.29;     % Steady-state internal voltage of armature
    Vt = 1.02;      % Terminal Voltage

    %% Hitung Steady-State
    init_state = calculateSteadyState(Tm, Efd, Vt);

    %% Validasi Input
    validateInputArguments(ekf, ukf, pf);

    %% Perbarui deviasi awal dengan menambahkan steady-state
    % obj.Deviations yang sebelumnya adalah scalar, menjadi vector
    % berukuran 4x1
    ekf_Deviations = ekf.Deviations + init_state;
    ukf_Deviations = ukf.Deviations + init_state;
    pf_Deviations  = pf.Deviations  + init_state;

    %% Titles untuk State Space
    titles = {
        'x1: Rotor Angle Estimation ($\delta$)', ...
        'x2: Rotor Speed Estimation ($\Delta\omega$)', ...
        'x3: Transient voltage -q axis Estimation ($e''_{q}$)', ...
        'x4: Transient voltage -d axis Estimation ($e''_{d}$)', ...
        'y: Output Power (Pt)'
    };

    %% Hitung log10 MSE
    numStates = ekf.NumOfState;
    logMSEStructEKF = computeLogMSE(ekf, numStates);
    logMSEStructUKF = computeLogMSE(ukf, numStates);
    logMSEStructPF  = computeLogMSE(pf, numStates);

    %% Plot Data
    fig = initializeFigure();
    t = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    plotLegend = [];

    for j = 1:numStates
        nexttile;
        hold on;
        grid on;
        grid minor;

        % Plot data untuk setiap metode
        plotEKF = plot(ekf_Deviations(j, :), logMSEStructEKF(:, j), 'Color', 'red', 'Marker', 'o', ...
            'DisplayName', 'EKF Method');
        plotUKF = plot(ukf_Deviations(j, :), logMSEStructUKF(:, j), 'Color', 'blue', 'Marker', 'square', ...
            'DisplayName', 'UKF Method');
        plotPF = plot(pf_Deviations(j, :), logMSEStructPF(:, j), 'Color', 'green', 'Marker', 'v', ...
            'DisplayName', 'PF Method');

        % Set judul subplot
        title(titles{j}, 'Interpreter', 'latex', 'FontSize', 16);

        % Simpan legenda untuk subplot kedua
        if j == 2
            plotLegend = [plotEKF, plotUKF, plotPF];
        end

        % Set label sumbu
        xlabel('Initial Value MC', 'Interpreter', 'latex', 'FontSize', 14);
        ylabel('$\log_{10}(MSE)$', 'Interpreter', 'latex', 'FontSize', 14);

        % Set axis limit agar menempati seluruh area plot
        xData = [ekf_Deviations(j, :), ukf_Deviations(j, :), pf_Deviations(j, :)];
        yData = [logMSEStructEKF(:, j)', logMSEStructUKF(:, j)', logMSEStructPF(:, j)'];
        xlim([min(xData) - 0.05 * range(xData), max(xData) + 0.05 * range(xData)]);
        ylim([min(yData) - 0.05 * range(yData), max(yData) + 0.05 * range(yData)]);
    end

    % Tambahkan legenda di luar plot
    legend(plotLegend, 'Location', 'northeastoutside', 'Interpreter', 'latex', ...
        'FontSize', 10);

    %% Ekspor Grafik
    exportFigureIfRequested(fig);
    fprintf('Grafik berhasil diplot.\n');
end

%% Fungsi Bantuan
function init_state = calculateSteadyState(Tm, Efd, Vt)
    % Hitung nilai steady-state menggunakan ODE dan fsolve
    x0 = [0.6; 0; 0; 0]; % Nilai awal state
    odefun = @(x) odeFun(x, [Tm, Efd, Vt]); % Definisi ODE
    options = optimoptions('fsolve', 'Display', 'none'); % Opsi fsolve
    init_state = fsolve(odefun, x0, options);
end

function validateInputArguments(ekf, ukf, pf)
    % Validasi input dan properti object
    if nargin ~= 3
        error('Tiga input argumen (ekf, ukf, pf) diperlukan.');
    end

    requiredFields = {'MSEStruct', 'Deviations'};
    inputObjects = {ekf, ukf, pf};
    for i = 1:numel(inputObjects)
        obj = inputObjects{i};
        for j = 1:numel(requiredFields)
            if ~isprop(obj, requiredFields{j})
                error('Input object harus memiliki properti "%s".', requiredFields{j});
            end
        end
    end
end

function logMSEStruct = computeLogMSE(object, numStates)
    % Hitung log10 MSE untuk setiap state
    logMSEStruct = zeros(size(object.MSEStruct.x, 1), numStates);
    for i = 1:numStates
        logMSEStruct(:, i) = log10(object.MSEStruct.x(:, i));
    end
end

function fig = initializeFigure()
    % Inisialisasi figure dengan ukuran tertentu
    fig = figure('Position', [100, 100, 1200, 600]);
end

function exportFigureIfRequested(fig)
    % Tanya apakah grafik ingin diekspor
    Print = input('Mau di print atau ngga? (Y/N): ', 's');
    if strcmpi(Print, 'y')
        set(fig, 'Units', 'normalized', 'OuterPosition', [0, 0, 1, 1]); % Fullscreen
        try
            exportFileName = input("Nama file: ", 's');
            exportFileName = [exportFileName '.pdf'];
            exportgraphics(fig, exportFileName, 'ContentType', 'vector');
            fprintf('Figure telah diekspor ke: %s\n', exportFileName);
        catch ME
            warning('Gagal mengekspor grafik ke file PDF: %s', ME.message);
        end
    else
        fprintf("Nothing happen!\n");
    end
end

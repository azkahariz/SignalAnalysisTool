% plot_MSE_comparison_sampling_times.m
%
% Deskripsi: Kode ini digunakan untuk membandingkan nilai MSE (Mean Squared
% Error) yang dinyatakan dalam skala log10 untuk tiga metode filtering
% (EKF, UKF, dan PF) pada berbagai tingkat frekuensi sampling (100 Hz, 500
% Hz, dan 1000 Hz).
%
% Setelah data simulasi dari tiga metode ini diload ke workspace (dalam
% bentuk variable struktur ekf_ts100, ekf_ts500, ekf_ts1000, ukf_ts100,
% ukf_ts500, ukf_ts1000, pf_ts100, pf_ts500, pf_ts1000), kode ini akan
% menghitung rata-rata log10(MSE) serta standar deviasi untuk masing-masing
% metode pada tiap frekuensi sampling, lalu mem-plot hasilnya dalam bentuk
% errorbar.
%
% Setiap plot menggambarkan hasil estimasi untuk empat state yang berbeda:
% x1 (Rotor Angle), x2 (Rotor Speed), x3 (Transient voltage q-axis), x4
% (Transient voltage d-axis).
%
% Hasil plot kemudian diekspor menjadi file PDF vektor dengan nama yang
% ditentukan oleh pengguna.

% Pastikan data berikut telah diload sebelum menjalankan script ini:
% ekf_ts100, ekf_ts500, ekf_ts1000 ukf_ts100, ukf_ts500, ukf_ts1000
% pf_ts100, pf_ts500, pf_ts1000

% Nama sumbu-x untuk kategori scenario (sampling time)
scenario_names = 'Uniform Noise';

% Kategori yang akan digunakan pada sumbu-x
categories = {'Case 1', 'Case 2', 'Case 3'};
x_categories = 1:length(categories); % posisinya adalah 1, 2, 3 untuk 100Hz, 500Hz, 1000Hz

% Kumpulan Data untuk masing-masing metode filter dalam bentuk cell
ekf = {ekf_uniformnoise03, ekf_uniformnoise02, ekf_uniformnoise01};
ukf = {ukf_uniformnoise03, ukf_uniformnoise02, ukf_uniformnoise01};
pf  = {pf_uniformnoise03,  pf_uniformnoise02,  pf_uniformnoise01};

% Menghitung mean log10(MSE) untuk tiap skenario dengan menggunakan fungsi
% calculateMean serta standard deviasi dengan calculateStd. cellfun
% memudahkan kita untuk mengaplikasikan fungsi tersebut ke setiap elemen
% cell array tanpa loop manual.
logMSE_EKF = cell2mat(cellfun(@(x) calculateMean(x), ekf, 'UniformOutput', false));
logMSE_UKF = cell2mat(cellfun(@(x) calculateMean(x), ukf, 'UniformOutput', false));
logMSE_PF  = cell2mat(cellfun(@(x) calculateMean(x), pf,  'UniformOutput', false));

error_EKF = cell2mat(cellfun(@(x) calculateStd(x), ekf, 'UniformOutput', false));
error_UKF = cell2mat(cellfun(@(x) calculateStd(x), ukf, 'UniformOutput', false));
error_PF  = cell2mat(cellfun(@(x) calculateStd(x), pf,  'UniformOutput', false));

% Offset untuk memisahkan posisi plot masing-masing metode pada sumbu-x
offset = 0.1; 
x_EKF = x_categories - offset; % Posisi EKF sedikit bergeser ke kiri
x_UKF = x_categories;          % Posisi UKF tetap di tengah
x_PF  = x_categories + offset; % Posisi PF sedikit bergeser ke kanan

% Judul untuk masing-masing state yang akan diplot
titles = {
    'x1: Rotor Angle Estimation ($\delta$)', ...
    'x2: Rotor Speed Estimation ($\Delta\omega$)', ...
    'x3: Transient voltage -q axis Estimation ($e''_{q}$)', ...
    'x4: Transient voltage -d axis Estimation ($e''_{d}$)'
};

% Meminta input dari user untuk nama file ekspor (PDF)
exportFileName = input("Nama file: ", 's');

% Parameter plotting
markerSize = 10;
lineWidth = 1.2;
fontSize = 12;

% Mengatur ukuran figure (dalam satuan inci) agar cocok dengan format
% two-column paper
width = 7;   % Lebar (sekitar 18 cm)
height = 3.5; % Tinggi (sekitar 9 cm)

% Looping untuk setiap state agar dibuatkan plotnya
for i=1:length(titles)
    % Membuat figure baru dengan ukuran yang telah ditentukan
    fig = figure('Units', 'inches', 'Position', [1, 1, width, height]);
    hold on;

    % Plot metode EKF
    errorbar(x_EKF, logMSE_EKF(i,:), error_EKF(i,:), 'o', ...
        'MarkerSize', markerSize, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r', ...
        'Color', 'k', 'LineWidth', lineWidth);

    % Plot metode UKF
    errorbar(x_UKF, logMSE_UKF(i,:), error_UKF(i,:), 'o', ...
        'MarkerSize', markerSize, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g', ...
        'Color', 'k', 'LineWidth', lineWidth);

    % Plot metode PF
    errorbar(x_PF, logMSE_PF(i,:), error_PF(i,:), 'o', ...
        'MarkerSize', markerSize, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'b', ...
        'Color', 'k', 'LineWidth', lineWidth);

    % Pemberian label sumbu dan judul plot
    xlabel(scenario_names,'FontSize',fontSize,'Interpreter','latex');
    ylabel('$\log_{10}$MSE','FontSize',fontSize, 'Interpreter','latex');
    title(titles{i},'Interpreter','latex','FontSize',fontSize);

    % Mengatur sumbu-x agar menampilkan kategori dengan benar
    xticks(x_categories);
    xticklabels(categories);
    ax = gca; 
    ax.XAxis.FontSize = fontSize;
    ax.XAxis.TickLabelInterpreter = 'latex';

    % Menambahkan legenda untuk ketiga metode
    legend({'EKF method', 'UKF method', 'PF method'}, ...
        'Location', 'northeastoutside', 'FontSize', 10, ...
        'Interpreter', 'latex');

    % Membuat batas sumbu x sedikit lebih lebar dari data
    xlim([min(x_EKF) - 0.5, max(x_PF) + 0.5]);

    % Menghitung batas sumbu y agar dapat menampung nilai errorbar secara
    % penuh
    all_logMSE = [logMSE_EKF(i,:), logMSE_UKF(i,:), logMSE_PF(i,:)]; 
    all_error = [error_EKF(i,:), error_UKF(i,:), error_PF(i,:)]; 

    % Batas atas dan bawah Y berdasarkan data plot
    y_max = max(all_logMSE + all_error);
    y_min = min(all_logMSE - all_error);
    ylim([y_min - 0.1, y_max + 0.1]);

    % Mengaktifkan grid untuk memudahkan pembacaan plot
    grid on; 
    grid minor;
    hold off;

    % Mengatur ukuran untuk penyimpanan dalam format PDF
    set(fig, 'PaperUnits', 'inches');
    set(fig, 'PaperPosition', [0, 0, width, height]);

    % Mengekspor gambar ke dalam format PDF vektor berkualitas tinggi
    exportgraphics(fig, [exportFileName '_state0' num2str(i) '.pdf'], 'ContentType', 'vector');
end

% Memberikan notifikasi bahwa proses ekspor selesai
fprintf('Figure telah diekspor ke: %s\n', exportFileName);


%% Fungsi Pendukung

% Fungsi calculateMean: Input: object (struktur data berisi x_true dan
% xMean) Tujuan: Menghitung mean dari nilai log10(error) pada setiap state
% Output: meanObject (rata-rata log10(MSE) dalam vektor kolom)
function meanObject = calculateMean(object)
    x_true = object.SimulationData.x.values{1};
    error = abs(object.StatOfMonteCarlo.xMean - x_true);
    log10error = log10(error);
    meanObject = reshape(mean(log10error), [], 1);
end

% Fungsi calculateStd: Input: object (struktur data berisi x_true dan
% xMean) Tujuan: Menghitung standard deviasi dari nilai log10(error) pada
% setiap state Output: stdObject (standar deviasi log10(MSE) dalam vektor
% kolom)
function stdObject = calculateStd(object)
    x_true = object.SimulationData.x.values{1};
    error = abs(object.StatOfMonteCarlo.xMean - x_true);
    log10error = log10(error);
    stdObject = reshape(std(log10error), [], 1);
end

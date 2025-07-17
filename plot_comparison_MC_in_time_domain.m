clc;
clear;
close all;

% Prompt user to select filter methods
disp('Select scenario methods:');
disp('1. Sampling Time');
disp('2. Process Noise');
disp('3. Measurement Noise');
disp('4. Uniform Noise');
selectedScenearioInput = input('Enter your choice: ', 's');

switch selectedScenearioInput
    case '1'
        % Case sampling time
        ekf = {loadMatFile('ekf_ts100'), loadMatFile('ekf_ts500'), loadMatFile('ekf_ts1000')};
        ukf = {loadMatFile('ukf_ts100'), loadMatFile('ukf_ts500'), loadMatFile('ukf_ts1000')};
        pf  = {loadMatFile('pf_ts100'),  loadMatFile('pf_ts500'),  loadMatFile('pf_ts1000')};
        caseLabels = {{'Case' '100 Hz'},{'Case' '500 Hz'},{'Case' '1000 Hz'}};
        saveName = 'plotComparisonMCTimeDomain_ts_state0';
    case '2'
        % Case process noise
        ekf = {loadMatFile('ekf_process0001'), loadMatFile('ekf_process001'), loadMatFile('ekf_process01')};
        ukf = {loadMatFile('ukf_process0001'), loadMatFile('ukf_process001'), loadMatFile('ukf_process01')};
        pf  = {loadMatFile('pf_process0001'),  loadMatFile('pf_process001'),  loadMatFile('pf_process01')};
        caseLabels = {{'Process Noise' 'Case 1'},{'Process Noise' 'Case 2'},{'Process Noise' 'Case 3'}};
        saveName = 'plotComparisonMCTimeDomain_processnoise_state0';
    case '3'
        % Case measurement noise
        ekf = {loadMatFile('ekf_meas01'), loadMatFile('ekf_meas055'), loadMatFile('ekf_meas1')};
        ukf = {loadMatFile('ukf_meas01'), loadMatFile('ukf_meas055'), loadMatFile('ukf_meas1')};
        pf  = {loadMatFile('pf_meas01'),  loadMatFile('pf_meas055'),  loadMatFile('pf_meas1')};
        caseLabels = {{'Meas. Noise' 'Case 1'},{'Meas. Noise' 'Case 2'},{'Meas. Noise' 'Case 3'}};
        saveName = 'plotComparisonMCTimeDomain_measnoise_state0';
    case '4'
        % Case uniform noise
        ekf = {loadMatFile('ekf_uniformnoise03'), loadMatFile('ekf_uniformnoise02'), loadMatFile('ekf_uniformnoise01')};
        ukf = {loadMatFile('ukf_uniformnoise03'), loadMatFile('ukf_uniformnoise02'), loadMatFile('ukf_uniformnoise01')};
        pf  = {loadMatFile('pf_uniformnoise03'),  loadMatFile('pf_uniformnoise02'),  loadMatFile('pf_uniformnoise01')};
        caseLabels = {{'Uniform Noise' 'Case 1'},{'Uniform Noise' 'Case 2'},{'Uniform Noise' 'Case 3'}};
        saveName = 'plotComparisonMCTimeDomain_uniformnoise_state0';
end

%% Define titles, labels, and display name for plots
titles = {'x1: Rotor Angle Estimation ($\delta$)', ...
          'x2: Rotor Speed Estimation ($\Delta\omega$)', ...
          'x3: Transient voltage -q axis Estimation ($e''_{q}$)', ...
          'x4: Transient voltage -d axis Estimation ($e`_{d}$)', ...
          'y: Output Power (Pt)'
          };

yLabels = {'$Elec. Rad (pu)$', ...
           '$Elec. Rad/s (pu)$', ...
           '$Voltage (pu)$', ...
           '$Voltage (pu)$', ...
           '$Active Power (pu)$'
          };

dispNames = {'EKF est.', ...
             'UKF est.', ...
             'PF est', ...
             'True State'
            };

x_true = ekf{1}.SimulationData.x.values{1};
x_time = ekf{1}.SimulationData.x.time{1};

width_cm = 18;   % Lebar dalam cm
height_cm = 11;  % Tinggi dalam cm
fontSize = 12;

% for i = 1:4 untuk state dari 1 sampai 4
for i =  1:4
    % Gunakan tiledlayout agar subplot memiliki ukuran yang sama dan pengaturan lebih mudah
    fig = figure('Units', 'centimeters', 'Position', [0, 0, width_cm, height_cm]);
    t = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % Create axes
    axes_position = {[0.13 0.673920940456227 0.596764705882353 0.243867519981841], ...
                     [0.13 0.395053420474387 0.595588235294118 0.243867519981841], ...
                     [0.13 0.116185900492546 0.597941176470588 0.243867519981841]};
    % for j = 1:3 untuk masing-masing skenario yang jumlahnya ada 3
    for j = 1:3
        est_time = ekf{j}.SimulationData.x_hat.time{1};
        ekf_est  = ekf{j}.StatOfMonteCarlo.xMean;
        ukf_est  = ukf{j}.StatOfMonteCarlo.xMean;
        pf_est   = pf{j}.StatOfMonteCarlo.xMean;

        plotAxes = axes('Parent', fig, 'Position', axes_position{j});
        hold(plotAxes, 'on');

        % Plot sinyal asli
        signalTrue = plot(x_time, x_true(:,i), 'b', ...
                    'LineWidth', 6, 'Parent', plotAxes);
        set(signalTrue, 'DisplayName', dispNames{4});

        % Plot EKF
        plotMeanEKF = plot(est_time, ekf_est(:,i), ...
                    'LineWidth', 3, 'Color', [1 0 0 0.9], ...
                    'LineStyle',":", 'Parent', plotAxes);
        set(plotMeanEKF, 'DisplayName', dispNames{1});
        
        % Plot UKF
        plotMeanUKF = plot(est_time, ukf_est(:,i), ...
                    'LineWidth', 3, 'Color', [0 1 0 0.6], ...
                    'LineStyle',":",'Parent', plotAxes);
        set(plotMeanUKF, 'DisplayName', dispNames{2});
        
        % Plot PF
        plotMeanPF  = plot(est_time,  pf_est(:,i), ...
                    'LineWidth', 3, 'Color', [1 1 0 0.3], ...
                    'LineStyle',":",'Parent', plotAxes);
        set(plotMeanPF, 'DisplayName', dispNames{3});

        % Hanya beri ylabel di setiap subplot, berganti untuk tiap state
        % variabel
        ylabel(yLabels{i},'Interpreter','latex','FontSize',fontSize);

        % Pengaturan axis y agar mencakup data dengan margin
        widthAxes = max(x_true(:,i)) - min(x_true(:,i));
        axis([-Inf, inf, min(x_true(:,i)) - widthAxes, max(x_true(:,i)) + widthAxes ])

        % Grid
        grid on;

        % Hanya contoh teks sesuai dengan kasus:
        caseLabel = caseLabels{j};
        
        % Tambahkan teks rotasi 90 derajat di sebelah kiri subplot
        % Koordinat normalized (-0.1, 0.5) berarti sedikit ke kiri dari area plot
        % dan di tengah ketinggian subplot secara vertikal.
        ax = gca;
        text(ax, 1.05, 0.5, caseLabel, 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'Rotation', 270, ...
            'Interpreter', 'latex', 'FontSize', fontSize, 'Clipping', 'on');

        % Legend hanya muncul di subplot pertama (j == 1)
        if j == 1
            legend([signalTrue plotMeanEKF plotMeanUKF plotMeanPF], ...
                'Position', [0.794999841816281 0.70655181363098 0.184215844458228 0.197836542863112], ...
                'Interpreter', 'latex', ...
                'FontSize', fontSize);
        end

        % Hanya subplot terakhir (j == 3) yang diberi xlabel waktu
        if j == 3
            xlabel('Time ($s$)', 'Interpreter', 'latex','FontSize', fontSize);
        else
            % Subplot 1 dan 2: Hilangkan label sumbu x (tick label), 
            % agar tidak menumpuk dan lebih bersih.
            set(gca,'XTickLabel',[]);
        end
        hold(plotAxes, 'off');
    end

    % Beri judul utama untuk keseluruhan figure menggunakan sgtitle
    sgtitle(titles{i},'Interpreter','latex','FontSize', fontSize);
    
    % Ekspor figure ke dalam format PDF vektor berkualitas tinggi
    exportgraphics(fig, [saveName num2str(i) '.pdf'], 'ContentType','vector');
end
%%
function object = loadMatFile(name_var)
    load(name_var)
    object = eval(name_var);
end
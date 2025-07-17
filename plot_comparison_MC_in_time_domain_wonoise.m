clc;
clear;
close all;

ekf = loadMatFile('ekf_wonoise');
ukf = loadMatFile('ukf_wonoise');
pf  = loadMatFile('pf_wonoise');

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

x_true = ekf.SimulationData.x.values{1};
x_time = ekf.SimulationData.x.time{1};

width_cm = 18;   % Lebar dalam cm
height_cm = 14;  % Tinggi dalam cm
fontSize = 12;

fig = figure('Units', 'centimeters', 'Position', [0, 0, width_cm, height_cm]);
t = tiledlayout(4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
for j = 1:4
    nexttile
    est_time = ekf.SimulationData.x_hat.time{1};
    ekf_est  = ekf.StatOfMonteCarlo.xMean;
    ukf_est  = ukf.StatOfMonteCarlo.xMean;
    pf_est   = pf.StatOfMonteCarlo.xMean;

    % Plot data
    signalTrue = plot(x_time, x_true(:,j), 'b', ...
        'LineWidth', 6);
    set(signalTrue, 'DisplayName', dispNames{4});
    hold on

    plotMeanEKF = plot(est_time, ekf_est(:,j), ...
        'LineWidth', 3, 'Color', [1 0 0 0.9], ...
        'LineStyle',":");
    set(plotMeanEKF, 'DisplayName', dispNames{1});

    plotMeanUKF = plot(est_time, ukf_est(:,j), ...
        'LineWidth', 3, 'Color', [0 1 0 0.6], ...
        'LineStyle',":");
    set(plotMeanUKF, 'DisplayName', dispNames{2});

    plotMeanPF  = plot(est_time,  pf_est(:,j), ...
        'LineWidth', 3, 'Color', [1 1 0 0.3], ...
        'LineStyle',":");
    set(plotMeanPF, 'DisplayName', dispNames{3});

    % Hanya beri ylabel di setiap subplot (jika mau sama)
    ylabel(yLabels{j},'Interpreter','latex','FontSize',fontSize);

    % Pengaturan axis y agar mencakup data dengan margin
    widthAxes = max(x_true(:,j)) - min(x_true(:,j));
    axis([-Inf, inf, min(x_true(:,j)) - widthAxes, max(x_true(:,j)) + widthAxes ])

    % Grid
    grid on;

    % Legend hanya muncul di subplot pertama (j == 1)
    if j == 1
        legend([signalTrue plotMeanEKF plotMeanUKF plotMeanPF], ...
            'Interpreter', 'latex', ...
            'FontSize', fontSize, ...
            'Location','northeastoutside');
    end

    % Hanya subplot terakhir (j == 4) yang diberi xlabel waktu
    if j == 4
        xlabel('Time ($s$)', 'Interpreter', 'latex','FontSize', fontSize);
    else
        % Subplot 1 dan 2: Hilangkan label sumbu x (tick label),
        % agar tidak menumpuk dan lebih bersih.
        set(gca,'XTickLabel',[]);
    end
    title(titles{j}, 'Interpreter', 'latex', 'FontSize', fontSize);
    hold off
end

% Ekspor figure ke dalam format PDF vektor berkualitas tinggi
exportgraphics(fig, 'plotComparisonMCTimeDomain_wonoise.pdf' , 'ContentType','vector');

%%
function object = loadMatFile(name_var)
load(name_var)
object = eval(name_var);
end
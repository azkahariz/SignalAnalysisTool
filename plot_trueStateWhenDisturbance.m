clc
clear;
close all;

load('ResultMatFile\wonoise\ekf_without_noise.mat');

data = analysisResults_ekf;

x = data.SimulationData.x.values{1};
y = data.SimulationData.y.values{1};
time = data.SimulationData.x.time{1};

titles = {'x1: Rotor Angle Estimation ($\delta$)',...
          'x2: Rotor Speed Estimation ($\Delta\omega$)', ...
          'x3: Transient voltage -q axis Estimation ($e''_{q}$)',...
          'x4: Transient voltage -d axis Estimation ($e`_{d}$)', ...
          'y: Output Power (Pt)'};
yLabels = {'$Elec. Rad (pu)$', '$Elec. Rad/s (pu)$', '$Voltage (pu)$', '$Voltage (pu)$', '$Active Power (pu)$'};

widthAxes = max(x) - min(x);


figure;
for i=1:4
    subplot(3,2,i)
    plot(time, x(:,i),'LineStyle','-','Color','blue', 'LineWidth',2);
    title(titles(i),'Interpreter','latex');
    xlabel('$Time(s)$', 'Interpreter', 'latex');
    ylabel(yLabels{i},'Interpreter','latex')
    axis([-Inf Inf min(x(:,i)) - widthAxes(i) max(x(:,i)) + widthAxes(i)])
    grid on
    grid minor
end

widthAxes = max(y) - min(y);
% Subplot untuk Measurement (dilebarkan menjadi 2 kolom)
subplot(3, 1, 3); % Baris ke-3, mengambil 1 kolom penuh
plot(time, y,'LineStyle','-','Color','blue', 'LineWidth',2);
title(titles(5),'Interpreter','latex');
xlabel('$Time(s)$', 'Interpreter', 'latex');
ylabel(yLabels{5},'Interpreter','latex')
axis([-Inf Inf min(y) - widthAxes max(y) + widthAxes])
grid on
grid minor



% % Buat figure
% figure;
% 
% % Subplot untuk State 1
% subplot(3, 2, 1); % Baris ke-1, Kolom ke-2, Plot ke-1
% title('State 1');
% axis off; % Hilangkan axis untuk template
% 
% % Subplot untuk State 2
% subplot(3, 2, 2); % Baris ke-1, Kolom ke-2, Plot ke-2
% title('State 2');
% axis off; % Hilangkan axis untuk template
% 
% % Subplot untuk State 3
% subplot(3, 2, 3); % Baris ke-2, Kolom ke-2, Plot ke-3
% title('State 3');
% axis off; % Hilangkan axis untuk template
% 
% % Subplot untuk State 4
% subplot(3, 2, 4); % Baris ke-2, Kolom ke-2, Plot ke-4
% title('State 4');
% axis off; % Hilangkan axis untuk template
% 
% % Subplot untuk Measurement (dilebarkan menjadi 2 kolom)
% subplot(3, 1, 3); % Baris ke-3, mengambil 1 kolom penuh
% title('Measurement');
% axis off; % Hilangkan axis untuk template

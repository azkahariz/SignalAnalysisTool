clc;clear;close all;
%% Generate data dan transformasi sigmoid
n = 3000;           % Jumlah data
mu = [0; 0];        % Rata-rata
sigma = [1 0; 0 1]; % Matriks kovarians

% Generate data dari distribusi normal multivariat
X = mvnrnd(mu, sigma, n);

% Fungsi sigmoid 2D
sigmoid_2d = @(x1, x2) 1./(1 + exp(-(x1 + x2)));

% Hitung nilai sigmoid untuk setiap titik
Y = sigmoid_2d(X(:,1), X(:,2));

%% Figure 1: Plot 3D Transformasi Sigmoid
figure(1);
hold on;

% Grid untuk permukaan sigmoid
x1_grid = linspace(-4, 4, 50);
x2_grid = linspace(-4, 4, 50);
[X1_grid, X2_grid] = meshgrid(x1_grid, x2_grid);
Y_grid = sigmoid_2d(X1_grid, X2_grid);

% Plot permukaan
surf(X1_grid, X2_grid, Y_grid, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
colormap('parula');

% Plot data transformasi
scatter3(X(:,1), X(:,2), Y, 40, Y, 'filled', 'MarkerFaceAlpha', 0.7);

% Tandai rata-rata
plot3(mu(1), mu(2), sigmoid_2d(mu(1), mu(2)), 'xk', 'MarkerSize', 15, 'LineWidth', 2);

xlabel('X1');
ylabel('X2');
zlabel('y = Sigmoid(X1+X2)');
title('Transformasi Sigmoid 2D dari Data Normal Bivariat');
grid on;
axis tight;
view(-30, 30);
colorbar;
hold off;

%% Figure 2: Distribusi Y
figure(2);

% Histogram dengan kernel density estimate
histogram(Y, 30, 'Normalization', 'pdf', 'FaceColor', [0.5 0.7 1], 'EdgeColor', 'none');
hold on;

% Plot kernel density estimate
[f, xi] = ksdensity(Y);
plot(xi, f, 'LineWidth', 2, 'Color', 'b');

% Garis vertikal di mean
y_mean = mean(Y);
xlims = xlim;
plot([y_mean y_mean], [0 xlims(2)], '--r', 'LineWidth', 1.5);

% Anotasi statistik
text(0.7, max(ylim)*0.9, sprintf('Mean = %.3f\nStd = %.3f', mean(Y), std(Y)), ...
    'BackgroundColor', 'white', 'EdgeColor', 'black');

xlabel('Nilai y');
ylabel('Densitas Probabilitas');
title('Distribusi Variabel Hasil Transformasi Sigmoid');
grid on;
legend({'Histogram', 'Kernel Density', 'Rata-rata'}, 'Location', 'northwest');
hold off;
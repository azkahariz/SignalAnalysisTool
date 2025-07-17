% parametrized_sphere_visualization.m
%
% Deskripsi:
%   Program ini memvisualisasikan parameterisasi permukaan bola dalam 3D
%   menggunakan dua parameter sudut: u (polar/elevasi) dan v (azimuth).
%   Program menghasilkan:
%     - Figure 1: Visualisasi permukaan bola lengkap dengan efek pencahayaan
%     - Figure 2–4: Visualisasi komponen X, Y, dan Z terhadap (u,v)
%                  dengan mesh grid yang lebih jarang agar kontur terlihat jelas
%
% Parameter:
%   R   : Jari-jari bola (default = 1)
%   u   : Sudut polar (0 sampai pi)
%   v   : Sudut azimuth (0 sampai 2pi)
%
% Catatan:
%   Komponen dipetakan dari sistem koordinat bola ke kartesian:
%     x = R*sin(u)*cos(v)
%     y = R*sin(u)*sin(v)
%     z = R*cos(u)
%
%   Grid dikurangi kerapatannya dengan 'skip' agar garis mesh lebih jarang
%   dan sumbu ditampilkan dalam notasi simbolik pi.
%
% Ditulis oleh : Azka Hariz Sartono
% Tanggal      : [Tanggal Hari Ini]

clc;clear;close all;
% Parameter bola
R = 1;

% Buat grid parameter
u = linspace(0, pi, 100);
v = linspace(0, 2*pi, 100);
[U, V] = meshgrid(u, v);

% Komponen vektor parametrik
X = R * sin(U) .* cos(V);
Y = R * sin(U) .* sin(V);
Z = R * cos(U);

%% Figure 1 – Bola lengkap (sudah dijelaskan sebelumnya)
figure(1);
s = surf(X, Y, Z);
s.FaceColor = [1, 0.5, 0];
s.EdgeColor = 'k';
s.FaceLighting = 'gouraud';
s.AmbientStrength = 0.3;
s.DiffuseStrength = 0.8;
s.SpecularStrength = 0.9;
s.SpecularExponent = 25;
camlight headlight;
axis equal;
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Figure 1: Parametrized Sphere');

% Downsample untuk mengurangi kerapatan grid
skip = 5; % Ambil setiap 5 titik
uu = u(1:skip:end);
vv = v(1:skip:end);
[UU, VV] = meshgrid(uu, vv);

% Hitung ulang komponen pada grid yang lebih jarang
XX = R * sin(UU) .* cos(VV);
YY = R * sin(UU) .* sin(VV);
ZZ = R * cos(UU);

% Tick marks dan label dalam pi
tick_vals_u = [0, pi/2, pi];
tick_vals_v = [0, pi/2, pi, 3*pi/2, 2*pi];
tick_labels = {'0', '\pi/2', '\pi', '3\pi/2', '2\pi'};

%% Figure 2 – Komponen X
figure(2);
s2 = surf(uu, vv, XX);
xlabel('u');
ylabel('v');
zlabel('x(u,v)');
title('Figure 2: Komponen X (Reduced Mesh)');
colormap('autumn');
s2.EdgeColor = [0.3, 0.3, 0.3];
s2.FaceColor = 'interp';
s2.LineStyle = '-';
grid on;
colorbar;
xticks(tick_vals_u);
yticks(tick_vals_v);
xticklabels(tick_labels(1:3)); % hanya sampai pi
yticklabels(tick_labels);      % sampai 2pi

%% Figure 3 – Komponen Y
figure(3);
s3 = surf(uu, vv, YY);
xlabel('u');
ylabel('v');
zlabel('y(u,v)');
title('Figure 3: Komponen Y (Reduced Mesh)');
colormap('winter');
s3.EdgeColor = [0.3, 0.3, 0.3];
s3.FaceColor = 'interp';
s3.LineStyle = '-';
grid on;
colorbar;
xticks(tick_vals_u);
yticks(tick_vals_v);
xticklabels(tick_labels(1:3));
yticklabels(tick_labels);

%% Figure 4 – Komponen Z
figure(4);
s4 = surf(uu, vv, ZZ);
xlabel('u');
ylabel('v');
zlabel('z(u,v)');
title('Figure 4: Komponen Z (Reduced Mesh)');
colormap('spring');
s4.EdgeColor = [0.3, 0.3, 0.3];
s4.FaceColor = 'interp';
s4.LineStyle = '-';
grid on;
colorbar;
xticks(tick_vals_u);
yticks(tick_vals_v);
xticklabels(tick_labels(1:3));
yticklabels(tick_labels);

close all;
clc;

load('mean_errors_q1.mat');
load('mean_errors_q2.mat');
load('mean_errors_q3.mat');
load('mean_errors_q4.mat');

plot_myerror(mean_errors_q1);
plot_myerror(mean_errors_q2);
plot_myerror(mean_errors_q3);
plot_myerror(mean_errors_q4);


function plot_myerror(mean_errors)
% figure
% semilogy(mean_errors);
% grid on
% grid minor
% legend show

figure
log10_mean_errors = log10(mean_errors);
plot(log10_mean_errors);
title('Log 10')
grid on
grid minor
legend show

% figure
% log_mean_errors = log(mean_errors);
% plot(log_mean_errors);
% title('Log exp')
% grid on
% grid minor
% legend show
end

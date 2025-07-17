%% Parameter Generator Sinkron
d   = 0.05;                 % Damping Factor
J   = 10;                   % Inertia
Tdo = 0.13;                 % d-axis transient open circuit time constant
Tqo = 0.01;                 % q-axis transient open circuit time constant
xd  = 2.06;                 % d-axis reactance
xq  = 1.21;                 % q-axis reactance
xd1 = 0.37;                 % d-axis transient reactance
xq1 = 0.37;                 % q-axis transient reactance
Tm  = 0.8;                  % Mechanical Torque Input
w0  = 377;                  % Nominal synch. speed
Efd = 2.29;                 % Steady state internal voltage of armature
Vt  = 1.02;                 % Terminal Voltage
x0  = [0.6;0;0;0];    	    % Initial values states

%% Step size
% Sampling Generator
% Sampling metode estimasi ada di model Simulink
% Untuk PMU biasanya dapat melakukan sampling sebanyak
% PMU 20 - 50 x per detik
% dt = 5e-02;      % (sampling = 20/detik)
% dt = 2e-02;      % (sampling = 50/detik)
  dt = 5e-04;      % (sampling = 2000/detik)
% dt = 1e-06;      % (sampling = 1jt/detik)
% dt = 5e-06;      % (sampling = 200000/detik)

%% Jumlah frekuensi sampling masing-masing metode
ekf_freq = 100;
ekf_dt = 1/ekf_freq;

ukf_freq = 100;
ukf_dt = 1/ukf_freq;

pf_freq = 100;
pf_dt = 1/pf_freq;

%% ODE Function
odefun = @(x)odeFun(x, [Tm, Efd, Vt]);

%% Measurement function
ym = @measFun;

%% Find inital condition integrator
% Initial Condition integrator
init_state = fsolve(odefun,x0);

%% Determine process noise and measurement noise for system
% Qsys = [0.001^2, 0.0001^2, 0.001^2, 0.001^2];
% Qsys = [1e-1^2 1e-2^2 1e-1^2 1e-1^2];
Qsys = [1e-2^2 1e-2^2 1e-2^2 1e-2^2];
Rsys = 0;
% Rsys = 0.001^2;

noise_min  =-0.1;
noise_max  = 0.1;

%% Determine initial process noise and measurement noise for method
% Initial Covariance matrix untuk EKF
P0_ekf = 10^3*eye(4);
% Initial Covariance matrix untuk UKF
P0_ukf = 0.0001^2*eye(4);
% P0 = eye(4);
% P0 = 10*eye(4);
% P0 = 0.1*eye(4);

% Process Noise Covariance matrix
% Q0=0.08^2*eye(4);
% Q0 = diag([0.01^2, 0.001^2, 0.01^2, 0.01^2]);
% Q0 = 0.01^2 * eye(4); % (Default process noise yang digunakan)
% Q0 = (1/50)^2 * eye(4) %% (Dari hasil uji empiris PF terbaik)
% Q0 = diag(10*[0.000158477870294882   2.448576952494e-06   6.5111462744571e-05   8.9480570589606e-05])
% Q0 = diag([15 * 0.000158477870294882   150 * 2.448576952494e-06   15 * 6.5111462744571e-05   15 * 8.9480570589606e-05]); % (sejauh ini yang bagus)
% Q0 = diag([ 36 * 0.000158477870294882   250 * 2.448576952494e-06 6.5111462744571e-05   8.9480570589606e-05]);   % (ini juga bagus)
Q0 = diag([ 36 * 0.000158477870294882   2.448576952494e-06 6.5111462744571e-05   8.9480570589606e-05]); 
% Q0 = diag([20 * 0.000158477870294882   150 * 2.448576952494e-06 ...
%     20 * 6.5111462744571e-05   20 * 8.9480570589606e-05]);
% Q0 = diag([0.0004   3.4280077334916e-05   0.000911560478423994   0.00125272798825448])
% Q0 = 0.1^2*eye(4);
% Q0 = 0.0001^2*eye(4);




% Measurement Noise Covariance matrix
% R0 = 0.2^2;
% R0 = 0.001^2;
% R0 = 0.1^2;
% R0 = (1/68)^2; % (Dari hasil uji empiris PF terbaik)
% R0 = 0.000153902173274822; % (terakhir bagus pakai ini)
R0 = 0.002;    % atau 0.003
% R0 = 0.01^2;
% R0 = 1;
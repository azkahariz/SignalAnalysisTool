function xdot = odeFun(x,u)
% format u = Tm, Efd, Vt
%% Machine Parameter
d   = 0.05;                 % Damping Factor
J   = 10;                   % Inertia
Tdo = 0.13;                 % d-axis transient open circuit time constant
Tqo = 0.01;                 % q-axis transient open circuit time constant
xd  = 2.06;                 % d-axis reactance
xq  = 1.21;                 % q-axis reactance
xd1 = 0.37;                 % d-axis transient reactance
xq1 = 0.37;                 % q-axis transient reactance
w0  = 377;              % Nominal synch. speed

%% Input
u1 = u(1);                  % Mechanical Torque Input
u2 = u(2);                 % Steady state internal voltage of armature
u3 = u(3);                 % Terminal Voltage

%% System Model
x1 = x(1);
x2 = x(2);
x3 = x(3);
x4 = x(4);


%% Versi xdot 1
xdot1 = w0*x2;
xdot2 = 1/J *(u1 - (u3/xd1 * x3 * sin(x1) + u3^2/2 * ( 1/xq - 1/xd1 ) * sin(2*x1) ) - d*x2);
xdot3 = 1/Tdo * (u2 - x3 - (xd - xd1) * ( ( x3 - u3*cos(x1) ) /xd1) );
xdot4 = 1/Tqo * (-x4 + (xq - xq1) * (u3/xq * sin(x1)));

xdot = [xdot1;
        xdot2;
        xdot3;
        xdot4];

%% Versi xdot 2
% xdot = [w0*x2;
%        (1/J)*(Tm - (Vt/xd1)*x3*sin(x1) + (Vt^2/2)*(1/xq - 1/xd1)*sin(2*x1) - d*x2);
%        (1/Tdo)*(Efd - x3 - (xd - xd1)*((x3 - Vt*cos(x1))/xd1));
%        (1/Tqo)*(-x4 + (xq - xq1)*((Vt*sin(x1))/xq))];

%% Versi xdot 3
% id = (x3 - u3*cos(x1)) /xd1;
% iq = (u3*sin(x1))/xq;
% ed = u3*sin(x1);
% eq = u3*cos(x1);
% Te = ed*id + eq*iq;
% 
% xdot = [w0*x2;
%         (1/J)*(u1 - Te - d*x2);
% 		1/Tdo*(u2 - x3 - (xd-xd1)*id)
% 		1/Tqo*(-x4 + (xq - xq1)*iq )];
% end
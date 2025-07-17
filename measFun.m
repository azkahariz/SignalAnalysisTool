function Te = measFun(x,u)
% format u = Tm, Efd, Vt
%% Machine Parameter
xq  = 1.21;                 % q-axis reactance
xd1 = 0.37;                 % d-axis transient reactance
Vt  = u(3);                 % Terminal Voltage

%% Measurement function
x1 = x(1);
x3 = x(3);

Te = (Vt/xd1)*x3*sin(x1) + ((Vt^2)/2)*(1/xq - 1/xd1)*sin(2*x1);
end
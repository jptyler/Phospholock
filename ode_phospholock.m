function dxdt = ode_phospholock(~,x,p)

% A function that outputs the ODEs of the phospholock model

% Parameters 

a1 = p(1);
b1 = p(2);
a2 = p(3);
b2 = p(4);
a3 = p(5);
b3 = p(6);
k1f = p(7);
k1r = p(8);
k2f = p(9);
k2r = p(10);
k3 = p(11);
A_T = p(12);

% New Parameters

K1 = (k2f+k2r+k3)/(k2r+k3);
K2 = (k1r*k2r+k1r*k3+k2f*k3)/(k1f*(k2r+k3));

% Solve for AR

AR = (K1*A_T+K1*x(3)+K2-sqrt((K1*A_T+K1*x(3)+K2)^2-4*K1^2*A_T*x(3)))/(2*K1^2);

A = A_T-K1*AR;

% ODEs

dx1 = (1.325/24)*(a1*A-b1*x(1));
dx2 = (1.325/24)*(a2*x(1)-b2*x(2));
dx3 = (1.325/24)*(a3*x(2)-b3*x(3));

% ODE System

dxdt = [dx1;dx2;dx3];

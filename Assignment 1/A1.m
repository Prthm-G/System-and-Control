%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solution to Assignment 1 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize

% Clean up workspace
clear all; clc;
SN    = 90300815;


% Parameters
R1 = 190; 
L1 = 10e-3;
C1 = 13E-6;
R2 = 100;
L2 = 10E-3;
C2 = 18E-6;
Vin = 11;
Iin = 50;
%   QUESTION 1  % 
%%%%%%%%%%%%%%%%
% Solve using tf
%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%
num = [5.162E0 5.371E4 2.090E7];
den = [9.292E-3 8.924E0 7.357E4 3.971];
Q1.G = tf(num, den);

figure(1); clf;
impulse(Q1.G);
grid on;

%%%%%%%%%%%%%%%%%
% Solve using zpk
%%%%%%%%%%%%%%%%%
z = roots(num);
p = roots(den);
xfer2 = zpk(z, p, 1);

[y, t] = impulse(xfer2);


figure(2); clf;
h = plot(t, y, 'k');
set(h, 'LineWidth', 3);
grid on;
title('Dotted Red Impulse Response');
xlabel('Time (sec)');
ylabel('OP Voltage (V)');
set(gca, 'FontSize', 14);

%   QUESTION 2  % 
%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%
% Solve
%%%%%%%%%%%%%%%%
s = tf('s');                % Laplace operator

% Impedances
zr1   = R1;
zr2   = R2;

zl1   = (s*L1);
zl2   = s*L2;

zc1   = 1/(s*C1);
zc2   = 1/(s*C2);

zt1 = R2 + s*L2 + 1/(s*C2);
zt2 = 1/(1/(R1+1/(s*C1))+1/(s*L1));
zt3 = R1 +1/(s*C1);

G11 = (1/zt1 + 1/zr2);
G12 = (-1/zr2);
G13 = (-1/zt1);
G22 = 1/zr2 + 1/zt3 + 1/zl1 +1/zl2;
G23 = -1/zl2;
G33 = 1/zl2 + 1/zt1 + 1/zc2;

% Y matrix

Y = [G11 G12 G13; G12 G22 G23; G13 G23 G33];                                     

% I vector
I = [Iin/s 0 0]';

% Solve for voltage vector
V = inv(Y) * I;

% Display the 2 transfer functions that result
% Could have just omitted the ";" from the previous command
V
V2 = V(2, 1);
Iout = V2/zt2;
Q2.G = Iout/(Iin/s);
Iout


% Compute voltages
[v, t] = impulse(V);

%%%%%%%%%%%%%%%%
% Plot
%%%%%%%%%%%%%%%%
% Plot both voltages together
figure(5); clf; hold on; grid on;
plot(t, v(:,1), 'k-', 'Linewidth', 3);
plot(t, v(:,2), 'b:', 'Linewidth', 3);
title('Node Voltages');
ylabel('Voltage (V)');
legend('Node 1', 'Node 2', 'Location', 'East');
set(gca, 'FontSize', 14);

% Re-Plot voltages at node 2
figure(6); clf; hold on; grid on;
plot(t, v(:,2), 'k-', 'Linewidth', 3);
title('Node  2 Voltage');
ylabel('Voltage (V)');
set(gca, 'FontSize', 14);


% Question 3%
a1DSPlot(SN,2);
y1 = 9.224E1*exp(-9.728*t);
y2 = -9.224E1*exp(-9.728*t);

Q3.K = 8.621e1;
Q3.sigma = -8.596E0;

% QUESTION 4 %

num3 = 3820; % from trail and error %
den3 = [1 -2*Q3.sigma 25*((Q3.sigma)^2)];
%G = -9/(s^2+2*Q3.sigma*s+25*(Q3.sigma)^2);
G = tf(num3,den3);
figure(7); clf;
impulse(G);
grid on;

num4 = 120;
den4 =[1 120];
Gn = tf(num4,den4);
Gp = Gn * G;
figure(8); clf;
impulse(Gp);
grid on;

Q4.p = 120;

% Question 5 %

num5 = 120;
den5 =[1 120];
Gn1 = tf(num5,den5);
Gpn = Gn1* G;
figure(9); clf;
impulse(Gpn);
grid on;

Q5.p = 120;
Q5.err = 0;

%%%%%%%%%%%%%%%%%%
% Create PNG files
%%%%%%%%%%%%%%%%%%
if 0                % saves you the trouble of commenting out many lines
end
